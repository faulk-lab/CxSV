"""Patient intake: bring any patient input to an hs1 CRAM made with the panel recipe.

The 1,619-genome background panel was aligned to T2T-CHM13v2 (UCSC hs1 names)
with one recipe, and a patient CxSV is only comparable with the panel when the
patient reads were aligned the same way (roadmap item 16k; docs/INTAKE.md):

    samtools fastq -T MM,ML,MN,qs | minimap2 2.24-r1122 -y -ax map-ont --eqx -k 17 -K <batch>
        | samtools sort | CRAM against hs1, indexed

Inputs accepted: unaligned BAM (uBAM / dorado modBAM), FASTQ (plain or gzip),
or a BAM/CRAM aligned by another aligner or against another reference naming
(only primary records are used, so every read enters exactly once).

This module is plain Python (no Snakemake objects) so it can be unit-tested;
``workflow/intake/Snakefile`` reads the intake sheet with it and
``workflow/scripts/intake_tool.py`` exposes the stream and check helpers that
``workflow/scripts/intake_align.sh`` calls.
"""
from __future__ import annotations

import csv
import gzip
import re
import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, Iterable

# ─── the panel recipe (docs/POSITIVE_CONTROL.md step 1, docs/INTAKE.md) ──────
RECIPE = {
    "minimap2_version": "2.24-r1122",
    "preset": "map-ont",
    "k": 17,
    "w": 10,                       # map-ont default; recorded in the .mmi header
    "required_args": ("-y", "--eqx"),
    "tags": ("MM", "ML", "MN", "qs"),
}
TAGS = RECIPE["tags"]
PLATFORMS = ("ONT",)               # map-ont recipe only; see docs/INTAKE.md
SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
FASTQ_EXT = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
ALIGNMENT_EXT = (".bam", ".cram")
# a SAM optional field: TAG:TYPE:VALUE (SAM spec 1.4)
_SAM_TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]):([AifZHB]):(.*)$")


class IntakeSheetError(ValueError):
    """The intake sheet cannot be turned into a valid list of samples."""


class IntakeCheckError(RuntimeError):
    """An intake output failed one of its integrity checks."""


@dataclass(frozen=True)
class IntakeRecord:
    sample_id: str
    inputs: tuple[str, ...]
    individual_id: str
    platform: str = "ONT"
    input_reference: str = ""      # FASTA needed to decode CRAM inputs (another reference)
    description: str = ""          # free text for the @RG DS field (basecaller, kit, ...)
    formats: tuple[str, ...] = field(default=())

    @property
    def format_summary(self) -> str:
        return ",".join(sorted(set(self.formats)))


def input_format(path: str | Path) -> str:
    """'fastq', 'bam' or 'cram' from the file name (case-insensitive)."""
    name = str(path).lower()
    if name.endswith(FASTQ_EXT):
        return "fastq"
    for ext in ALIGNMENT_EXT:
        if name.endswith(ext):
            return ext[1:]
    raise IntakeSheetError(f"unsupported input type (need .bam, .cram, .fastq/.fq[.gz]): {path}")


def read_intake_sheet(sheet: str | Path, check_files: bool = True) -> list[IntakeRecord]:
    """Parse the tab-separated intake sheet, listing every problem at once.

    Required columns: sample_id, input (one path, or several separated by ','
    for a sample sequenced on several flow cells). Optional: individual_id
    (default sample_id), platform (default ONT; only ONT is accepted),
    input_reference (FASTA for CRAM inputs whose @SQ M5 are not in the bundle's
    hs1 MD5 cache), description (text for the read group DS field).
    """
    with open(sheet, newline="") as fh:
        rows = list(csv.DictReader((l for l in fh if l.strip() and not l.startswith("#")), delimiter="\t"))
    if not rows:
        raise IntakeSheetError(f"intake sheet {sheet} has no rows")
    for col in ("sample_id", "input"):
        if col not in rows[0]:
            raise IntakeSheetError(f"intake sheet {sheet} lacks required column '{col}'")
    problems: list[str] = []
    seen: set[str] = set()
    out: list[IntakeRecord] = []
    for i, row in enumerate(rows, start=2):
        sid = (row.get("sample_id") or "").strip()
        raw = (row.get("input") or "").strip()
        inputs = tuple(p.strip() for p in raw.split(",") if p.strip())
        platform = (row.get("platform") or "ONT").strip() or "ONT"
        ref = (row.get("input_reference") or "").strip()
        if not sid:
            problems.append(f"line {i}: empty sample_id")
            continue
        if not SAMPLE_ID_RE.match(sid):
            problems.append(f"line {i}: sample_id '{sid}' may only contain letters, digits, '.', '_' and '-'")
        if sid in seen:
            problems.append(f"line {i}: duplicate sample_id '{sid}'")
        seen.add(sid)
        if platform.upper() not in PLATFORMS:
            problems.append(f"line {i}: {sid}: platform '{platform}' not supported by the map-ont recipe "
                            f"(accepted: {', '.join(PLATFORMS)})")
        if not inputs:
            problems.append(f"line {i}: {sid}: no input")
        formats = []
        for p in inputs:
            try:
                formats.append(input_format(p))
            except IntakeSheetError as e:
                problems.append(f"line {i}: {sid}: {e}")
                continue
            if check_files and not Path(p).is_file():
                problems.append(f"line {i}: {sid}: input not found: {p}")
        if len(set(inputs)) != len(inputs):
            problems.append(f"line {i}: {sid}: the same input is listed twice")
        if ref and check_files and not Path(ref).is_file():
            problems.append(f"line {i}: {sid}: input_reference not found: {ref}")
        out.append(IntakeRecord(sample_id=sid, inputs=inputs,
                                individual_id=(row.get("individual_id") or "").strip() or sid,
                                platform=platform.upper(), input_reference=ref,
                                description=(row.get("description") or "").strip(),
                                formats=tuple(formats)))
    if problems:
        shown = "\n  ".join(problems[:50])
        more = f"\n  ... and {len(problems) - 50} more" if len(problems) > 50 else ""
        raise IntakeSheetError(f"intake sheet {sheet} has {len(problems)} problem(s):\n  {shown}{more}")
    return out


# ─── FASTQ comments -> SAM tags minimap2 -y can copy ─────────────────────────
def keep_tags(comment: str, tags: Iterable[str] = TAGS) -> str:
    """Tab-joined SAM tags from a FASTQ comment, keeping only `tags`.

    dorado writes tags tab-separated; other tools separate tokens with spaces
    or add non-SAM tokens (runid=..., ch=...). minimap2 -y copies the comment
    verbatim into the SAM record, so anything else would make invalid SAM.
    MM/ML/MN/qs values never contain whitespace, so splitting on whitespace is safe.
    """
    wanted = set(tags)
    kept = []
    for tok in comment.split():
        m = _SAM_TAG_RE.match(tok)
        if m and m.group(1) in wanted:
            kept.append(tok)
    return "\t".join(kept)


def _open_text(path: str) -> IO[str]:
    if path == "-":
        import sys
        return sys.stdin
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def normalise_fastq(paths: Iterable[str], out: IO[str], tags: Iterable[str] = TAGS) -> int:
    """Copy FASTQ records to `out` with comments reduced to the recipe tags; returns
    the number of reads. Raises IntakeCheckError on a truncated or malformed record."""
    n = 0
    for path in paths:
        with _open_text(path) as fh:
            while True:
                head = fh.readline()
                if not head:
                    break
                seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
                if not head.startswith("@") or not plus.startswith("+") or not qual:
                    raise IntakeCheckError(f"{path}: malformed or truncated FASTQ record {n + 1} ({head[:60]!r})")
                seq, qual = seq.rstrip("\n"), qual.rstrip("\n")
                if len(seq) != len(qual):
                    raise IntakeCheckError(f"{path}: record {n + 1}: SEQ and QUAL lengths differ")
                parts = head[1:].rstrip("\n").split(None, 1)
                if not parts:
                    raise IntakeCheckError(f"{path}: record {n + 1} has no read name")
                name, comment = parts[0], (parts[1] if len(parts) > 1 else "")
                tagstr = keep_tags(comment, tags)
                out.write(f"@{name}\t{tagstr}\n" if tagstr else f"@{name}\n")
                out.write(f"{seq}\n+\n{qual}\n")
                n += 1
    return n


# ─── minimap2 index header ───────────────────────────────────────────────────
def mmi_header(path: str | Path) -> dict:
    """w, k, b, n_seq and flag from a minimap2 .mmi index (format magic MMI\\2)."""
    with open(path, "rb") as fh:
        head = fh.read(24)
    if len(head) < 24 or head[:4] != b"MMI\x02":
        raise IntakeCheckError(f"{path} is not a minimap2 index (magic {head[:4]!r})")
    w, k, b, n_seq, flag = struct.unpack("<5I", head[4:24])
    return {"w": w, "k": k, "b": b, "n_seq": n_seq, "flag": flag}


def check_mmi(path: str | Path, n_contigs: int | None = None) -> dict:
    h = mmi_header(path)
    problems = []
    if h["k"] != RECIPE["k"] or h["w"] != RECIPE["w"]:
        problems.append(f"k={h['k']} w={h['w']}, recipe needs k={RECIPE['k']} w={RECIPE['w']}")
    if n_contigs is not None and h["n_seq"] != n_contigs:
        problems.append(f"{h['n_seq']} sequences, reference .fai has {n_contigs}")
    if h["flag"] & 0x1:   # MM_I_HPC
        problems.append("homopolymer-compressed index (map-pb style), recipe is map-ont")
    if problems:
        raise IntakeCheckError(f"{path}: " + "; ".join(problems))
    return h


# ─── output checks ───────────────────────────────────────────────────────────
def parse_flagstat(text: str) -> dict:
    """primary, primary_mapped, secondary, supplementary from `samtools flagstat` (QC-passed column)."""
    want = {"primary": r" primary$", "primary_mapped": r" primary mapped \(",
            "secondary": r" secondary$", "supplementary": r" supplementary$", "total": r" in total \("}
    out = {}
    for line in text.splitlines():
        for key, pat in want.items():
            if key not in out and re.search(pat, line):
                a, b = line.split(" + ", 1)
                out[key] = int(a) + int(b.split()[0])
    missing = [k for k in ("primary", "primary_mapped") if k not in out]
    if missing:
        raise IntakeCheckError(f"flagstat output lacks {missing}")
    out["unmapped"] = out["primary"] - out["primary_mapped"]
    return out


def parse_samtools_fastq_count(stderr_text: str) -> int:
    """Reads written by `samtools fastq` ('[M::bam2fq_mainloop] processed N reads')."""
    m = re.findall(r"processed (\d+) reads", stderr_text)
    if not m:
        raise IntakeCheckError("samtools fastq did not report a read count")
    return int(m[-1])


def check_header(header_text: str, sample_id: str) -> list[str]:
    """Problems with an intake output header: exactly one @RG with SM == sample_id,
    and a minimap2 @PG of the recipe version with the recipe arguments.

    Written for intake outputs, whose @PG carries an explicit '-k 17'. Panel CRAMs
    made before the intake (k from the prebuilt .mmi only) have no '-k 17' token
    and would be reported here although their alignments are the same."""
    problems = []
    rgs = [l for l in header_text.splitlines() if l.startswith("@RG\t")]
    sms = [f.split(":", 1)[1] for l in rgs for f in l.split("\t") if f.startswith("SM:")]
    if len(rgs) != 1:
        problems.append(f"{len(rgs)} @RG lines, expected 1")
    if sms != [sample_id]:
        problems.append(f"@RG SM {sms} != sample_id {sample_id!r}")
    pgs = [dict(f.split(":", 1) for f in l.split("\t")[1:] if ":" in f)
           for l in header_text.splitlines() if l.startswith("@PG\t")]
    mm = [p for p in pgs if p.get("PN") == "minimap2" or p.get("ID", "").startswith("minimap2")]
    if len(mm) != 1:
        problems.append(f"{len(mm)} minimap2 @PG lines, expected 1 (a re-used alignment?)")
    else:
        p = mm[0]
        if p.get("VN") != RECIPE["minimap2_version"]:
            problems.append(f"minimap2 {p.get('VN')} != {RECIPE['minimap2_version']}")
        argv = p.get("CL", "").split()
        has_preset = any(argv[i] in ("-x", "-ax") and i + 1 < len(argv) and argv[i + 1] == RECIPE["preset"]
                         for i in range(len(argv)))
        if not has_preset:
            problems.append(f"minimap2 CL lacks -x {RECIPE['preset']}")
        if "-a" not in argv and "-ax" not in argv:
            problems.append("minimap2 CL lacks -a")
        for a in RECIPE["required_args"]:
            if a not in argv:
                problems.append(f"minimap2 CL lacks {a}")
        if not any(argv[i] == "-k" and i + 1 < len(argv) and argv[i + 1] == str(RECIPE["k"]) for i in range(len(argv))):
            problems.append(f"minimap2 CL lacks -k {RECIPE['k']}")
    return problems


def verify(sample_id: str, reads_in: int, flagstat_text: str, header_text: str,
           mod_records_in: int | None = None, mod_records_out: int | None = None) -> dict:
    """Integrity verdict of one intake output; raises IntakeCheckError listing every failure.

    Read count: every primary input read must come out exactly once, as a
    primary alignment or an unmapped record (input primary == output
    primary mapped + unmapped). Modified-base tags: when the input sample of
    reads carries MM/ML, the output sample must too (minimap2 -y).
    """
    fs = parse_flagstat(flagstat_text)
    problems = check_header(header_text, sample_id)
    if reads_in != fs["primary_mapped"] + fs["unmapped"]:
        problems.append(f"read count: input primary {reads_in} != output primary mapped "
                        f"{fs['primary_mapped']} + unmapped {fs['unmapped']}")
    if reads_in == 0:
        problems.append("no reads in the input")
    if mod_records_in and not mod_records_out:
        problems.append(f"input carries MM/ML tags ({mod_records_in} of the sampled records) "
                        f"but the output sample has none")
    result = {"sample_id": sample_id, "reads_in": reads_in, **{f"out_{k}": v for k, v in fs.items()},
              "mod_records_in_sample": mod_records_in, "mod_records_out_sample": mod_records_out,
              "verdict": "PASS" if not problems else "FAIL", "problems": problems}
    if problems:
        raise IntakeCheckError(f"{sample_id}: " + "; ".join(problems))
    return result


def patient_sheet_rows(records: Iterable[IntakeRecord], cram_of) -> list[dict]:
    """Rows of the patient sample sheet (workflow/lib/cxsv_inputs.py format) for the
    intake outputs; `cram_of(sample_id)` gives the CRAM path."""
    rows = []
    for r in records:
        cram = str(cram_of(r.sample_id))
        rows.append({"sample_id": r.sample_id, "alignment": cram, "index": cram + ".crai",
                     "individual_id": r.individual_id, "platform": r.platform})
    return rows


def write_tsv(rows: list[dict], path: str | Path, columns: list[str]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

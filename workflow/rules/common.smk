# ═══════════════════════════════════════════════════════════════════════════════
# COMMON — sample discovery, alignment inputs, and CRAM reference handling
# Included first by the Snakefile; everything below is plain Python + helpers.
# ═══════════════════════════════════════════════════════════════════════════════
#
# A cohort's samples come from ONE of two places (config/config.yaml -> cohorts):
#
#   Directory mode (default)
#       Every *.bam and *.cram in bamfiles/<cohort>/. Sample ID = the filename
#       minus .bam/.cram. A missing index (.bam.bai / .cram.crai) is created
#       next to the file by rules/indexing.smk.
#
#   Sample-sheet mode   cohorts.<name>.samples: <path/to/sheet.tsv>
#       Tab-separated, header required, columns:
#           sample_id   alignment   [index]   [coverage_x]  ...extra columns ignored
#       Alignments may live anywhere (e.g. a shared, read-only datasets tree).
#       Their indexes must already exist: the pipeline never writes into a
#       sample sheet's data directories. Optional cohorts.<name>.min_coverage
#       keeps only rows whose coverage_x is at or above that value.
#       The background-set sheets are produced by the acquisition code's
#       build_registry.py (registry/sample_sheets/<cohort>.tsv).
#
# CRAM reference handling
#   A CRAM stores reads as differences from the reference, so every tool that
#   reads one must find the exact same reference sequence. htslib looks it up by
#   each @SQ line's M5 (MD5) tag via REF_PATH, else by the UR path baked into the
#   header -- which is only valid on the machine that wrote the CRAM. We build an
#   MD5-keyed cache of config['reference'] once (rule build_ref_cache), point
#   REF_PATH/REF_CACHE at it for every tool (HTS_ENV), and check each CRAM's M5
#   tags against it before calling (rule check_cram_reference). With REF_PATH
#   set to the local cache only, htslib never falls back to a network download.

import csv
import re
from pathlib import Path

ALIGNMENT_INDEX_EXT = {".bam": ".bai", ".cram": ".crai"}


def _default_index(alignment):
    return alignment + ALIGNMENT_INDEX_EXT[Path(alignment).suffix]


def _samples_from_dir(data_dir):
    table = {}
    for ext in ALIGNMENT_INDEX_EXT:
        for p in sorted(Path(data_dir).glob(f"*{ext}")):
            sid = p.name[: -len(ext)]
            if sid in table:
                raise ValueError(
                    f"sample '{sid}' appears as both {table[sid]['alignment']} and {p} "
                    f"in {data_dir}/ -- keep one alignment per sample"
                )
            table[sid] = {"alignment": str(p), "index": _default_index(str(p))}
    return table


def _samples_from_sheet(sheet, min_coverage=None):
    with open(sheet, newline="") as fh:
        rows = list(csv.DictReader((ln.replace("\r", "") for ln in fh), delimiter="\t"))
    if not rows:
        raise ValueError(f"sample sheet {sheet} has no rows")
    for col in ("sample_id", "alignment"):
        if col not in rows[0]:
            raise ValueError(f"sample sheet {sheet} lacks required column '{col}'")
    if min_coverage is not None and "coverage_x" not in rows[0]:
        raise ValueError(f"min_coverage is set but {sheet} has no coverage_x column")

    table, problems, dropped = {}, [], 0
    for i, row in enumerate(rows, start=2):
        sid, aln = row["sample_id"].strip(), row["alignment"].strip()
        if min_coverage is not None:
            try:
                if float(row["coverage_x"]) < float(min_coverage):
                    dropped += 1
                    continue
            except ValueError:
                problems.append(f"line {i}: {sid}: coverage_x '{row['coverage_x']}' is not a number")
                continue
        if Path(aln).suffix not in ALIGNMENT_INDEX_EXT:
            problems.append(f"line {i}: {sid}: alignment must end in .bam or .cram: {aln}")
            continue
        idx = (row.get("index") or "").strip() or _default_index(aln)
        if sid in table:
            problems.append(f"line {i}: duplicate sample_id '{sid}'")
        if not Path(aln).exists():
            problems.append(f"line {i}: {sid}: alignment not found: {aln}")
        if not Path(idx).exists():
            problems.append(f"line {i}: {sid}: index not found: {idx}")
        table[sid] = {"alignment": aln, "index": idx}
    if problems:
        shown = "\n  ".join(problems[:20])
        more = f"\n  ... and {len(problems) - 20} more" if len(problems) > 20 else ""
        raise ValueError(f"sample sheet {sheet} has {len(problems)} problem(s):\n  {shown}{more}")
    if dropped:
        print(f"[cohort={COHORT}] min_coverage={min_coverage}: {dropped} sample(s) below the floor "
              f"left out, {len(table)} kept", file=sys.stderr)
    return table


_cohort_cfg = COHORTS_CFG[COHORT] or {}
SAMPLE_SHEET = _cohort_cfg.get("samples")
if SAMPLE_SHEET:
    SAMPLE_TABLE = _samples_from_sheet(SAMPLE_SHEET, _cohort_cfg.get("min_coverage"))
    INPUT_SOURCE = f"sample sheet {SAMPLE_SHEET}"
else:
    SAMPLE_TABLE = _samples_from_dir(DATA)
    INPUT_SOURCE = f"{DATA}/"

SAMPLES = sorted(SAMPLE_TABLE)
CRAM_SAMPLES = {s for s in SAMPLES if SAMPLE_TABLE[s]["alignment"].endswith(".cram")}

# Pin {sample} to the known IDs so dotted IDs (HG00096.ERR15142589...) can
# never be split ambiguously against other wildcard patterns. (The Snakefile
# stops with a clear error when SAMPLES is empty.)
SAMPLE_REGEX = "|".join(re.escape(s) for s in SAMPLES) or "(?!)"

wildcard_constraints:
    sample=SAMPLE_REGEX,

# ─── Reference cache (shared by every cohort) ────────────────────────────────
REF_CACHE = os.path.abspath(config.get("ref_cache", "resources/ref_cache"))
REF_MD5_MANIFEST = f"{REF_CACHE}/md5_manifest.tsv"
# Prefixed to every shell command that may read a CRAM. '%2s/%2s/%s' is
# htslib's MD5 path template (first two hex chars / next two / remainder).
HTS_ENV = (
    f"export REF_PATH='{REF_CACHE}/%2s/%2s/%s' REF_CACHE='{REF_CACHE}/%2s/%2s/%s'; "
)


# ─── Input functions used by indexing.smk / calling.smk ──────────────────────
def alignment(wildcards):
    return SAMPLE_TABLE[wildcards.sample]["alignment"]


def alignment_index(wildcards):
    return SAMPLE_TABLE[wildcards.sample]["index"]


def caller_inputs(wildcards):
    """aln + idx for every sample; CRAMs also wait for their reference check."""
    s = wildcards.sample
    inputs = {"aln": SAMPLE_TABLE[s]["alignment"], "idx": SAMPLE_TABLE[s]["index"]}
    if s in CRAM_SAMPLES:
        inputs["refcheck"] = f"{RESULTS}/refcheck/{s}.ok"
    return inputs

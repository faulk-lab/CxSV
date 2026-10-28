#!/usr/bin/env python3
"""Rebuild the CxSV background-genome registry from what is actually on disk.

Writes three things under ``<datasets>/registry/``:

``sample_registry.tsv``
    One row per alignment file in the background set (columns and vocabulary in
    ``registry_schema.py``).  Index files (.bai/.crai) are not separate rows; a
    row is only ``verified`` when its index is present next to it.
``individuals.tsv``
    One row per person (normalised NA#####/HG#####), so the same individual
    sequenced in several cohorts is counted ONCE in background frequency.
``sample_sheets/<cohort>.tsv``
    Pipeline input: sample_id, alignment, index (absolute paths), plus
    coverage/floor columns.  Only verified rows are included.  Point a cohort's
    ``samples:`` key in the pipeline's config.yaml at one of these.

Why rebuild instead of append: every input below is already authoritative
(download ledgers, per-sample integrity records, the files themselves), so a
full rebuild is idempotent and cannot drift from the filesystem.  Re-run it
after any acquisition step.

Sources
-------
* manifests/manifest_tierA_A{1,2,3}.tsv + logs/tierA_A*_download.freya.tsv
  (1kGP-LRSC ONT BAMs, 1kGP PacBio HiFi BAMs)
* logs/relocate_ledger.tsv (the held 100 1kGP-LRSC BAMs)
* manifests/manifest_PRJEB89727_odp.tsv + logs/align/integrity/<ACC>.tsv
  (the 1,017 in-house CHM13 CRAMs)
* manifests/manifest_PRJEB89727.tsv (the 2 ENA-only runs, recorded as excluded)
* registry/sources/: 1kGP_LRSC_500_ONT_Metadata.tsv (pore, kit, sex, ancestry),
  ont500_cramino_coverage.tsv, pacbio_chr20_depth.tsv

Usage
-----
    python3 build_registry.py                       # default datasets root
    python3 build_registry.py --datasets /mnt/sargasso1/CxSV/datasets
    python3 build_registry.py --dry-run             # report counts, write nothing
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from registry_schema import COLUMNS, normalise_individual_id  # noqa: E402

COVERAGE_FLOOR_X = 17.0          # Chris, 2026-10-22: floor on raw yield / reported depth
REFERENCE = "T2T-CHM13v2"
LRSC_COHORT = "1kGP_LRSC_ONT_CHM13"
PACBIO_COHORT = "1kGP_PacBio_CHM13"
ONT1019_COHORT = "1kGP_ONT_1019_CHM13"
ENA_ONLY_RUNS = {   # PI decision 2026-10-27: FASTQ-only on ENA, not on the AWS mirror
    "ERR15142724": "HG02013",
    "ERR15142661": "HG00377",
}


# ─── small I/O helpers ───────────────────────────────────────────────────────
def read_tsv(path: Path) -> list[dict]:
    """Read a TSV written on Windows/SMB or Linux alike (strips CR)."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader((ln.replace("\r", "") for ln in fh), delimiter="\t"))
    return rows


def last_status(ledger: Path) -> dict[str, dict]:
    """filename -> last ledger row (ledgers are append-only; resumes add rows)."""
    out: dict[str, dict] = {}
    if ledger.exists():
        for row in read_tsv(ledger):
            out[row["filename"]] = row
    return out


def iso_date(ts: str) -> str:
    return ts[:10] if ts else ""


def fsize(p: Path) -> str:
    try:
        return str(p.stat().st_size)
    except FileNotFoundError:
        return ""


def write_atomic(path: Path, header: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, path)


# ─── metadata parsers ────────────────────────────────────────────────────────
def a2_basecaller(fname: str) -> str:
    """'...LSK110.dorado050.sup.5mCG...' -> 'dorado-0.5.0-sup';
    '...LSK110.guppy657.sup.5mc_5hmC...' -> 'guppy-6.5.7-sup'."""
    m = re.search(r"\.dorado(\d)(\d)(\d)\.(\w+)\.", fname)
    if m:
        return f"dorado-{m.group(1)}.{m.group(2)}.{m.group(3)}-{m.group(4)}"
    m = re.search(r"\.guppy(\d)(\d)(\d)\.(\w+)\.", fname)   # '.guppy657.sup.' -> guppy-6.5.7-sup
    if m:
        return f"guppy-{m.group(1)}.{m.group(2)}.{m.group(3)}-{m.group(4)}"
    m = re.search(r"\.guppy\.(\w+)\.", fname)
    return f"guppy-{m.group(1)}" if m else ""


def load_lrsc_metadata(src: Path) -> dict[str, dict]:
    meta = {}
    for row in read_tsv(src / "1kGP_LRSC_500_ONT_Metadata.tsv"):
        meta[normalise_individual_id(row["NHGRI_ID"])] = row
    return meta


def load_cramino(src: Path) -> dict[str, float]:
    return {normalise_individual_id(r["sample"]): float(r["mean_cov"])
            for r in read_tsv(src / "ont500_cramino_coverage.tsv") if r.get("mean_cov")}


def load_pacbio_depth(src: Path) -> dict[str, float]:
    out = {}
    for r in read_tsv(src / "pacbio_chr20_depth.tsv"):
        if r.get("meandepth_chr20_1Mb"):
            out[Path(r["key"]).name] = float(r["meandepth_chr20_1Mb"])
    return out


# ─── row builders, one per cohort ────────────────────────────────────────────
def base_row(**kw) -> dict:
    row = {c: "" for c in COLUMNS}
    row.update({k: ("" if v is None else str(v)) for k, v in kw.items()})
    return row


def lrsc_rows(ds: Path, meta, cramino) -> list[dict]:
    """Held 100 (relocated) + Tier A1 (R10 PMDV) + Tier A2 (R9 MINIMAP_2026)."""
    rows = []
    cdir = ds / LRSC_COHORT
    # held 100: relocate ledger
    for r in read_tsv(ds / "logs" / "relocate_ledger.tsv"):
        rows.append(("held", r["bam"], "1000g-ont:ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100 (held copy)",
                     r["status"], r["ts"], r["src_bytes"]))
    for tier in ("A1", "A2"):
        led = last_status(ds / "logs" / f"tierA_{tier}_download.freya.tsv")
        for m in read_tsv(ds / "manifests" / f"manifest_tierA_{tier}.tsv"):
            if not m["filename"].endswith(".bam"):
                continue
            L = led.get(m["filename"], {})
            rows.append((tier, m["filename"], m["bucket_key"], L.get("status", "pending"),
                         L.get("ts", ""), m["size_bytes"]))
    out = []
    for tier, fname, key, status, ts, exp_bytes in rows:
        iid = normalise_individual_id(fname)
        md = meta.get(iid, {})
        bam = cdir / fname
        have = bam.exists() and (cdir / (fname + ".bai")).exists()
        size = fsize(bam)
        ok = status == "verified" and have and size == str(exp_bytes)
        cov = cramino.get(iid)
        notes = [f"tier={tier}"]
        if cov is not None:
            notes.append("coverage from cramino on the published hg38 alignment of the same reads")
        out.append(base_row(
            file_id=f"{LRSC_COHORT}/{fname}", individual_id=iid,
            source_sample_id=re.match(r"[A-Z]{2}\d{5}", fname).group(0),
            cohort=LRSC_COHORT, accession=key, platform="ONT",
            chemistry=md.get("ONT_pore", ""), library=md.get("ONT_library", ""),
            basecaller=a2_basecaller(fname) if tier == "A2" else "",
            data_type="aligned_bam", reference=REFERENCE,
            coverage_x=f"{cov:.2f}" if cov is not None else "",
            coverage_source="cramino" if cov is not None else "unknown",
            rel_path=f"{LRSC_COHORT}/{fname}", size_bytes=size, md5="",
            status="verified" if ok else ("downloading" if (cdir / (fname + ".part")).exists() else status or "pending"),
            date_added=iso_date(ts), notes="; ".join(notes)))
    return out


def pacbio_rows(ds: Path, depth) -> list[dict]:
    led = last_status(ds / "logs" / "tierA_A3_download.freya.tsv")
    cdir = ds / PACBIO_COHORT
    out = []
    for m in read_tsv(ds / "manifests" / "manifest_tierA_A3.tsv"):
        fname = m["filename"]
        if not fname.endswith(".bam"):
            continue
        L = led.get(fname, {})
        bam = cdir / fname
        size = fsize(bam)
        ok = (L.get("status") == "verified" and (cdir / (fname + ".bai")).exists()
              and size == m["size_bytes"])
        if ok:
            status = "verified"
        elif (cdir / (fname + ".part")).exists():
            status = "downloading"
        else:
            status = L.get("status") or "pending"
        d = depth.get(fname)
        out.append(base_row(
            file_id=f"{PACBIO_COHORT}/{fname}", individual_id=normalise_individual_id(fname),
            source_sample_id=re.match(r"[A-Z]{2}\d{5}", fname).group(0),
            cohort=PACBIO_COHORT, accession=m["bucket_key"], platform="PACBIO",
            chemistry="", library="HiFi", basecaller="", data_type="aligned_bam",
            reference=REFERENCE, coverage_x=f"{d:.2f}" if d is not None else "",
            coverage_source="measured_chr20_depth" if d is not None else "unknown",
            rel_path=f"{PACBIO_COHORT}/{fname}", size_bytes=size or m["size_bytes"], md5="",
            status=status, date_added=iso_date(L.get("ts", "")),
            notes="tier=A3; one release per individual (W preferred over C where both exist)"))
    return out


def ont1019_rows(ds: Path) -> list[dict]:
    cdir = ds / ONT1019_COHORT
    integ_dir = ds / "logs" / "align" / "integrity"
    crams = {p.name.split(".")[1]: p for p in cdir.glob("*.cram")}
    out = []
    for m in read_tsv(ds / "manifests" / "manifest_PRJEB89727_odp.tsv"):
        acc, iid = m["run_accession"], m["individual_id"]
        cram = crams.get(acc)
        rec = None
        ip = integ_dir / f"{acc}.tsv"
        if ip.exists():
            f = ip.read_text().strip().split("\t")
            rec = {"sra_reads": f[2], "cram_reads": f[3], "verdict": f[4]}
        ok = (cram is not None and Path(str(cram) + ".crai").exists()
              and rec is not None and rec["verdict"] == "read_count_match"
              and rec["sra_reads"] == rec["cram_reads"])
        mt = dt.datetime.fromtimestamp(cram.stat().st_mtime).date().isoformat() if cram else ""
        notes = ["aligned in-house from the AWS SRA mirror .sra: minimap2 2.24-r1122 "
                 "-x map-ont -a --eqx -K 10g, samtools sort to CRAM"]
        if rec:
            notes.append(f"reads sra={rec['sra_reads']} cram={rec['cram_reads']}")
        out.append(base_row(
            file_id=f"{ONT1019_COHORT}/{cram.name if cram else acc}", individual_id=iid,
            source_sample_id=iid, cohort=ONT1019_COHORT, accession=acc, platform="ONT",
            chemistry="", library="", basecaller="", data_type="aligned_cram",
            reference=REFERENCE, coverage_x=m["est_coverage_x"], coverage_source=m["coverage_source"],
            rel_path=f"{ONT1019_COHORT}/{cram.name}" if cram else "",
            size_bytes=fsize(cram) if cram else "", md5="",
            status="verified" if ok else ("failed" if rec else "pending"),
            date_added=mt, notes="; ".join(notes)))
    # the two ENA-only runs: recorded, not acquired
    full = {r["run_accession"]: r for r in read_tsv(ds / "manifests" / "manifest_PRJEB89727.tsv")}
    for acc, iid in ENA_ONLY_RUNS.items():
        m = full.get(acc, {})
        out.append(base_row(
            file_id=f"{ONT1019_COHORT}/{acc}", individual_id=iid, source_sample_id=iid,
            cohort=ONT1019_COHORT, accession=acc, platform="ONT", data_type="fastq",
            reference="none", coverage_x=m.get("est_coverage_x", ""),
            coverage_source=m.get("coverage_source", "ena_base_count"),
            size_bytes=m.get("size_bytes", ""), md5=m.get("md5", ""), status="excluded",
            date_added="2026-10-27",
            notes="ENA-only FASTQ; not on the AWS SRA mirror; skipped by PI decision 2026-10-27"))
    return out


# ─── derived tables ──────────────────────────────────────────────────────────
def passes_floor(row: dict) -> bool:
    try:
        return float(row["coverage_x"]) >= COVERAGE_FLOOR_X
    except ValueError:
        return False


def individuals_table(rows: list[dict], meta) -> list[dict]:
    by = defaultdict(list)
    for r in rows:
        if r["status"] == "verified":
            by[r["individual_id"]].append(r)
    out = []
    for iid in sorted(by):
        rs = by[iid]
        covs = [float(r["coverage_x"]) for r in rs if r["coverage_x"]]
        md = meta.get(iid, {})
        out.append({
            "individual_id": iid,
            "n_verified_files": len(rs),
            "cohorts": ",".join(sorted({r["cohort"] for r in rs})),
            "platforms": ",".join(sorted({r["platform"] for r in rs})),
            "max_coverage_x": f"{max(covs):.2f}" if covs else "",
            "passes_17x_any": "yes" if any(passes_floor(r) for r in rs) else "no",
            "sex": md.get("Sex", ""), "subpopulation": md.get("SubPopulation", ""),
            "superpopulation": md.get("SuperPopulation", ""),
        })
    return out


def sample_sheets(rows: list[dict], ds: Path) -> dict[str, list[dict]]:
    """Per-cohort pipeline sheets. sample_id = individual_id, suffixed with the
    accession only if one person has two files in the same cohort."""
    sheets = defaultdict(list)
    for cohort in sorted({r["cohort"] for r in rows}):
        rs = [r for r in rows if r["cohort"] == cohort and r["status"] == "verified"]
        dup = Counter(r["individual_id"] for r in rs)
        for r in sorted(rs, key=lambda r: r["individual_id"]):
            aln = ds / r["rel_path"]
            idx = Path(str(aln) + (".crai" if aln.suffix == ".cram" else ".bai"))
            sid = r["individual_id"] if dup[r["individual_id"]] == 1 else \
                f"{r['individual_id']}_{Path(r['accession']).name.split('.')[0]}"
            sheets[cohort].append({
                "sample_id": sid, "alignment": str(aln), "index": str(idx),
                "individual_id": r["individual_id"], "platform": r["platform"],
                "coverage_x": r["coverage_x"], "coverage_source": r["coverage_source"],
                "passes_17x": "yes" if passes_floor(r) else "no",
            })
    return sheets


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--datasets", default="/mnt/sargasso1/CxSV/datasets", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    ds = a.datasets
    src = ds / "registry" / "sources"
    meta, cramino, depth = load_lrsc_metadata(src), load_cramino(src), load_pacbio_depth(src)

    rows = lrsc_rows(ds, meta, cramino) + pacbio_rows(ds, depth) + ont1019_rows(ds)
    ids = Counter(r["file_id"] for r in rows)
    dups = [k for k, n in ids.items() if n > 1]
    if dups:
        sys.exit(f"duplicate file_id(s): {dups[:5]}")
    inds = individuals_table(rows, meta)
    sheets = sample_sheets(rows, ds)

    print("cohort\tstatus\tfiles")
    for (c, s), n in sorted(Counter((r["cohort"], r["status"]) for r in rows).items()):
        print(f"{c}\t{s}\t{n}")
    v = [r for r in rows if r["status"] == "verified"]
    print(f"verified files: {len(v)}; distinct individuals: {len(inds)}; "
          f"individuals passing {COVERAGE_FLOOR_X:g}x in at least one file: "
          f"{sum(i['passes_17x_any'] == 'yes' for i in inds)}")
    multi = Counter(i["cohorts"].count(",") + 1 for i in inds)
    print("individuals by number of cohorts: " + ", ".join(f"{k}:{n}" for k, n in sorted(multi.items())))
    if a.dry_run:
        return 0
    reg = ds / "registry"
    write_atomic(reg / "sample_registry.tsv", COLUMNS, rows)
    write_atomic(reg / "individuals.tsv", list(inds[0].keys()) if inds else ["individual_id"], inds)
    for cohort, srows in sheets.items():
        write_atomic(reg / "sample_sheets" / f"{cohort}.tsv", list(srows[0].keys()), srows)
    print(f"wrote {reg}/sample_registry.tsv ({len(rows)} rows), individuals.tsv ({len(inds)}), "
          f"sample_sheets/ ({', '.join(f'{c}={len(s)}' for c, s in sheets.items())})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

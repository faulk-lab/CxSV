#!/usr/bin/env python3
"""Download every file in a manifest, verify it, and record the outcome.

One downloader serves every cohort in this project: the ENA FASTQ sets and the
S3 bucket BAM sets differ only in which manifest columns carry the URL, the
expected size and (where published) the checksum.

    python3 download_manifest.py manifest.tsv --dest /mnt/.../1kGP_ONT_1019_fastq \\
        --jobs 6 --ledger /mnt/.../logs/ena_1019.tsv

Design points that matter in practice
-------------------------------------
* **Resumable at two levels.** A file already in the ledger as ``verified`` is
  skipped entirely; a partially transferred file is resumed with an HTTP range
  request rather than restarted.
* **Verified, not assumed.** Size is always checked. Where the manifest carries
  an md5 it is recomputed from the bytes on disk; a size-only check cannot
  detect silent corruption, and ENA publishes md5s precisely so this is possible.
* **Parallel across files, not within one.** A single stream reached ~57 MB/s in
  testing while six concurrent streams reached ~285 MB/s, so concurrency is what
  fills the link - and per-file streams keep the ledger honest about which file
  is complete.
* **Interrupt-safe.** Each file is written to ``.part`` and only renamed into
  place after verification, so an interrupted run never leaves a truncated file
  looking like a finished one.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

LEDGER_COLUMNS = ["filename", "url", "expected_bytes", "actual_bytes",
                  "expected_md5", "actual_md5", "status", "ts", "note"]

_ledger_lock = threading.Lock()
FLOOR_BYTES = int(3e12)


def md5_of(path: Path, chunk: int = 16 << 20) -> str:
    """Streaming md5 so a 200 GB BAM does not need to fit in RAM."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_ledger(path: Path) -> dict[str, str]:
    """Map filename -> status for rows already recorded."""
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["filename"]] = row["status"]
    return out


def append_ledger(path: Path, row: dict) -> None:
    with _ledger_lock:
        new = not path.exists()
        with open(path, "a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=LEDGER_COLUMNS, delimiter="\t")
            if new:
                w.writeheader()
            w.writerow(row)


def fetch(url: str, part: Path, tries: int = 3) -> None:
    """Fetch ``url`` into ``part``, resuming a partial file if one exists."""
    for attempt in range(1, tries + 1):
        cmd = ["curl", "-sSL", "--fail", "--retry", "3", "--retry-delay", "5",
               "-C", "-", "-o", str(part), url]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode == 0:
            return
        # curl exits 33 when the server cannot honour a resume request; start over.
        if p.returncode == 33 and part.exists():
            part.unlink()
            continue
        if attempt == tries:
            raise RuntimeError(f"curl exit {p.returncode}: {p.stderr.strip()[:200]}")


_space_lock = threading.Lock()
_reserved = 0          # bytes promised to transfers currently in flight


def reserve_space(dest: Path, nbytes: int, floor_bytes: int) -> bool:
    """Reserve room for a transfer, or refuse if it would breach the floor.

    The share is used by other projects too, so free space is re-read from the
    filesystem for every file rather than budgeted once at start. In-flight
    transfers are counted as already spent: six concurrent 400 GB downloads
    would otherwise each see the same free space and all start.
    """
    global _reserved
    with _space_lock:
        st = os.statvfs(dest)
        free = st.f_bavail * st.f_frsize
        if free - _reserved - nbytes < floor_bytes:
            return False
        _reserved += nbytes
        return True


def release_space(nbytes: int) -> None:
    global _reserved
    with _space_lock:
        _reserved = max(0, _reserved - nbytes)


def one_file(row: dict, dest: Path, ledger: Path, verify_md5: bool,
             cols: dict[str, str]) -> tuple[str, str]:
    name = row[cols["filename"]]
    url = row[cols["url"]]
    exp_bytes = int(row[cols["size"]]) if row.get(cols["size"]) else 0
    exp_md5 = (row.get(cols["md5"]) or "").strip() if cols["md5"] else ""

    final = dest / name
    part = dest / (name + ".part")

    have = part.stat().st_size if part.exists() else 0
    need = max(0, exp_bytes - have)
    if not reserve_space(dest, need, FLOOR_BYTES):
        return name, "deferred_low_space"
    try:
        fetch(url, part)
        actual = part.stat().st_size
        if exp_bytes and actual != exp_bytes:
            status, note = "size_mismatch", f"expected {exp_bytes}, got {actual}"
            act_md5 = ""
        else:
            act_md5 = md5_of(part) if (verify_md5 and exp_md5) else ""
            if exp_md5 and act_md5 and act_md5 != exp_md5:
                status, note = "md5_mismatch", "checksum differs from ENA-published value"
            elif name.endswith(".bam") and not exp_md5 and subprocess.run(
                    ["samtools", "quickcheck", str(part)], capture_output=True).returncode != 0:
                # No published md5 for the S3 BAMs, so size alone cannot tell a
                # full-length but corrupt file from a good one. quickcheck reads
                # the header and the BGZF EOF block.
                status, note = "quickcheck_failed", "size matched but samtools quickcheck failed"
            else:
                part.rename(final)          # only now does it look complete
                status, note = "verified", ("" if exp_md5 else
                                            "size + samtools quickcheck (no published md5)"
                                            if name.endswith(".bam") else
                                            "size-only check (no published md5)")
        append_ledger(ledger, {
            "filename": name, "url": url, "expected_bytes": exp_bytes,
            "actual_bytes": actual, "expected_md5": exp_md5, "actual_md5": act_md5,
            "status": status, "ts": datetime.now(timezone.utc).isoformat(), "note": note,
        })
        release_space(need)
        return name, status
    except Exception as e:                   # network error, disk full, ...
        release_space(need)
        append_ledger(ledger, {
            "filename": name, "url": url, "expected_bytes": exp_bytes,
            "actual_bytes": part.stat().st_size if part.exists() else 0,
            "expected_md5": exp_md5, "actual_md5": "", "status": "failed",
            "ts": datetime.now(timezone.utc).isoformat(), "note": str(e)[:200],
        })
        return name, "failed"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest")
    p.add_argument("--dest", required=True)
    p.add_argument("--ledger", required=True)
    p.add_argument("--jobs", type=int, default=6,
                   help="concurrent file transfers (6 saturated the link in testing)")
    p.add_argument("--no-md5", action="store_true",
                   help="skip md5 recomputation (size check only) - not recommended")
    p.add_argument("--limit", type=int, default=0, help="download only the first N files")
    p.add_argument("--min-free-tb", type=float, default=3.0,
                   help="never let free space on the destination fall below this (TB, decimal)")
    p.add_argument("--only-passing", action="store_true",
                   help="restrict to rows with passes_coverage_floor == yes")
    p.add_argument("--col-filename", default="filename")
    p.add_argument("--col-url", default="url")
    p.add_argument("--col-size", default="size_bytes")
    p.add_argument("--col-md5", default="md5")
    a = p.parse_args()

    global FLOOR_BYTES
    FLOOR_BYTES = int(a.min_free_tb * 1e12)
    dest = Path(a.dest); dest.mkdir(parents=True, exist_ok=True)
    ledger = Path(a.ledger); ledger.parent.mkdir(parents=True, exist_ok=True)

    with open(a.manifest, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if a.only_passing:
        rows = [r for r in rows if r.get("passes_coverage_floor") == "yes"]

    done = load_ledger(ledger)
    outstanding = [r for r in rows
                   if done.get(r[a.col_filename]) != "verified"
                   and not (dest / r[a.col_filename]).exists()]
    # Report completed and limit-deferred separately - conflating them makes a
    # --limit test run look as though the whole cohort is already in hand.
    todo = outstanding[: a.limit] if a.limit else outstanding
    deferred = len(outstanding) - len(todo)

    planned = sum(int(r[a.col_size]) for r in todo if r.get(a.col_size))
    print(f"{len(rows)} in manifest, {len(rows) - len(outstanding)} already complete, "
          f"{len(todo)} to fetch now ({planned / 1e12:.2f} TB)"
          + (f", {deferred} deferred by --limit" if deferred else "")
          + f", {a.jobs} parallel", flush=True)

    cols = {"filename": a.col_filename, "url": a.col_url,
            "size": a.col_size, "md5": a.col_md5}
    tally: dict[str, int] = {}
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = [ex.submit(one_file, r, dest, ledger, not a.no_md5, cols) for r in todo]
        for i, f in enumerate(as_completed(futs), 1):
            name, status = f.result()
            tally[status] = tally.get(status, 0) + 1
            print(f"[{i}/{len(todo)}] {status}: {name}", flush=True)

    print("--- summary ---")
    for s, n in sorted(tally.items()):
        print(f"  {s}: {n}")
    if tally.get("deferred_low_space"):
        print(f"STOPPED SHORT: {tally['deferred_low_space']} file(s) deferred - free space "
              f"would have fallen below {a.min_free_tb} TB. Free space, then re-run; "
              "verified files are skipped automatically.")
    if set(tally) - {"verified"}:
        sys.exit(1)


if __name__ == "__main__":
    main()

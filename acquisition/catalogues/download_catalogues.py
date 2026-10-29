#!/usr/bin/env python3
"""Download external SV catalogues listed in sources.tsv, with checksums.

Each file lands in <root>/<catalogue>/<version>/<file>. For every catalogue
version a SOURCES.tsv is (re)written next to the files, one row per file:

    file  url  download_date_utc  size_bytes  published_checksum
    checksum_source  computed_md5  computed_sha256  verified

``verified`` is ``yes`` when a published checksum exists and matches,
``no-published-checksum`` when none is published (the computed sha256 is then
the record of what we hold), and the script exits non-zero on any mismatch.

Published-checksum keywords in sources.tsv:
  md5:<hex>     checksum quoted from the provider's page or manifest
  gcs           MD5 read from the Google Cloud Storage object metadata
  ucsc_md5sum   MD5 looked up in md5sum.txt of the same UCSC directory
  .             none published

Files already present with a matching published checksum are not
re-downloaded (``--force`` re-downloads). Standard library only, so it runs
in any Python >= 3.8. Downloads use curl (retries, resumable).

Usage:
  python3 download_catalogues.py --root /mnt/sargasso1/CxSV/datasets/external_catalogues \
      [--catalogue gnomAD-SV] [--dry-run] [--force]
"""
import argparse
import base64
import csv
import datetime as dt
import hashlib
import os
import subprocess
import sys
import urllib.request
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCES_COLS = ["file", "url", "download_date_utc", "size_bytes",
                "published_checksum", "checksum_source", "computed_md5",
                "computed_sha256", "verified"]


def read_sources(path):
    """Return the rows of sources.tsv as dicts (comment lines skipped)."""
    with open(path) as fh:
        lines = [l for l in fh if not l.startswith("#") and l.strip()]
    return list(csv.DictReader(lines, delimiter="\t"))


def head(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return {k.lower(): v for k, v in r.headers.items()}, r.headers.get_all("x-goog-hash") or []


def resolve_published(row):
    """Turn the sources.tsv checksum keyword into 'md5:<hex>' or '.'."""
    pc = row["published_checksum"]
    if pc == "gcs":
        _, goog = head(row["url"])
        for h in goog:
            for part in h.split(","):
                part = part.strip()
                if part.startswith("md5="):
                    return "md5:" + base64.b64decode(part[4:]).hex()
        return "."
    if pc == "ucsc_md5sum":
        md5_url = row["url"].rsplit("/", 1)[0] + "/md5sum.txt"
        with urllib.request.urlopen(md5_url, timeout=60) as r:
            for line in r.read().decode().splitlines():
                parts = line.split()
                if len(parts) == 2 and parts[1].lstrip("*") == row["file"]:
                    return "md5:" + parts[0]
        return "."
    return pc


def hashes(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            md5.update(chunk)
            sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


def fetch(url, dest):
    """curl with retries into dest.part, then atomic rename (NAS-safe)."""
    part = dest + ".part"
    cmd = ["curl", "-fsSL", "--retry", "5", "--retry-delay", "10",
           "-C", "-", "-o", part, url]
    subprocess.run(cmd, check=True)
    os.replace(part, dest)


def write_sources_tsv(path, rows):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SOURCES_COLS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


def read_existing_sources(path):
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return {r["file"]: r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True)
    ap.add_argument("--sources", default=os.path.join(HERE, "sources.tsv"))
    ap.add_argument("--catalogue", action="append",
                    help="limit to these catalogue names (repeatable)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print sizes (HTTP HEAD) and the total, download nothing")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    rows = read_sources(args.sources)
    if args.catalogue:
        rows = [r for r in rows if r["catalogue"] in args.catalogue]

    if args.dry_run:
        total = 0
        for r in rows:
            h, _ = head(r["url"])
            size = int(h.get("content-length", 0))
            total += size
            print(f"{r['catalogue']}\t{r['version']}\t{r['file']}\t{size}")
        print(f"TOTAL_BYTES\t{total}\t{total / 1e9:.2f} GB")
        return 0

    groups = OrderedDict()
    for r in rows:
        groups.setdefault((r["catalogue"], r["version"]), []).append(r)

    failures = 0
    for (cat, ver), grp in groups.items():
        outdir = os.path.join(args.root, cat, ver)
        os.makedirs(outdir, exist_ok=True)
        src_path = os.path.join(outdir, "SOURCES.tsv")
        existing = read_existing_sources(src_path)
        records = []
        for r in grp:
            dest = os.path.join(outdir, r["file"])
            published = resolve_published(r)
            prev = existing.get(r["file"])
            have = os.path.exists(dest)
            if have and not args.force:
                md5, sha = hashes(dest)
                ok = published == "." or published == "md5:" + md5
                if ok:
                    date = prev["download_date_utc"] if prev else "unknown"
                    print(f"[keep] {cat}/{ver}/{r['file']}")
                else:
                    have = False
            if not have or args.force:
                print(f"[get ] {cat}/{ver}/{r['file']}", flush=True)
                fetch(r["url"], dest)
                date = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
                md5, sha = hashes(dest)
            if published == ".":
                verified = "no-published-checksum"
            elif published == "md5:" + md5:
                verified = "yes"
            else:
                verified = "MISMATCH"
                failures += 1
                print(f"[FAIL] {r['file']}: published {published}, got md5:{md5}",
                      file=sys.stderr)
            records.append({
                "file": r["file"], "url": r["url"], "download_date_utc": date,
                "size_bytes": os.path.getsize(dest),
                "published_checksum": published,
                "checksum_source": r["checksum_source"],
                "computed_md5": md5, "computed_sha256": sha,
                "verified": verified})
        write_sources_tsv(src_path, records)
        print(f"[done] {src_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

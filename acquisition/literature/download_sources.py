#!/usr/bin/env python3
"""Download the literature-set source files listed in sources.tsv and verify their sha256.

    python3 acquisition/literature/download_sources.py --out /mnt/sargasso1/CxSV/resources/literature/_downloads/YYYYMMDD
        [--roles parsed,curated] [--chain-hg19]

Files whose role is 'parsed' are needed by build_literature_set.py; 'curated'
files are the published tables the curated/*.tsv rows were transcribed from (kept
for audit). Rows without a URL or checksum (context / not_obtained) are skipped.
A Europe PMC 'supplementaryFiles' URL returns a zip; the named member is extracted.
--chain-hg19 also fetches UCSC hg19ToHg38.over.chain.gz (checked against UCSC md5sum.txt).
Writes SOURCES.tsv (file, url, retrieved_utc, bytes, sha256, verified) in --out.
Exit status 1 if any checksum differs: a publisher replaced a file, so the
curation must be re-checked before the set is rebuilt.
"""
import argparse
import datetime as dt
import hashlib
import io
import os
import sys
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import literature_set as L  # noqa: E402

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/"


def fetch(url):
    with urllib.request.urlopen(url, timeout=300) as r:
        return r.read()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--roles", default="parsed,curated")
    ap.add_argument("--chain-hg19", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    roles = set(a.roles.split(","))
    rows, bad = [], 0
    for r in L.read_tsv(os.path.join(HERE, "sources.tsv")):
        if r["role"] not in roles or r["url"] == "." or r["sha256"] == ".":
            continue
        data = fetch(r["url"])
        if data[:2] == b"PK" and "supplementaryFiles" in r["url"]:
            data = zipfile.ZipFile(io.BytesIO(data)).read(r["file"])
        path = os.path.join(a.out, r["file"])
        with open(path, "wb") as fh:
            fh.write(data)
        got = hashlib.sha256(data).hexdigest()
        ok = got == r["sha256"]
        bad += not ok
        rows.append(dict(file=r["file"], url=r["url"], bytes=len(data), sha256=got, verified="yes" if ok else "NO",
                         retrieved_utc=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")))
        print(f"{r['file']}: {len(data)} bytes, sha256 {'ok' if ok else 'MISMATCH ' + got}", file=sys.stderr)
    if a.chain_hg19:
        name = "hg19ToHg38.over.chain.gz"
        data = fetch(UCSC + name)
        want = [ln.split()[0] for ln in fetch(UCSC + "md5sum.txt").decode().splitlines() if ln.endswith(name)]
        ok = bool(want) and hashlib.md5(data).hexdigest() == want[0]
        bad += not ok
        with open(os.path.join(a.out, name), "wb") as fh:
            fh.write(data)
        rows.append(dict(file=name, url=UCSC + name, bytes=len(data), sha256=hashlib.sha256(data).hexdigest(),
                         verified="yes (UCSC md5)" if ok else "NO",
                         retrieved_utc=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")))
    L.write_tsv(os.path.join(a.out, "SOURCES.tsv"), rows, ["file", "url", "retrieved_utc", "bytes", "sha256", "verified"])
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

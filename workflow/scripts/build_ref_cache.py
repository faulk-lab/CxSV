#!/usr/bin/env python3
"""Build an htslib MD5 reference cache from a FASTA (plain or gzip/bgzip).

For each contig, writes the upper-cased sequence (no newlines) to
``<cache>/<md5[0:2]>/<md5[2:4]>/<md5[4:]>``, the layout htslib resolves with
``REF_PATH=<cache>/%2s/%2s/%s``. The MD5 is computed exactly as the CRAM
``@SQ M5`` tag is (upper case, whitespace removed), so a CRAM made against this
FASTA decodes from the cache on any machine, whatever its header's UR path says.

Writes ``<cache>/md5_manifest.tsv`` (SN, LN, M5) last and atomically; its
presence marks a complete cache. Existing cache entries of the right length
are kept, so re-running is cheap.

Usage: build_ref_cache.py <reference.fa[.gz]> <cache_dir>
"""
import gzip
import hashlib
import os
import sys
from pathlib import Path


def contigs(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main(ref, cache):
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, seq in contigs(ref):
        seq = seq.upper().encode()
        md5 = hashlib.md5(seq).hexdigest()
        dest = cache / md5[:2] / md5[2:4] / md5[4:]
        if not (dest.exists() and dest.stat().st_size == len(seq)):
            dest.parent.mkdir(parents=True, exist_ok=True)
            tmp = dest.with_name(dest.name + ".tmp")
            tmp.write_bytes(seq)
            os.replace(tmp, dest)
        manifest.append((name, len(seq), md5))
        print(f"{name}\t{len(seq)}\t{md5}", file=sys.stderr)
    tmp = cache / "md5_manifest.tsv.tmp"
    with open(tmp, "w") as fh:
        fh.write("SN\tLN\tM5\n")
        for row in manifest:
            fh.write("\t".join(map(str, row)) + "\n")
    os.replace(tmp, cache / "md5_manifest.tsv")
    print(f"cached {len(manifest)} contigs from {ref} into {cache}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])

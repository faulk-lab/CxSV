#!/usr/bin/env python3
"""
fetch_sample_list.py
====================
Discovers all sample IDs in the 1000G-ONT S3 bucket under the CHM13
NAPU pipeline prefix and writes them to config/samples.txt.

Also probes one sample to auto-detect the exact BAM filename pattern
and prints the correct url_template / bai_template to paste into
config/remote_1000g.yaml.

This script requires NO AWS credentials — the bucket is public.
It uses the S3 XML listing API (HTTP GET with ?list-type=2 parameters).

Usage
-----
    python scripts/fetch_sample_list.py
    python scripts/fetch_sample_list.py --prefix ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/NAPU_PIPELINE/CHM13
    python scripts/fetch_sample_list.py --out config/samples.txt --max 100

Output
------
    config/samples.txt   — one sample ID per line
    Printed to stdout    — detected URL templates to copy into config/remote_1000g.yaml
"""

import argparse
import os
import sys
import urllib.request
import xml.etree.ElementTree as ET

BUCKET        = "1000g-ont"
BASE_HTTPS    = f"https://s3.amazonaws.com/{BUCKET}"
DEFAULT_PREFIX = "ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/NAPU_PIPELINE/CHM13"
NS            = "http://s3.amazonaws.com/doc/2006-03-01/"   # S3 XML namespace


def s3_list_objects(prefix: str, delimiter: str = "/",
                    max_keys: int = 1000) -> list[str]:
    """
    List S3 common prefixes (i.e. sub-directories) under prefix.
    Returns list of prefix strings, one per logical directory entry.
    Handles continuation tokens (pagination) automatically.
    """
    prefixes = []
    continuation = None

    while True:
        params = (
            f"?list-type=2"
            f"&prefix={urllib.parse.quote(prefix, safe='/')}"
            f"&delimiter={urllib.parse.quote(delimiter)}"
            f"&max-keys={max_keys}"
        )
        if continuation:
            params += f"&continuation-token={urllib.parse.quote(continuation)}"

        url = f"{BASE_HTTPS}/{params}"
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                body = resp.read()
        except Exception as e:
            sys.exit(f"ERROR fetching S3 listing: {e}\nURL: {url}")

        root = ET.fromstring(body)
        for cp in root.findall(f"{{{NS}}}CommonPrefixes"):
            p = cp.findtext(f"{{{NS}}}Prefix", "")
            if p:
                prefixes.append(p)

        is_truncated = root.findtext(f"{{{NS}}}IsTruncated", "false")
        if is_truncated.lower() == "true":
            continuation = root.findtext(f"{{{NS}}}NextContinuationToken")
        else:
            break

    return prefixes


def s3_list_keys(prefix: str, max_keys: int = 200) -> list[str]:
    """List object keys (files) directly under prefix."""
    import urllib.parse
    keys = []
    continuation = None

    while True:
        params = (
            f"?list-type=2"
            f"&prefix={urllib.parse.quote(prefix, safe='/')}"
            f"&max-keys={max_keys}"
        )
        if continuation:
            params += f"&continuation-token={urllib.parse.quote(continuation)}"

        url = f"{BASE_HTTPS}/{params}"
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                body = resp.read()
        except Exception as e:
            sys.stderr.write(f"WARN: could not list {prefix}: {e}\n")
            return keys

        root = ET.fromstring(body)
        for c in root.findall(f"{{{NS}}}Contents"):
            k = c.findtext(f"{{{NS}}}Key", "")
            if k:
                keys.append(k)

        is_truncated = root.findtext(f"{{{NS}}}IsTruncated", "false")
        if is_truncated.lower() == "true":
            continuation = root.findtext(f"{{{NS}}}NextContinuationToken")
        else:
            break

    return keys


def detect_bam_pattern(sample_id: str, prefix: str) -> tuple[str, str]:
    """
    Look inside one sample directory and return (bam_key, bai_key).
    Returns (None, None) if no BAM found.
    """
    sample_prefix = f"{prefix.rstrip('/')}/{sample_id}/"
    keys = s3_list_keys(sample_prefix)
    bam_key = next((k for k in keys if k.endswith(".bam")), None)
    bai_key = next((k for k in keys if k.endswith(".bam.bai")), None)
    return bam_key, bai_key


def main():
    import urllib.parse

    parser = argparse.ArgumentParser(
        description="Discover 1000G-ONT sample IDs from S3 and write samples.txt"
    )
    parser.add_argument("--prefix",  default=DEFAULT_PREFIX,
                        help="S3 key prefix (no leading/trailing slash)")
    parser.add_argument("--out",     default="config/samples.txt",
                        help="Output sample list file")
    parser.add_argument("--max",     type=int, default=1000,
                        help="Max samples to retrieve (default 1000 = all 100)")
    args = parser.parse_args()

    prefix = args.prefix.strip("/") + "/"
    print(f"Listing S3 directories under: s3://{BUCKET}/{prefix}")

    sub_prefixes = s3_list_objects(prefix, delimiter="/", max_keys=args.max)
    if not sub_prefixes:
        sys.exit(
            f"No sub-directories found under {prefix}.\n"
            f"Check that the prefix is correct in config/remote_1000g.yaml and this script."
        )

    # Sample ID = last component of each sub-prefix
    samples = []
    for p in sub_prefixes:
        # p looks like: ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/.../CHM13/HG00096/
        parts = p.rstrip("/").split("/")
        sid   = parts[-1]
        if sid:
            samples.append(sid)

    samples.sort()
    print(f"\nFound {len(samples)} samples.")

    # Auto-detect BAM filename pattern from the first sample
    print(f"\nProbing first sample to detect BAM filename pattern...")
    bam_key, bai_key = detect_bam_pattern(samples[0], prefix.rstrip("/"))

    if bam_key:
        # Derive the template by replacing the sample ID with {sample}
        bam_template = f"{BASE_HTTPS}/{bam_key}".replace(samples[0], "{sample}")
        bai_template = (f"{BASE_HTTPS}/{bai_key}".replace(samples[0], "{sample}")
                        if bai_key else bam_template + ".bai")
        print(f"\n✓ Detected BAM pattern. Add these to config/remote_1000g.yaml:\n")
        print(f"  url_template: \"{bam_template}\"")
        print(f"  bai_template: \"{bai_template}\"")
    else:
        print(
            f"\nWARN: Could not find a .bam file under {prefix}{samples[0]}/\n"
            f"       Check the bucket structure manually and update url_template\n"
            f"       in config/remote_1000g.yaml."
        )

    # Write samples.txt
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(f"# 1000G-ONT samples — {prefix}\n")
        fh.write(f"# Generated by scripts/fetch_sample_list.py\n")
        for s in samples:
            fh.write(s + "\n")

    print(f"\n✓ Wrote {len(samples)} sample IDs to: {args.out}")
    print(f"\nNext steps:")
    print(f"  1. Confirm url_template and bai_template in config/remote_1000g.yaml")
    print(f"  2. bash scripts/download_1000g_bams.sh")
    print(f"  3. bash scripts/run_pipeline.sh")


if __name__ == "__main__":
    main()

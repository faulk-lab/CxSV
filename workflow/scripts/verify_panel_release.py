#!/usr/bin/env python3
"""Verify a panel release against its release_manifest.tsv; refuse it if anything differs.

    verify_panel_release.py <release_dir> [--out check.txt]

Checks (workflow/lib/cxsv_release.py: verify): every listed file exists with
the listed size and sha256, no file is present that the manifest does not
list, and release.json is listed. On success prints the report and writes it to
--out; on failure prints every problem and exits 1 without writing --out.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_release as cr  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("release_dir")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    res = cr.verify(a.release_dir)
    lines = [f"panel release:     {os.path.abspath(a.release_dir)}"]
    meta = {}
    if os.path.isfile(os.path.join(a.release_dir, cr.METADATA)):
        try:
            with open(os.path.join(a.release_dir, cr.METADATA)) as fh:
                meta = json.load(fh)
        except ValueError as e:
            res["problems"].append(f"{cr.METADATA} is not valid JSON: {e}")
            res["ok"] = False
    lines += [f"release id:        {meta.get('release_id', 'unknown')}",
              f"test release:      {meta.get('test_release', 'unknown')}",
              f"source commit:     {meta.get('source', {}).get('commit', 'unknown')}",
              f"fingerprint:       {meta.get('calling_fingerprint', {}).get('sha256', 'unknown')}",
              f"manifest sha256:   {res['manifest_sha256']}",
              f"files verified:    {res['n_files']} ({res['bytes']} bytes)",
              # where this check ran (set by the container image / the bundle's run_cxsv.sh)
              f"runtime image:     {os.environ.get('CXSV_IMAGE') or 'none (host environment)'}",
              f"runtime image id:  {os.environ.get('CXSV_IMAGE_ID') or '-'}",
              f"runtime SIF:       {os.environ.get('CXSV_SIF_SHA256') or '-'}",
              f"bundle:            {os.environ.get('CXSV_BUNDLE_ID') or '-'}"]
    if res["ok"]:
        lines.append("status: VERIFIED")
    else:
        lines.append(f"status: FAILED ({len(res['problems'])} problem(s)) -- this release must not be used")
        lines += [f"  {p}" for p in res["problems"]]
    text = "\n".join(lines) + "\n"
    print(text, end="")
    if not res["ok"]:
        sys.exit(1)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(text)


if __name__ == "__main__":
    main()

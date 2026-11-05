#!/usr/bin/env python3
"""Reuse a built panel's upstream outputs in a new build root, as hard links.

panel-v1.1 rebuild (roadmap item 13). The new build re-runs only the rules
downstream of the calling change; everything upstream (caller VCFs, filtered
VCFs, merged_sv.vcf, reference checks, callability masks, .snf files, cohort
call sets, the Sniffles2 population VCF) is the old build's output, reused.
The lab NAS (CIFS) has no symlinks, so the files are HARD-linked: one inode,
no copy, the old file's mtime.

    link_upstream.py --old OLD/results/panel --new NEW/results/panel \
        --dryrun DRYRUN.txt --inventory NEW/provenance/reused_upstream.tsv [--apply]
    link_upstream.py --verify NEW/provenance/reused_upstream.tsv

What is linked: every file under OLD except
  - release/ (the old release is never touched; the new one is built fresh),
  - any file under a logs/ or benchmarks/ directory (a rule that ran in the
    new build could write to such a path in place, through the shared inode),
    except the logs named by --include-log: logs of held rules that a re-run
    rule reads (logs/svim.log, read by check_svim_output),
  - every output, log and benchmark of a job the new build's dry run plans
    (DRYRUN.txt: `snakemake -n` of the new code against the OLD tree, so the
    paths are OLD paths), and anything under a planned directory output.
So no job of the new build has a linked file as output or log: Snakemake
never writes to a linked inode, and the old tree cannot change.

The inventory (relative path, bytes, mtime_ns, inode, origin) is the
provenance record of the reused outputs. --verify re-stats the OLD files after
the build and fails if any size or mtime changed or a link is broken.
"""
import argparse
import csv
import os
import re
import sys

FIELDS = ("output:", "log:", "benchmark:")


def planned_paths(dryrun):
    """Outputs, logs and benchmarks of every job in a snakemake dry run."""
    out, rules = set(), {}
    with open(dryrun) as fh:
        for line in fh:
            m = re.match(r"^(?:local)?rule (\S+):", line)
            if m:
                rules[m.group(1)] = rules.get(m.group(1), 0) + 1
                continue
            s = line.strip()
            for f in FIELDS:
                if s.startswith(f):
                    out.update(os.path.normpath(p.strip()) for p in s[len(f):].split(",") if p.strip())
    return out, rules


def origin(rel):
    """Coarse origin of a reused file, for the provenance record."""
    parts = rel.split("/")
    if parts[0] == "samples":
        sub = "/".join(parts[3:])
        if sub.startswith("callers/"):
            return "callers (run_*/filter_*)"
        if sub.startswith("callability/"):
            return "callable_mask"
        if sub in ("merged_sv.vcf", "merged_vcf_list.txt"):
            return "merge_sv_callers"
        if sub.startswith("logs/"):
            return "held-rule log read downstream: " + sub
        if sub == "refcheck.ok":
            return "check_cram_reference"
        return "per-sample other: " + sub
    if parts[0] == "cohorts":
        return "cohort call sets/QC (rules not planned)"
    if parts[:2] == ["panel", "sniffles_population"]:
        return "panel_sniffles_population"
    if parts[:2] == ["panel", "callability"]:
        return "panel_callability_masks"
    return "other: " + rel


def under_planned_dir(p, planned):
    d = os.path.dirname(p)
    while len(d) > 1:
        if d in planned:
            return True
        d = os.path.dirname(d)
    return False


def build(a):
    old, new = os.path.abspath(a.old), os.path.abspath(a.new)
    planned, rules = planned_paths(a.dryrun)
    if not planned:
        sys.exit(f"no planned outputs found in {a.dryrun}")
    rows, skipped = [], {"release": 0, "logs/benchmarks": 0, "planned": 0}
    for root, dirs, files in os.walk(old):
        rel_root = os.path.relpath(root, old)
        if rel_root == "release" or rel_root.startswith("release/"):
            dirs[:] = []
            skipped["release"] += len(files)
            continue
        dirs.sort()
        for f in sorted(files):
            p = os.path.join(root, f)
            rel = os.path.relpath(p, old)
            comps = rel.split("/")[:-1]
            keep_log = comps[-1:] == ["logs"] and f in a.include_log
            if ("logs" in comps or "benchmarks" in comps) and not keep_log:
                skipped["logs/benchmarks"] += 1
            elif p in planned or under_planned_dir(p, planned):
                skipped["planned"] += 1
            else:
                st = os.stat(p)
                rows.append((rel, st.st_size, st.st_mtime_ns, st.st_ino, origin(rel)))
    print(f"planned jobs by rule: {dict(sorted(rules.items()))}")
    print(f"files to link: {len(rows)} ({sum(r[1] for r in rows) / 1e9:.1f} GB); skipped {skipped}")
    if not a.apply:
        return
    os.makedirs(os.path.dirname(os.path.abspath(a.inventory)), exist_ok=True)
    made = 0
    for rel, *_ in rows:
        src, dst = os.path.join(old, rel), os.path.join(new, rel)
        if os.path.exists(dst):
            if os.stat(dst).st_ino != os.stat(src).st_ino:
                sys.exit(f"{dst} exists and is not a link to the old file")
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        os.link(src, dst)
        made += 1
    tmp = a.inventory + ".tmp"
    with open(tmp, "w", newline="") as fh:
        fh.write(f"# old\t{old}\n# new\t{new}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(("path", "bytes", "mtime_ns", "inode", "origin"))
        w.writerows(rows)
    os.replace(tmp, a.inventory)
    print(f"linked {made} new ({len(rows)} total); inventory {a.inventory}")


def verify(inv):
    old = new = None
    with open(inv) as fh:
        for line in fh:
            if line.startswith("# old\t"):
                old = line.rstrip("\n").split("\t")[1]
            elif line.startswith("# new\t"):
                new = line.rstrip("\n").split("\t")[1]
    bad = n = 0
    with open(inv) as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            n += 1
            so, sn = os.stat(os.path.join(old, r["path"])), os.stat(os.path.join(new, r["path"]))
            if (so.st_size, so.st_mtime_ns) != (int(r["bytes"]), int(r["mtime_ns"])) or so.st_ino != sn.st_ino:
                bad += 1
                if bad <= 10:
                    print("CHANGED", r["path"])
    print(f"verify: {n} reused files, {bad} changed or unlinked")
    return bad == 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--old")
    ap.add_argument("--new")
    ap.add_argument("--dryrun")
    ap.add_argument("--inventory")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--include-log", action="append", default=[], metavar="NAME",
                    help="log file name of a held rule that a re-run rule reads (e.g. svim.log, read by "
                         "check_svim_output); linked although it is under logs/ (never written: its rule is held)")
    ap.add_argument("--verify", metavar="INVENTORY")
    a = ap.parse_args()
    if a.verify:
        sys.exit(0 if verify(a.verify) else 1)
    build(a)

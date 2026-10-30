#!/usr/bin/env python3
"""Per-rule cost table and per-platform cost model from measured test genomes.

    cost_model.py --cost-root ROOT --registry sample_registry.tsv --ref-speed S \
        --out-rules per_rule_cost.tsv --out-model cost_model.tsv

ROOT is a measure_cost.sh directory: results/panel/samples/<cohort>/<sample>/
(benchmarks, panel_build/cost.tsv, callability/mosdepth.mosdepth.summary.txt)
and logs/<sample>.{time.txt,stage.tsv,snakemake.log}.

per_rule_cost.tsv  one row per test sample and rule: wall s, threads, CPU s,
                   peak RSS, plus each sample's coverage (sheet coverage_x and
                   mosdepth mean depth at MAPQ >= 20), bytes read and kept.
cost_model.tsv     one row per platform group, consumed by make_shards.py:
                   cpu_h_per_x, wall_h_per_x (sample critical path), busy_threads
                   (CPU / wall), cpu_h_per_gb / wall_h_per_gb (alignment
                   decompression, fitted on the two PacBio genomes of equal
                   depth and 7x different file size, shared by BAM groups,
                   0 for CRAM), peak_rss_gb (SVIM + the largest of the other
                   readers + 0.5 GB, which is what can run at once with
                   --cores 6), kept_bytes_per_x, stage_mb_s, test_coverage_x,
                   ref_speed. Per-x values divide by the sheet's coverage_x,
                   the same number the planner multiplies by. Groups with two
                   test genomes average their per-x values and take the larger
                   peak RSS.
Assumption (stated in docs/PANEL_BUILD.md): cost is linear in coverage with
no intercept; SVIM's clustering may grow faster than linearly in deep data.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import shardlib  # noqa: E402
from make_shards import platform_group  # noqa: E402

READERS = ("run_sniffles", "run_cutesv", "run_svim", "callable_mask")


def rule_threads(snakemake_log: str) -> dict:
    out, cur = {}, None
    if not os.path.exists(snakemake_log):
        return out
    for ln in open(snakemake_log):
        m = re.match(r"^(?:local)?rule (\w+):", ln.strip())
        if m:
            cur = m.group(1)
        m = re.match(r"^threads: (\d+)", ln.strip())
        if m and cur:
            out[cur] = int(m.group(1))
    return out


def time_v(path: str) -> dict:
    d = {}
    if not os.path.exists(path):
        return d
    for ln in open(path):
        k, _, v = ln.strip().partition(": ")
        d[k] = v
    def hms(x):
        p = [float(y) for y in x.split(":")]
        return sum(v * 60 ** i for i, v in enumerate(reversed(p)))
    return {"wall_s": hms(d.get("Elapsed (wall clock) time (h:mm:ss or m:ss)", "0")),
            "cpu_s": float(d.get("User time (seconds)", 0)) + float(d.get("System time (seconds)", 0)),
            "max_rss_largest_proc_mb": float(d.get("Maximum resident set size (kbytes)", 0)) / 1024}


def mosdepth_mean(path: str) -> float:
    if not os.path.exists(path):
        return float("nan")
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["chrom"] == "total":
                return float(r["mean"])
    return float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cost-root", required=True)
    ap.add_argument("--registry", required=True)
    ap.add_argument("--ref-speed", type=float, required=True)
    ap.add_argument("--out-rules", required=True)
    ap.add_argument("--out-model", required=True)
    a = ap.parse_args()
    _, reg = shardlib.read_sheet(a.registry)
    chem = {r["file_id"]: r.get("chemistry", "") for r in reg}
    size = {r["file_id"]: int(r.get("size_bytes") or 0) for r in reg}
    rows, samples = [], []
    for f in sorted(glob.glob(f"{a.cost_root}/results/panel/samples/*/*/panel_build/cost.tsv")):
        sdir = os.path.dirname(os.path.dirname(f))
        recs = list(csv.DictReader(open(f), delimiter="\t"))
        if not recs:
            continue
        cohort, sample = recs[0]["cohort"], recs[0]["sample_id"]
        sheet = f"{a.cost_root}/sheets/{sample}/{cohort}.tsv"
        aln = shardlib.read_sheet(sheet)[1][0]["alignment"] if os.path.exists(sheet) else ""
        fid = f"{cohort}/{os.path.basename(aln)}"
        group = platform_group(cohort, chem.get(fid, ""))
        thr = rule_threads(f"{a.cost_root}/logs/{sample}.snakemake.log")
        tv = time_v(f"{a.cost_root}/logs/{sample}.time.txt")
        st = shardlib.read_sheet(f"{a.cost_root}/logs/{sample}.stage.tsv")[1] \
            if os.path.exists(f"{a.cost_root}/logs/{sample}.stage.tsv") else []
        depth = mosdepth_mean(f"{sdir}/callability/mosdepth.mosdepth.summary.txt")
        cov = float(recs[0]["coverage_x"])
        info = {"cohort": cohort, "sample_id": sample, "group": group, "coverage_x": cov,
                "mosdepth_mean_depth": round(depth, 2), "alignment_bytes": size.get(fid, 0),
                "kept_bytes": int(recs[0]["sample_bytes_kept"]),
                "stage_s": float(st[0]["stage_s"]) if st else float("nan"),
                "instance_wall_s": tv.get("wall_s", float("nan")), "instance_cpu_s": tv.get("cpu_s", float("nan"))}
        rss = {}
        for r in recs:
            cpu, wall, mrss = float(r["cpu_time"] or 0), float(r["s"] or 0), float(r["max_rss"] or 0)
            rss[r["rule"]] = mrss
            rows.append({**{k: info[k] for k in ("cohort", "sample_id", "group", "coverage_x", "mosdepth_mean_depth")},
                         "rule": r["rule"], "threads": thr.get(r["rule"], 1), "wall_s": round(wall, 1),
                         "cpu_s": round(cpu, 1), "max_rss_mb": round(mrss, 1),
                         "io_in_mb": r.get("io_in", ""), "io_out_mb": r.get("io_out", "")})
        info["bench_cpu_s"] = sum(float(r["cpu_time"] or 0) for r in recs)
        info["peak_rss_gb"] = (rss.get("run_svim", 0) + max(rss.get(k, 0) for k in READERS if k != "run_svim")) / 1024 + 0.5
        samples.append(info)
    if not samples:
        sys.exit(f"no cost records under {a.cost_root}")
    for s in samples:   # sample totals as rule "ALL"
        rows.append({"cohort": s["cohort"], "sample_id": s["sample_id"], "group": s["group"],
                     "coverage_x": s["coverage_x"], "mosdepth_mean_depth": s["mosdepth_mean_depth"],
                     "rule": "ALL (instance)", "threads": "", "wall_s": round(s["instance_wall_s"], 1),
                     "cpu_s": round(s["instance_cpu_s"], 1), "max_rss_mb": round(s["peak_rss_gb"] * 1024, 1),
                     "io_in_mb": round(s["alignment_bytes"] / 1e6), "io_out_mb": round(s["kept_bytes"] / 1e6)})
    shardlib.write_sheet(a.out_rules, list(rows[0].keys()), rows)
    by = defaultdict(list)
    for s in samples:
        by[s["group"]].append(s)

    def cpu_s(s):
        return max(s["instance_cpu_s"], s["bench_cpu_s"])

    # Cost = coverage term (per-read work) + bytes term (decompressing the
    # alignment, done once per reader). A group measured on two genomes of the
    # same depth but different file size (PacBio: with / without kinetics
    # tags) gives the bytes term; other BAM groups share it; CRAM groups
    # (one genome, file size ~ coverage) get none.
    def fit(ss, y):
        pts = sorted({(s["alignment_bytes"] / 1e9, y(s), s["coverage_x"]) for s in ss})
        if len(pts) >= 2 and abs(pts[-1][0] - pts[0][0]) > 1:
            (g0, y0, _), (g1, y1, _) = pts[0], pts[-1]
            return (y1 - y0) / (g1 - g0)
        return None
    bam_per_gb = {}
    for key, y in (("cpu", cpu_s), ("wall", lambda s: s["instance_wall_s"])):
        fits = [f for g, ss in by.items() for f in [fit(ss, y)] if f is not None and f > 0]
        bam_per_gb[key] = sum(fits) / len(fits) if fits else 0.0
    model = []
    for g, ss in sorted(by.items()):
        n = len(ss)
        is_cram = all(s["cohort"] == "1kGP_ONT_1019_CHM13" for s in ss)
        per_gb = {k: (0.0 if is_cram else v) for k, v in bam_per_gb.items()}
        cpu_x = [(cpu_s(s) - per_gb["cpu"] * s["alignment_bytes"] / 1e9) / s["coverage_x"] for s in ss]
        wall_x = [(s["instance_wall_s"] - per_gb["wall"] * s["alignment_bytes"] / 1e9) / s["coverage_x"] for s in ss]
        model.append({"group": g, "test_samples": ",".join(f"{s['cohort']}/{s['sample_id']}" for s in ss),
                      "test_coverage_x": round(sum(s["coverage_x"] for s in ss) / n, 2),
                      "cpu_h_per_x": round(sum(cpu_x) / n / 3600, 6),
                      "cpu_h_per_gb": round(per_gb["cpu"] / 3600, 7),
                      "wall_h_per_x": round(sum(wall_x) / n / 3600, 6),
                      "wall_h_per_gb": round(per_gb["wall"] / 3600, 7),
                      "busy_threads": round(sum(cpu_s(s) / s["instance_wall_s"] for s in ss) / n, 2),
                      "peak_rss_gb": round(max(s["peak_rss_gb"] for s in ss), 2),
                      "kept_bytes_per_x": round(sum(s["kept_bytes"] / s["coverage_x"] for s in ss) / n),
                      "stage_mb_s": round(sum(s["alignment_bytes"] / 1e6 / s["stage_s"] for s in ss
                                              if s["stage_s"] == s["stage_s"] and s["stage_s"] > 1) /
                                          max(1, sum(1 for s in ss if s["stage_s"] == s["stage_s"] and s["stage_s"] > 1)), 1),
                      "ref_speed": a.ref_speed})
    shardlib.write_sheet(a.out_model, list(model[0].keys()), model)
    for m in model:
        print(m)


if __name__ == "__main__":
    main()

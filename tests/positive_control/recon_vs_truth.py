#!/usr/bin/env python3
"""Validate the structural reconstruction (workflow/lib/cxsv_recon.py) against the
positive-control truth: does the derivative drawn from a run's outputs reproduce the planted
segment order and orientation?

For every Tier 1 locus of a positive-control run that overlaps a planted event with a
derivative path in the truth (events.yaml, add-on addon.yaml files), the locus is
reconstructed exactly as in the clinician report and compared with the truth:

  exact            the drawn derivative(s) equal the truth (segment order and orientation)
  in_alternatives  not the drawn one, but among the equally consistent orders listed
  mismatch         neither

Both sides are first reduced to the same resolution: reference segments that the truth
never separates (always adjacent in reference order, same copy number; e.g. a planted
segment A that stays joined to the left flank) are merged, and reconstructed segments are
named by the merged truth segment they overlap most. Truth events without a derivative path
(L03, windowed chromothripsis) are skipped and listed.

  recon_vs_truth.py --run <run>/results/patient/samples/<cohort>/<sample> \
      --truth events.yaml [--truth addon.yaml ...] --genes genes.final.bed \
      --out-tsv recon_vs_truth.tsv [--out-figure recon_vs_truth.png]

The comparison functions are also used by tests/unit/test_recon.py on truth-derived inputs.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "workflow", "lib"))
import cxsv_recon as R        # noqa: E402

FLANK_NAMES = ("p", "q")


def is_flank(n):
    return n in FLANK_NAMES or n.startswith(("p_", "q_"))


def parse(el):
    return el[:-1], el[-1]


def flip(s):
    return "-" if s == "+" else "+"


def truth_events(paths):
    """Events with segments and paths from events.yaml / addon.yaml files."""
    import yaml
    out = []
    for p in paths:
        d = yaml.safe_load(open(p))
        evs = d["events"] if "events" in d else [d]
        out += [e for e in evs if e.get("segments") and e.get("paths")]
    seen, keep = set(), []
    for e in out:
        if e["id"] not in seen:
            seen.add(e["id"])
            keep.append(e)
    return keep


def _ref_order(segs):
    return sorted(segs, key=lambda n: (R._natural(segs[n][0]), segs[n][1]))


def merge_unbroken(segs, paths):
    """(merged segments {name: [chrom, start, end]}, paths in merged names, name map).

    Neighbours x|y (x.end == y.start) are merged when every occurrence of x in every path is
    followed by y+ (x+) or preceded by y- (x-), every y is the mirror, and both occur equally
    often: the truth never separates them."""
    segs = {k: list(v) for k, v in segs.items()}
    paths = {k: [parse(e) for e in v] for k, v in paths.items()}
    name_of = {k: k for k in segs}
    changed = True
    while changed:
        changed = False
        order = _ref_order(segs)
        for x, y in zip(order, order[1:]):
            X, Y = segs[x], segs[y]
            if X[0] != Y[0] or X[2] != Y[1]:
                continue
            if (is_flank(x) and is_flank(y)):
                continue
            ok, nx, ny = True, 0, 0
            for p in paths.values():
                for k, (n, s) in enumerate(p):
                    if n == x:
                        nx += 1
                        nxt = p[k + 1] if k + 1 < len(p) else None
                        prv = p[k - 1] if k > 0 else None
                        if not ((s == "+" and nxt == (y, "+")) or (s == "-" and prv == (y, "-"))):
                            ok = False
                    if n == y:
                        ny += 1
                        nxt = p[k + 1] if k + 1 < len(p) else None
                        prv = p[k - 1] if k > 0 else None
                        if not ((s == "+" and prv == (x, "+")) or (s == "-" and nxt == (x, "-"))):
                            ok = False
            if not ok or nx != ny:
                continue
            keep = x if is_flank(x) else (y if is_flank(y) else x)
            drop = y if keep == x else x
            segs[keep] = [X[0], X[1], Y[2]]
            del segs[drop]
            for k in name_of:
                if name_of[k] == drop:
                    name_of[k] = keep
            for d, p in paths.items():
                q = []
                for n, s in p:
                    if n == drop:
                        continue
                    q.append((n, s))
                paths[d] = q
            changed = True
            break
    return segs, paths, name_of


def canonical_paths(paths, drop_normal=False):
    """Set of paths, each read from its p flank in + if possible (else the smaller reading)."""
    out = []
    for p in paths:
        fw = tuple(p)
        rv = tuple((n, flip(s)) for n, s in reversed(p))
        def key(t):
            return (not (t and t[0][0].startswith("p") and t[0][1] == "+"), t)
        out.append(min(fw, rv, key=key))
    return tuple(sorted(out))


def to_truth_names(res_segments, path, merged):
    """A reconstructed path (['A+', ...]) in merged-truth names; consecutive pieces of one truth
    segment are collapsed."""
    by_name = {s["name"]: s for s in res_segments}
    def best(n):
        s = by_name[n]
        if s["kind"] != "internal":
            fl = s["kind"]
            cand = [k for k in merged if is_flank(k) and merged[k][0] == s["chrom"] and k.startswith(fl)]
            return cand[0] if cand else n
        bestk, bestov = None, 0
        for k, (c, a, b) in merged.items():
            if c != s["chrom"]:
                continue
            ov = min(b, s["end"]) - max(a, s["start"])
            if ov > bestov:
                bestk, bestov = k, ov
        return bestk or f"?{n}"
    out, prev = [], None
    for e in path:
        n, st = parse(e)
        t = best(n)
        if out and out[-1] == (t, st) and prev != n:
            prev = n
            continue          # two reconstructed pieces of one truth segment (an extra boundary)
        out.append((t, st))
        prev = n
    return out


def compare(res, event):
    """Outcome of one reconstruction against one truth event (dict)."""
    msegs, mpaths, _ = merge_unbroken(event["segments"], event["paths"])
    truth = canonical_paths(mpaths.values())
    chroms = {s["chrom"] for s in res.get("segments", [])}
    # a single-locus reconstruction of a multi-chromosome event covers only its chromosomes
    truth_here = tuple(p for p in truth if {msegs[n][0] for n, _ in p if n in msegs} & chroms) or truth
    drawn = canonical_paths([to_truth_names(res["segments"], p, msegs) for p in res.get("paths", {}).values()])
    alts = [canonical_paths([to_truth_names(res["segments"], p, msegs) for p in alt])
            for alt in res.get("alternatives", [])]

    def same(a, b):
        return set(a) == set(b)
    if same(drawn, truth_here):
        outcome = "exact"
    elif any(same(a, truth_here) for a in alts):
        outcome = "in_alternatives"
    else:
        outcome = "mismatch"
    fmt = lambda P: " | ".join(" ".join(f"{n}{s}" for n, s in p) for p in P)
    return {"outcome": outcome, "truth": fmt(truth_here), "drawn": fmt(drawn),
            "alternatives": " || ".join(fmt(a) for a in alts), "status": res.get("status"),
            "n_equivalent": res.get("n_equivalent", 1), "reasons": "; ".join(res.get("reasons", []))}


def event_for_locus(row, events):
    c, s, e = row["Chr"], int(row["Start"]), int(row["End"])
    best, bov = None, 0
    for ev in events:
        for n, (ch, a, b) in ev["segments"].items():
            if is_flank(n) and not n.startswith(("p_", "q_")):
                continue
            if ch == c:
                ov = min(b, e) - max(a, s)
                if ov > bov:
                    best, bov = ev, ov
    return best


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run", required=True, help="results/patient/samples/<cohort>/<sample> of a PC run")
    ap.add_argument("--truth", action="append", required=True)
    ap.add_argument("--genes", default=None)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json", default=None)
    ap.add_argument("--out-figure", default=None)
    a = ap.parse_args(argv)
    d = a.run
    with open(os.path.join(d, "scoring", "cxsv_scores.tsv"), newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    t1 = sorted([r for r in rows if r.get("Tier") == "Tier1_high_priority"], key=lambda r: -float(r["Total_Score"]))
    params = json.load(open(os.path.join(d, "scoring", "scoring_params.json")))
    stats = json.load(open(os.path.join(d, "cn", "stats.json")))
    sex = R.sc.resolve_sex(stats)["sex"]
    callers = {c: os.path.join(d, "callers", f) for c, f in R.CALLER_FILES.items()}
    recs = R.reconstruct_run(t1, os.path.join(d, "junctions", "junctions.tsv"), callers,
                             os.path.join(d, "cxsv", "cxsv_only.vcf"), os.path.join(d, "scoring", "locus_depth.tsv"),
                             os.path.join(d, "cn", "segments.tsv"), a.genes, params, sex)
    events = truth_events(a.truth)
    out, pairs = [], []
    for r, res in zip(t1, recs):
        ev = event_for_locus(r, events)
        row = {"locus": r["Locus_ID"], "score": r["Total_Score"], "class": r["Cmplx_Class"],
               "event": ev["id"] if ev else "none", "event_name": ev.get("name", "") if ev else ""}
        if ev is None:
            row.update(outcome="no_truth_path", truth="", drawn=" | ".join(" ".join(p) for p in res["paths"].values()),
                       alternatives="", status=res["status"], n_equivalent=res.get("n_equivalent", 1),
                       reasons="; ".join(res.get("reasons", [])))
        else:
            row.update(compare(res, ev))
            pairs.append((r, res, ev, row))
        out.append(row)
    cols = ["locus", "score", "class", "event", "event_name", "outcome", "status", "n_equivalent", "truth", "drawn",
            "alternatives", "reasons"]
    with open(a.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    n = {k: sum(1 for x in out if x["outcome"] == k) for k in ("exact", "in_alternatives", "mismatch", "no_truth_path")}
    print(json.dumps({"tier1_loci": len(out), **n}))
    if a.out_json:
        json.dump({"rows": out, "counts": n}, open(a.out_json, "w"), indent=1)
    if a.out_figure and pairs:
        import cxsv_recon_figure as RF
        RF.truth_vs_recon_figure(pairs, a.out_figure)
    return out


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Score the patient workflow's output on the positive control against the truth set.

Junction level (truth_junctions.tsv, one row per novel adjacency):
  exact    a merged call has one breakpoint within --slop of each junction end
           (either order; BND mates parsed from ALT, others use INFO/END);
  partial  a merged call has a breakpoint within --slop of either junction end
           (e.g. an INS or a one-sided call);
  callers  SUPP_VEC of the matching call(s), in SURVIVOR file order:
           cuteSV, SVIM, Sniffles2 (so 001 = Sniffles2 only).
  panel    the patient junction (private/junction_panel_status.tsv) matching the
           truth junction on both ends within --slop, and its panel carriers
           (breakpoint-level matching against the panel's junction sites).
Event level (truth_events.bed + events.yaml expectations):
  loci     patient CxSV loci (private/cxsv_panel_status.tsv) overlapping any
           span of the event, with class, criteria, carriers, panel status and
           region-overlap (hotspot) carriers;
  verdict  PASS if every expectation in events.yaml holds:
             expect.cxsv true, panel_status private (default) -> >=1 overlapping
               CxSV locus, none of them public (private, or shared_junctions:
               junctions seen in the panel but no individual carries the
               rearrangement), 0 rearrangement carriers;
             expect.cxsv true, panel_status public -> >=1 overlapping CxSV locus,
               all public with exactly expect.panel_carriers carriers (spike-in E13);
             expect.cxsv false -> no overlapping CxSV locus (simple SV controls)
               AND every truth junction matched on breakpoints by a patient
               junction with 0 panel carriers (not region overlap);
           and >=1 junction detected exactly (or partially for templated INS).
Add-on events (E13) keep their truth outside the release: --extra-truth DIR adds
DIR/truth_junctions.tsv and DIR/truth_events.bed.
Outside the events: CxSV loci not overlapping any event span are listed as the
host's own CxSV-like calls (expected; a 10-genome panel cannot make most public).

Scoring (roadmap item 5; SAMPLE_DIR/scoring/cxsv_scores.tsv, when present): each
event's best-scoring overlapping locus gives score, tier and top-3 contributions;
expect.score 'top' -> that locus is in the top tier (the first tier of
config/scoring.yaml), 'below_top' -> not in the top tier and below every 'top'
event, 'not_scored' -> no scored locus. The score expectation is part of the verdict.
class_match compares with expect.ideal_class; class_match_before with
expect.ideal_class_before (the label before roadmap item 5) when given.

Known failures (--known-failures TSV: event<TAB>reason) turn FAIL into XFAIL; a
listed event that passes is XPASS (update the list). The regression counts FAIL
and XPASS as failures, so it catches new problems without masking old ones.
class_match reports whether the pipeline's class equals expect.ideal_class.

Usage: score_pos_control.py TRUTH_DIR SAMPLE_DIR events.yaml OUT_PREFIX [--slop 1000]
                            [--known-failures known_failures.tsv] [--extra-truth DIR ...]
Writes OUT_PREFIX.junctions.tsv, OUT_PREFIX.events.tsv, OUT_PREFIX.offtarget.tsv.
"""
import argparse
import csv
import re

import yaml

BND_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def parse_calls(vcf):
    calls = []
    for line in open(vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        info = dict(kv.split("=", 1) if "=" in kv else (kv, True) for kv in f[7].split(";"))
        c1, p1 = f[0], int(f[1])
        m = BND_RE.search(f[4])
        if m:
            c2, p2 = m.group(1), int(m.group(2))
        else:
            c2, p2 = info.get("CHR2", c1), int(info.get("END", p1))
        calls.append(dict(id=f[2], svtype=info.get("SVTYPE", "?"), c1=c1, p1=p1, c2=c2, p2=p2,
                          supp=info.get("SUPP_VEC", "?"), svlen=info.get("SVLEN", ".")))
    return calls


def verdict(ok, is_known):
    if is_known:
        return "XPASS" if ok else "XFAIL"
    return "PASS" if ok else "FAIL"


def near(c, p, c0, p0, slop):
    return c == c0 and abs(p - p0) <= slop


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("truth"); ap.add_argument("sample"); ap.add_argument("events"); ap.add_argument("out")
    ap.add_argument("--slop", type=int, default=1000)
    ap.add_argument("--known-failures", default=None)
    ap.add_argument("--extra-truth", action="append", default=[])
    a = ap.parse_args()

    calls = parse_calls(f"{a.sample}/merged_sv.vcf")
    evs = {e["id"]: e for e in yaml.safe_load(open(a.events))["events"]}
    truth_dirs = [a.truth] + a.extra_truth
    spans = {}
    for d in truth_dirs:
        for line in open(f"{d}/truth_events.bed"):
            c, s, e, eid = line.split("\t")[:4]
            spans.setdefault(eid, []).append((c, int(s), int(e)))
    loci = list(csv.DictReader(open(f"{a.sample}/private/cxsv_panel_status.tsv"), delimiter="\t"))
    pj = list(csv.DictReader(open(f"{a.sample}/private/junction_panel_status.tsv"), delimiter="\t"))
    truth_j = [j for d in truth_dirs for j in csv.DictReader(open(f"{d}/truth_junctions.tsv"), delimiter="\t")]
    missing = sorted(set(evs) - set(spans))
    if missing:
        raise SystemExit(f"events without truth: {missing} (add-on truth dirs: --extra-truth)")

    def bp_match(c1, p1, c2, p2):
        """Patient junctions matching a truth junction on both ends (either order)."""
        return [k for k in pj if (near(k["chrom1"], int(k["pos1"]), c1, p1, a.slop)
                                  and near(k["chrom2"], int(k["pos2"]), c2, p2, a.slop))
                or (near(k["chrom1"], int(k["pos1"]), c2, p2, a.slop)
                    and near(k["chrom2"], int(k["pos2"]), c1, p1, a.slop))]

    jrows, jsum, jpanel = [], {}, {}
    for j in truth_j:
        c1, p1, c2, p2 = j["chrom1"], int(j["pos1"]), j["chrom2"], int(j["pos2"])
        exact = [k for k in calls if (near(k["c1"], k["p1"], c1, p1, a.slop) and near(k["c2"], k["p2"], c2, p2, a.slop))
                 or (near(k["c1"], k["p1"], c2, p2, a.slop) and near(k["c2"], k["p2"], c1, p1, a.slop))]
        part = [k for k in calls if k not in exact and any(near(c, p, cc, pp, a.slop)
                for c, p in ((k["c1"], k["p1"]), (k["c2"], k["p2"])) for cc, pp in ((c1, p1), (c2, p2)))]
        status = "exact" if exact else ("partial" if part else "missed")
        best = exact or part
        bm = bp_match(c1, p1, c2, p2)
        pc = max((int(k["panel_carriers"]) for k in bm), default=None)
        jpanel.setdefault(j["event"], []).append(pc)
        jrows.append(dict(event=j["event"], junction=f"{j['left']}|{j['right']}",
                          truth=f"{c1}:{p1}-{c2}:{p2}", status=status,
                          calls=";".join(f"{k['svtype']}:{k['c1']}:{k['p1']}-{k['c2']}:{k['p2']}" for k in best[:3]) or ".",
                          supp_vec=";".join(sorted({k["supp"] for k in best})) or ".",
                          patient_junction=";".join(f"{k['svtype']}:{k['chrom1']}:{k['pos1']}-{k['chrom2']}:{k['pos2']}"
                                                    for k in bm[:3]) or ".",
                          panel_carriers="." if pc is None else pc))
        s = jsum.setdefault(j["event"], {"exact": 0, "partial": 0, "missed": 0})
        s[status] += 1

    def overlaps(r, eid, pad=0):
        return any(r["Chr"] == c and int(r["Start"]) <= e + pad and int(r["End"]) >= s - pad for c, s, e in spans[eid])

    known = {}
    if a.known_failures:
        for line in open(a.known_failures):
            if line.strip() and not line.startswith("#"):
                eid, reason = (line.rstrip("\n").split("\t", 1) + [""])[:2]
                known[eid] = reason
    import os
    score_path = f"{a.sample}/scoring/cxsv_scores.tsv"
    scores = list(csv.DictReader(open(score_path), delimiter="\t")) if os.path.exists(score_path) else None
    top_tier = "Tier1_high_priority"        # first tier of config/scoring.yaml

    def best_score(eid):
        if scores is None:
            return None
        ov = [r for r in scores if overlaps(r, eid)]
        return max(ov, key=lambda r: float(r["Total_Score"])) if ov else None

    top_scores = [float(b["Total_Score"]) for e, ev in evs.items()
                  if ev.get("expect", {}).get("score") == "top" and (b := best_score(e))]
    erows = []
    for eid, ev in evs.items():
        ov = [r for r in loci if overlaps(r, eid)]
        exp = ev.get("expect", {})
        js = jsum.get(eid, {"exact": 0, "partial": 0, "missed": 0})
        n = sum(js.values())
        detected = js["exact"] > 0 or js["partial"] > 0
        want = exp.get("panel_status", "private")
        jp = jpanel.get(eid, [])
        bp_ok = bool(jp) and all(x is not None for x in jp)
        if exp.get("cxsv") and want == "public":
            ok_cx = bool(ov) and all(r["Panel_Status"] == "public"
                                     and int(r["Panel_Carriers"]) == int(exp.get("panel_carriers", 1)) for r in ov)
        elif exp.get("cxsv"):
            ok_cx = bool(ov) and all(r["Panel_Status"] != "public" and int(r["Panel_Carriers"]) == 0 for r in ov)
        else:
            # simple SV controls: not a CxSV, and matched on breakpoints (not region overlap):
            # every truth junction has a patient junction with 0 panel carriers
            ok_cx = not ov and bp_ok and all(x == 0 for x in jp)
        b = best_score(eid)
        want_score = exp.get("score")
        if scores is None or want_score is None:
            ok_score = True
        elif want_score == "top":
            ok_score = b is not None and b["Tier"] == top_tier
        elif want_score == "below_top":
            ok_score = b is not None and b["Tier"] != top_tier and bool(top_scores) \
                and float(b["Total_Score"]) < min(top_scores)
        else:
            ok_score = b is None
        before = exp.get("ideal_class_before", exp.get("ideal_class"))
        erows.append(dict(event=eid, name=ev["name"], gene=ev.get("gene") or "intergenic",
                          expect_cxsv=exp.get("cxsv"), ideal_class=exp.get("ideal_class", "."),
                          junctions=n, exact=js["exact"], partial=js["partial"], missed=js["missed"],
                          cxsv_loci=len(ov),
                          classes=";".join(sorted({r["Cmplx_Class"] for r in ov})) or ".",
                          criteria=";".join(sorted({r["CxSV_Criteria"] for r in ov})) or ".",
                          expect_status=want if exp.get("cxsv") else "not_cxsv",
                          carriers=";".join(sorted({str(r["Panel_Carriers"]) for r in ov})) or ".",
                          panel_status=";".join(sorted({r["Panel_Status"] for r in ov})) or ".",
                          region_carriers=";".join(sorted({str(r["Panel_Region_Carriers"]) for r in ov})) or ".",
                          hotspot=";".join(sorted({r["Panel_Hotspot"] for r in ov})) or ".",
                          junctions_matched_bp=f"{sum(x is not None for x in jp)}/{len(jp)}",
                          junction_panel_carriers=",".join("." if x is None else str(x) for x in jp) or ".",
                          class_match=("yes" if exp.get("ideal_class") in {r["Cmplx_Class"] for r in ov}
                                       else ("n/a" if not exp.get("cxsv") else "no")),
                          ideal_class_before=before or ".",
                          class_match_before=("yes" if before in {r["Cmplx_Class"] for r in ov}
                                              else ("n/a" if not exp.get("cxsv") else "no")),
                          expect_score=want_score or ".",
                          score=b["Total_Score"] if b else ".", tier=b["Tier"] if b else ".",
                          rank=b["Rank"] if b else ".", top_features=b["Top_Features"] if b else ".",
                          verdict=verdict(detected and ok_cx and ok_score, eid in known)))
    off = [r for r in loci if not any(overlaps(r, eid, pad=a.slop) for eid in spans)]

    def write(path, rows):
        with open(path, "w", newline="") as fh:
            if not rows:
                fh.write("none\n")
                return
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(rows)

    write(a.out + ".junctions.tsv", jrows)
    write(a.out + ".events.tsv", erows)
    write(a.out + ".offtarget.tsv", [{k: r[k] for k in ("Chr", "Start", "End", "SV_Types", "Cmplx_Class",
                                                         "CxSV_Criteria", "Panel_Junctions", "Panel_Junctions_Seen",
                                                         "Panel_Carriers", "Panel_Region_Carriers", "Panel_Hotspot",
                                                         "Panel_Status")}
                                     for r in off])
    for r in erows:
        print(f"{r['event']}\t{r['verdict']}{' (' + known[r['event']] + ')' if r['event'] in known else ''}\tjunctions {r['exact']}+{r['partial']}p/{r['junctions']}\t"
              f"cxsv_loci={r['cxsv_loci']} class={r['classes']} (expected {r['ideal_class']}, match {r['class_match']}) "
              f"status={r['panel_status']} (expected {r['expect_status']}) carriers={r['carriers']} "
              f"region_carriers={r['region_carriers']} bp_matched={r['junctions_matched_bp']} "
              f"junction_panel_carriers={r['junction_panel_carriers']} score={r['score']} tier={r['tier']} "
              f"(expected {r['expect_score']})")
    cx = [r for r in erows if r["class_match"] != "n/a"]
    print(f"class concordance: {sum(r['class_match'] == 'yes' for r in cx)}/{len(cx)} with expect.ideal_class, "
          f"{sum(r['class_match_before'] == 'yes' for r in cx)}/{len(cx)} with the labels before roadmap item 5")
    print(f"junctions: {sum(r['status'] == 'exact' for r in jrows)} exact, "
          f"{sum(r['status'] == 'partial' for r in jrows)} partial, "
          f"{sum(r['status'] == 'missed' for r in jrows)} missed of {len(jrows)}")
    print(f"off-target CxSV loci in the slices: {len(off)} ("
          + ", ".join(f"{sum(r['Panel_Status'] == s for r in off)} {s}" for s in ("private", "shared_junctions", "public"))
          + ")")


if __name__ == "__main__":
    main()

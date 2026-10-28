#!/usr/bin/env python3
"""Score the patient workflow's output on the positive control against the truth set.

Junction level (truth_junctions.tsv, one row per novel adjacency):
  exact    a merged call has one breakpoint within --slop of each junction end
           (either order; BND mates parsed from ALT, others use INFO/END);
  partial  a merged call has a breakpoint within --slop of either junction end
           (e.g. an INS or a one-sided call);
  callers  SUPP_VEC of the matching call(s), in SURVIVOR file order:
           cuteSV, SVIM, Sniffles2 (so 001 = Sniffles2 only).
Event level (truth_events.bed + events.yaml expectations):
  loci     patient CxSV loci (private/cxsv_panel_status.tsv) overlapping any
           span of the event, with class, criteria, carriers and panel status;
  verdict  PASS if every expectation in events.yaml holds:
             expect.cxsv true  -> >=1 overlapping CxSV locus, all of them private
             expect.cxsv false -> no overlapping CxSV locus (simple SV controls)
           and >=1 junction detected exactly (or partially for templated INS).
Outside the events: CxSV loci not overlapping any event span are listed as the
host's own CxSV-like calls (expected; a 10-genome panel cannot make most public).

Known failures (--known-failures TSV: event<TAB>reason) turn FAIL into XFAIL; a
listed event that passes is XPASS (update the list). The regression counts FAIL
and XPASS as failures, so it catches new problems without masking old ones.
class_match reports whether the pipeline's class equals expect.ideal_class.

Usage: score_pos_control.py TRUTH_DIR SAMPLE_DIR events.yaml OUT_PREFIX [--slop 1000]
                            [--known-failures known_failures.tsv]
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
    a = ap.parse_args()

    calls = parse_calls(f"{a.sample}/merged_sv.vcf")
    evs = {e["id"]: e for e in yaml.safe_load(open(a.events))["events"]}
    spans = {}
    for line in open(f"{a.truth}/truth_events.bed"):
        c, s, e, eid = line.split("\t")[:4]
        spans.setdefault(eid, []).append((c, int(s), int(e)))
    loci = list(csv.DictReader(open(f"{a.sample}/private/cxsv_panel_status.tsv"), delimiter="\t"))

    jrows, jsum = [], {}
    for j in csv.DictReader(open(f"{a.truth}/truth_junctions.tsv"), delimiter="\t"):
        c1, p1, c2, p2 = j["chrom1"], int(j["pos1"]), j["chrom2"], int(j["pos2"])
        exact = [k for k in calls if (near(k["c1"], k["p1"], c1, p1, a.slop) and near(k["c2"], k["p2"], c2, p2, a.slop))
                 or (near(k["c1"], k["p1"], c2, p2, a.slop) and near(k["c2"], k["p2"], c1, p1, a.slop))]
        part = [k for k in calls if k not in exact and any(near(c, p, cc, pp, a.slop)
                for c, p in ((k["c1"], k["p1"]), (k["c2"], k["p2"])) for cc, pp in ((c1, p1), (c2, p2)))]
        status = "exact" if exact else ("partial" if part else "missed")
        best = exact or part
        jrows.append(dict(event=j["event"], junction=f"{j['left']}|{j['right']}",
                          truth=f"{c1}:{p1}-{c2}:{p2}", status=status,
                          calls=";".join(f"{k['svtype']}:{k['c1']}:{k['p1']}-{k['c2']}:{k['p2']}" for k in best[:3]) or ".",
                          supp_vec=";".join(sorted({k["supp"] for k in best})) or "."))
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
    erows = []
    for eid, ev in evs.items():
        ov = [r for r in loci if overlaps(r, eid)]
        exp = ev.get("expect", {})
        js = jsum.get(eid, {"exact": 0, "partial": 0, "missed": 0})
        n = sum(js.values())
        detected = js["exact"] > 0 or js["partial"] > 0
        if exp.get("cxsv"):
            ok_cx = bool(ov) and all(r["Panel_Status"] == "private" for r in ov)
        else:
            ok_cx = not ov
        erows.append(dict(event=eid, name=ev["name"], gene=ev.get("gene") or "intergenic",
                          expect_cxsv=exp.get("cxsv"), ideal_class=exp.get("ideal_class", "."),
                          junctions=n, exact=js["exact"], partial=js["partial"], missed=js["missed"],
                          cxsv_loci=len(ov),
                          classes=";".join(sorted({r["Cmplx_Class"] for r in ov})) or ".",
                          criteria=";".join(sorted({r["CxSV_Criteria"] for r in ov})) or ".",
                          carriers=";".join(sorted({str(r["Panel_Carriers"]) for r in ov})) or ".",
                          panel_status=";".join(sorted({r["Panel_Status"] for r in ov})) or ".",
                          class_match=("yes" if exp.get("ideal_class") in {r["Cmplx_Class"] for r in ov}
                                       else ("n/a" if not exp.get("cxsv") else "no")),
                          verdict=verdict(detected and ok_cx, eid in known)))
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
                                                         "CxSV_Criteria", "Panel_Carriers", "Panel_Status")}
                                     for r in off])
    for r in erows:
        print(f"{r['event']}\t{r['verdict']}{' (' + known[r['event']] + ')' if r['event'] in known else ''}\tjunctions {r['exact']}+{r['partial']}p/{r['junctions']}\t"
              f"cxsv_loci={r['cxsv_loci']} class={r['classes']} (expected {r['ideal_class']}, match {r['class_match']}) status={r['panel_status']} carriers={r['carriers']}")
    print(f"junctions: {sum(r['status'] == 'exact' for r in jrows)} exact, "
          f"{sum(r['status'] == 'partial' for r in jrows)} partial, "
          f"{sum(r['status'] == 'missed' for r in jrows)} missed of {len(jrows)}")
    print(f"off-target CxSV loci in the slices: {len(off)} "
          f"({sum(r['Panel_Status'] == 'private' for r in off)} private)")


if __name__ == "__main__":
    main()

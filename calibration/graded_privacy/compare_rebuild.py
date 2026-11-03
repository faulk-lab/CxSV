#!/usr/bin/env python3
"""Per-locus comparison of the leave-out shortcut with a true panel rebuild.

Writes one row per locus with both values of every status column the scorer
reads (and the total / tier), and prints a one-line summary: loci compared,
loci with any status difference, and loci whose tier differs."""
import argparse
import csv

STATUS = ["Panel_Status", "Panel_Junctions", "Panel_Junctions_Seen", "Panel_Carriers", "Panel_Callable",
          "Panel_Individuals", "Panel_Best_Junction_Fraction", "Panel_Region_Carriers"]
SCORE = ["Total_Score", "Tier", "panel_privacy_points"]


def rows(path, key):
    with open(path, newline="") as fh:
        return {key(r): r for r in csv.DictReader(fh, delimiter="\t")}


def main():
    ap = argparse.ArgumentParser()
    for k in ("rebuild-status", "rebuild-scores", "shortcut-status", "shortcut-scores", "label", "out"):
        ap.add_argument(f"--{k}", required=True)
    a = ap.parse_args()
    lid = lambda r: f"{r['Chr']}:{r['Start']}-{r['End']}"  # noqa: E731
    rs, ss = rows(a.rebuild_status, lid), rows(a.shortcut_status, lid)
    rc, sc = rows(a.rebuild_scores, lambda r: r["Locus_ID"]), rows(a.shortcut_scores, lambda r: r["Locus_ID"])
    assert set(rs) == set(ss) == set(rc) == set(sc), "locus sets differ"
    n_diff = n_tier = n_status = 0
    diffs_by_col = {c: 0 for c in STATUS + SCORE}
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Locus_ID"] + [f"{c}_{t}" for c in STATUS + SCORE for t in ("rebuild", "shortcut")] + ["differs"])
        for L in sorted(rs):
            vals, d = [], False
            for c in STATUS:
                x, y = rs[L][c], ss[L][c]
                vals += [x, y]
                if x != y:
                    d = True
                    diffs_by_col[c] += 1
            for c in SCORE:
                x, y = rc[L][c], sc[L][c]
                vals += [x, y]
                if x != y:
                    d = True
                    diffs_by_col[c] += 1
            n_diff += d
            n_tier += rc[L]["Tier"] != sc[L]["Tier"]
            n_status += rs[L]["Panel_Status"] != ss[L]["Panel_Status"]
            w.writerow([L] + vals + ["yes" if d else "no"])
    print(f"{a.label}\tloci {len(rs)}\tany_difference {n_diff}\tstatus_differs {n_status}\ttier_differs {n_tier}\t"
          + ";".join(f"{k}:{v}" for k, v in diffs_by_col.items() if v))


if __name__ == "__main__":
    main()

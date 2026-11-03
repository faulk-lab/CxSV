#!/usr/bin/env python3
"""Family groups of the panel individuals (leave-one-family-out, item 10).

Relatedness sources (all under <registry>/sources):
  igsr_3202_samples_ped_population.txt  1kGP 3,202 pedigree: FamilyID, FatherID, MotherID
  igsr_20130606_sample_info.txt         phase-3 sample info: Family ID, Relationship and the
                                        cryptic-relatedness columns (Unexpected Parent/Child,
                                        Siblings, Grandparents, Avuncular, Half Siblings,
                                        Unknown Second Order, Third Order)
Every link (shared family id, parent-child, any listed relative) joins two ids
in a union-find over ALL ids in the sources, so two panel individuals related
through someone outside the panel (e.g. half-siblings via an unsequenced
parent) are in one family. A panel family is a connected component restricted
to the panel individuals.

Output (--out): individual_id, family_id, n_panel_family, panel_family_members
(the OTHER panel individuals of the family, comma list or empty),
relatedness_known (yes when the individual appears in at least one source),
sources, superpopulation, platforms, sex.
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gp_common import read_tsv, write_tsv  # noqa: E402

REL_COLS = ("Unexpected Parent/Child ", "Unexpected Parent/Child", "Non Paternity", "Siblings", "Grandparents",
            "Avuncular", "Half Siblings", "Unknown Second Order", "Third Order")
ID_RE = re.compile(r"\b(?:HG|NA|GM)\d{5}\b")


def norm(x):
    x = x.strip()
    return "NA" + x[2:] if x.startswith("GM") else x


class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--panel-individuals", required=True, help="release panel_individuals.tsv")
    ap.add_argument("--ped", required=True)
    ap.add_argument("--sample-info", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    panel = {r["individual_id"]: r for r in read_tsv(a.panel_individuals)}
    uf = UF()
    seen = defaultdict(set)
    links = defaultdict(int)
    # 1kGP 3,202 pedigree (space separated)
    fam_members = defaultdict(set)
    with open(a.ped) as fh:
        header = fh.readline().split()
        for line in fh:
            f = dict(zip(header, line.split()))
            s = norm(f["SampleID"])
            seen[s].add("ped3202")
            uf.find(s)
            fam_members["ped:" + f["FamilyID"]].add(s)
            for k in ("FatherID", "MotherID"):
                if f[k] not in ("0", ""):
                    uf.union(s, norm(f[k]))
                    seen[norm(f[k])].add("ped3202")
                    links["parent_child"] += 1
    # phase-3 sample info (tab separated)
    with open(a.sample_info, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = norm(r["Sample"])
            seen[s].add("sample_info_20130606")
            uf.find(s)
            if (r.get("Family ID") or "").strip():
                fam_members["info:" + r["Family ID"].strip()].add(s)
            for col in REL_COLS:
                for rel in ID_RE.findall(r.get(col) or ""):
                    uf.union(s, norm(rel))
                    links[col.strip()] += 1
    for fid, members in fam_members.items():
        m = sorted(members)
        for x in m[1:]:
            uf.union(m[0], x)
        if len(m) > 1:
            links["shared_family_id"] += len(m) - 1

    comp = defaultdict(list)
    for i in panel:
        comp[uf.find(i)].append(i)
    rows = []
    for root, members in comp.items():
        members.sort()
        fid = "FAM_" + members[0]
        for i in members:
            p = panel[i]
            rows.append({"individual_id": i, "family_id": fid, "n_panel_family": len(members),
                         "panel_family_members": ",".join(x for x in members if x != i),
                         "relatedness_known": "yes" if seen.get(i) else "no",
                         "sources": ",".join(sorted(seen.get(i, []))) or "none",
                         "superpopulation": p.get("superpopulation", ""), "platforms": p.get("platforms", ""),
                         "sex": p.get("sex", "")})
    rows.sort(key=lambda r: r["individual_id"])
    write_tsv(a.out, rows, ["individual_id", "family_id", "n_panel_family", "panel_family_members",
                            "relatedness_known", "sources", "superpopulation", "platforms", "sex"])
    sizes = defaultdict(int)
    for members in comp.values():
        sizes[len(members)] += 1
    unknown = [r["individual_id"] for r in rows if r["relatedness_known"] == "no"]
    print(f"{len(panel)} panel individuals in {len(comp)} panel families; family sizes "
          + ", ".join(f"{k}:{v}" for k, v in sorted(sizes.items()))
          + f"; relatedness unknown (in no source) for {len(unknown)}: {','.join(unknown[:30])}"
          + f"; links {dict(links)}", file=sys.stderr)


if __name__ == "__main__":
    main()

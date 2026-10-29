"""Panel junction sites and breakpoint-level matching of patient CxSV loci.

Panel side (`cluster_sites`): every background sample's junctions (the
collapsed per-sample tables written by scripts/extract_junctions.py) are
clustered across individuals into *sites*. Clustering is type-aware (svtype
and both contigs must agree) and leader-based: observations are visited in
position order and join the nearest open site whose leader (first member) has
both ends within `tolerance` bp, otherwise they found a new site. A site
therefore never spans more than `tolerance` bp at either end, and two sites
cannot chain into one. INS additionally need comparable lengths
(shorter/longer >= `ins_len_ratio`).

Each site is written once in the shared junction-site schema
(docs/JUNCTIONS.md; SITE_COLUMNS + af_<group>) with its representative
coordinates (lower median of its members), and once in the carrier table
(CARRIER_COLUMNS) with its coordinate envelope and carrier individuals, which is
what the patient matcher reads.

Patient side (`SiteIndex`, `classify_locus`): a patient junction matches a site
when svtype and contigs agree and each end lies within `tolerance` bp of the
site's envelope. A panel individual carries a patient CxSV locus when it
carries at least `min_junction_fraction` of the locus's junctions. Statuses:

  public            >= min_carriers individuals carry the rearrangement
  shared_junctions  some of the locus's junctions occur in the panel, but no
                    (or too few) individuals carry the rearrangement — e.g. a
                    different rearrangement, or a common simple SV, at a hotspot
  private           none of the locus's junctions occurs in the panel
"""
from __future__ import annotations

import bisect
import csv
from collections import defaultdict

SITE_COLUMNS = ("chrom1", "pos1", "chrom2", "pos2", "svtype", "svlen", "ac", "an", "af", "source", "version")
CARRIER_COLUMNS = ("chrom1", "pos1", "chrom2", "pos2", "svtype", "svlen",
                   "min_pos1", "max_pos1", "min_pos2", "max_pos2", "min_svlen", "max_svlen",
                   "n_observations", "n_carriers", "carriers")
STATUSES = ("private", "shared_junctions", "public")


# ─── schema conventions shared with the catalogue track (roadmap item 4,
# acquisition/catalogues/junction_sites.py on branch item4-catalogue-acquisition)
#   svlen unsigned; BND svlen '.' (these junction BNDs are breakend pairs, never
#   the outer span of a complex SV); af from the source when given, else ac/an;
#   sites of non-BND SVs shorter than the pipeline's sv_params.min_sv_len (50 bp)
#   are dropped; orientation and sort order by hs1 .fai contig order.
def svlen_field(svtype, svlen):
    return "." if svtype == "BND" else str(int(svlen))


def parse_svlen(x):
    return 0 if x in (".", "", None) else int(x)


def too_short(svtype, svlen, min_svlen):
    return svtype != "BND" and int(svlen) < min_svlen


# ─── individuals and ploidy ──────────────────────────────────────────────────
def ploidy(chrom, sex):
    """Copies of a contig in one individual (PAR not modelled: chrX/chrY are
    treated as hemizygous in males along their whole length)."""
    c = chrom.lower().removeprefix("chr")
    if c == "y":
        return 1 if sex in ("male", "unknown", "") else 0
    if c == "x":
        return 1 if sex == "male" else 2
    if c in ("m", "mt"):
        return 1
    return 2


def site_ploidy(chrom1, chrom2, sex):
    return min(ploidy(chrom1, sex), ploidy(chrom2, sex))


def groups_of(ind):
    """Frequency groups an individual counts in: each platform it was sequenced
    on and its superpopulation (when labelled)."""
    g = [f"platform:{p}" for p in sorted(ind["platforms"])]
    sp = ind.get("superpopulation") or ""
    if sp and sp.lower() not in ("na", "unknown", "."):
        g.append(f"superpopulation:{sp}")
    return g


def group_column(group):
    return "af_" + group.split(":", 1)[1]


# ─── clustering ──────────────────────────────────────────────────────────────
def compatible(lead, obs, tolerance, ins_len_ratio):
    if abs(lead["pos1"] - obs["pos1"]) > tolerance or abs(lead["pos2"] - obs["pos2"]) > tolerance:
        return False
    if lead["svtype"] == "INS":
        lo, hi = sorted((lead["svlen"], obs["svlen"]))
        if hi > 0 and lo / hi < ins_len_ratio:
            return False
    return True


def _lower_median(vals):
    s = sorted(vals)
    return s[(len(s) - 1) // 2]


def cluster_sites(observations, tolerance, ins_len_ratio, order):
    """Cluster junction observations (dicts: svtype chrom1 pos1 chrom2 pos2 svlen
    individual dosage ...) into sites. Deterministic for a given input set."""
    by_key = defaultdict(list)
    for o in observations:
        by_key[(o["svtype"], o["chrom1"], o["chrom2"])].append(o)
    sites = []
    for key in sorted(by_key, key=lambda k: (order.get(k[1], 1 << 30), order.get(k[2], 1 << 30), k[0])):
        obs = sorted(by_key[key], key=lambda o: (o["pos1"], o["pos2"], o["svlen"], o["individual"],
                                                 o.get("sample", "")))
        active = []                                  # sites whose leader is within tolerance of pos1
        for o in obs:
            active = [s for s in active if o["pos1"] - s["leader"]["pos1"] <= tolerance]
            best, best_d = None, None
            for s in active:
                if compatible(s["leader"], o, tolerance, ins_len_ratio):
                    d = abs(s["leader"]["pos1"] - o["pos1"]) + abs(s["leader"]["pos2"] - o["pos2"])
                    if best is None or d < best_d:
                        best, best_d = s, d
            if best is None:
                best = {"leader": o, "members": []}
                active.append(best)
                sites.append(best)
            best["members"].append(o)
    out = []
    for s in sites:
        m = s["members"]
        L = s["leader"]
        out.append({
            "chrom1": L["chrom1"], "chrom2": L["chrom2"], "svtype": L["svtype"],
            "pos1": _lower_median([x["pos1"] for x in m]), "pos2": _lower_median([x["pos2"] for x in m]),
            "svlen": _lower_median([x["svlen"] for x in m]),
            "min_pos1": min(x["pos1"] for x in m), "max_pos1": max(x["pos1"] for x in m),
            "min_pos2": min(x["pos2"] for x in m), "max_pos2": max(x["pos2"] for x in m),
            "min_svlen": min(x["svlen"] for x in m), "max_svlen": max(x["svlen"] for x in m),
            "n_observations": len(m), "members": m,
        })
    # representative coordinates are unique per (type, contigs): merge exact collisions
    merged = {}
    for s in out:
        k = (s["svtype"], s["chrom1"], s["pos1"], s["chrom2"], s["pos2"], s["svlen"])
        if k in merged:
            t = merged[k]
            for f, fn in (("min_pos1", min), ("max_pos1", max), ("min_pos2", min), ("max_pos2", max),
                          ("min_svlen", min), ("max_svlen", max)):
                t[f] = fn(t[f], s[f])
            t["members"] += s["members"]
            t["n_observations"] += s["n_observations"]
        else:
            merged[k] = s
    res = list(merged.values())
    res.sort(key=lambda s: (order.get(s["chrom1"], 1 << 30), s["pos1"], order.get(s["chrom2"], 1 << 30),
                            s["pos2"], s["svtype"], s["svlen"]))
    return res


def site_counts(site, individuals, groups):
    """Carriers, allele counts and allele numbers of one site, overall and per group.

    ac sums each carrier's dosage (1 het, 2 hom; capped by ploidy), an sums the
    ploidy of every panel individual, so af is an allele frequency. Carriers are
    distinct individuals (one person sequenced twice counts once)."""
    dose = {}
    for m in site["members"]:
        dose[m["individual"]] = max(dose.get(m["individual"], 0), int(m.get("dosage", 1)))
    res = {"carriers": sorted(dose)}
    tot = {"ac": 0, "an": 0}
    per = {g: {"ac": 0, "an": 0, "carriers": 0} for g in groups}
    for ind_id, ind in individuals.items():
        p = site_ploidy(site["chrom1"], site["chrom2"], ind.get("sex", ""))
        a = min(dose.get(ind_id, 0), p)
        if ind_id in dose and p == 0:
            a = dose[ind_id]          # observed where ploidy says none (e.g. sex mislabel): count it
            p = a
        tot["ac"] += a
        tot["an"] += p
        for g in groups_of(ind):
            if g in per:
                per[g]["ac"] += a
                per[g]["an"] += p
                per[g]["carriers"] += ind_id in dose
    res.update(tot)
    res["groups"] = per
    return res


def fmt_af(ac, an):
    return f"{ac / an:.6f}" if an else "."


# ─── matching ────────────────────────────────────────────────────────────────
def read_carriers(path):
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k in ("pos1", "pos2", "min_pos1", "max_pos1", "min_pos2", "max_pos2",
                      "min_svlen", "max_svlen", "n_observations", "n_carriers"):
                r[k] = int(r[k])
            r["svlen"] = parse_svlen(r["svlen"])
            r["carrier_set"] = set(r["carriers"].split(",")) if r["carriers"] else set()
            rows.append(r)
    return rows


class SiteIndex:
    """Panel sites by (svtype, chrom1, chrom2), sorted by envelope start."""

    def __init__(self, sites):
        self.by_key = defaultdict(list)
        for s in sites:
            self.by_key[(s["svtype"], s["chrom1"], s["chrom2"])].append(s)
        self.starts, self.maxspan = {}, {}
        for k, v in self.by_key.items():
            v.sort(key=lambda s: (s["min_pos1"], s["min_pos2"]))
            self.starts[k] = [s["min_pos1"] for s in v]
            self.maxspan[k] = max(s["max_pos1"] - s["min_pos1"] for s in v)

    def query(self, j, tolerance, ins_len_ratio):
        """Sites matching junction j (dict with svtype chrom1 pos1 chrom2 pos2 svlen)."""
        k = (j["svtype"], j["chrom1"], j["chrom2"])
        if k not in self.by_key:
            return []
        v, starts = self.by_key[k], self.starts[k]
        lo = bisect.bisect_left(starts, j["pos1"] - tolerance - self.maxspan[k])
        hi = bisect.bisect_right(starts, j["pos1"] + tolerance)
        hits = []
        for s in v[lo:hi]:
            if not (s["min_pos1"] - tolerance <= j["pos1"] <= s["max_pos1"] + tolerance):
                continue
            if not (s["min_pos2"] - tolerance <= j["pos2"] <= s["max_pos2"] + tolerance):
                continue
            if j["svtype"] == "INS":
                near = min(max(j["svlen"], s["min_svlen"]), s["max_svlen"])
                lo_l, hi_l = sorted((j["svlen"], near))
                if hi_l > 0 and lo_l / hi_l < ins_len_ratio:
                    continue
            hits.append(s)
        return hits


def junction_carriers(j, index, tolerance, ins_len_ratio):
    """(set of panel individuals carrying junction j, matched sites)."""
    hits = index.query(j, tolerance, ins_len_ratio)
    inds = set()
    for s in hits:
        inds |= s["carrier_set"]
    return inds, hits


def classify_locus(locus_junctions, index, tolerance, ins_len_ratio, min_junction_fraction, min_carriers):
    """Breakpoint-level panel status of one patient CxSV locus.

    Returns n_junctions, n_seen (junctions with >= 1 panel carrier), carriers
    (individuals carrying >= min_junction_fraction of the junctions), best
    fraction over panel individuals, per-junction carrier sets, and status."""
    per_j = [junction_carriers(j, index, tolerance, ins_len_ratio)[0] for j in locus_junctions]
    n = len(per_j)
    count = defaultdict(int)
    for inds in per_j:
        for i in inds:
            count[i] += 1
    frac = {i: c / n for i, c in count.items()} if n else {}
    carriers = sorted(i for i, f in frac.items() if f >= min_junction_fraction - 1e-9)
    n_seen = sum(1 for inds in per_j if inds)
    if n and len(carriers) >= min_carriers:
        status = "public"
    elif n_seen:
        status = "shared_junctions"
    else:
        status = "private"
    return {"n_junctions": n, "n_seen": n_seen, "carriers": carriers,
            "best_fraction": max(frac.values()) if frac else 0.0, "per_junction": per_j,
            "status": status}


# ─── point sites: external catalogues and family members (roadmap items 4, 5) ─
# The same matcher (SiteIndex.query) is used for the panel, the external SV
# catalogues and a patient's parents. A catalogue site or a parent's junction
# is a single point: its envelope is the point itself, and `tolerance` is
# applied around it exactly as around a panel site's envelope.
def point_site(svtype, chrom1, pos1, chrom2, pos2, svlen, **extra):
    s = {"svtype": svtype, "chrom1": chrom1, "pos1": pos1, "chrom2": chrom2, "pos2": pos2, "svlen": svlen,
         "min_pos1": pos1, "max_pos1": pos1, "min_pos2": pos2, "max_pos2": pos2,
         "min_svlen": svlen, "max_svlen": svlen, "carrier_set": set()}
    s.update(extra)
    return s


def _near_any(sorted_pos, p, tol):
    i = bisect.bisect_left(sorted_pos, p - tol)
    return i < len(sorted_pos) and sorted_pos[i] <= p + tol


def read_catalogue_sites(path, query_pos1=None, tolerance=0):
    """Rows of a shared-schema junction-site table (plain or gzip) as point sites.

    With query_pos1 ({chrom1: sorted positions}), only rows whose pos1 lies
    within `tolerance` of a query position are kept (a 1.2 M-row catalogue is
    streamed once, not held in memory). Each site keeps af (float or None) and
    the per-group af_<group> values (dict)."""
    import gzip
    opener = gzip.open if path.endswith((".gz", ".bgz")) else open
    out = []
    with opener(path, "rt", newline="") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        groups = [c for c in header if c.startswith("af_")]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            c1, p1 = f[ix["chrom1"]], int(f[ix["pos1"]])
            if query_pos1 is not None:
                qp = query_pos1.get(c1)
                if not qp or not _near_any(qp, p1, tolerance):
                    continue
            af = f[ix["af"]]
            out.append(point_site(f[ix["svtype"]], c1, p1, f[ix["chrom2"]], int(f[ix["pos2"]]),
                                  parse_svlen(f[ix["svlen"]]),
                                  af=None if af in (".", "") else float(af),
                                  group_af={g[3:]: float(f[ix[g]]) for g in groups if f[ix[g]] not in (".", "")}))
    return out


def junction_point_sites(junction_rows, carrier=None):
    """A sample's junction-table rows (cxsv_junctions.read_table) as point sites;
    `carrier` (e.g. the parent's id) is put in each site's carrier_set."""
    return [point_site(j["svtype"], j["chrom1"], j["pos1"], j["chrom2"], j["pos2"], j["svlen"],
                       carrier_set={carrier} if carrier else set()) for j in junction_rows]


def match_fraction(locus_junctions, index, tolerance, ins_len_ratio):
    """(n matched, n, matched sites per junction) of a locus against one point-site index."""
    per = [index.query(j, tolerance, ins_len_ratio) for j in locus_junctions]
    return sum(1 for h in per if h), len(per), per

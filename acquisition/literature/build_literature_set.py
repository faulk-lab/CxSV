#!/usr/bin/env python3
"""Build the curated literature set of pathogenic germline complex SVs on hs1 (roadmap item 11).

    python3 acquisition/literature/build_literature_set.py \\
        --sources /mnt/sargasso1/CxSV/resources/literature/_downloads/20261104 \\
        --chain-hg19 .../hg19ToHg38.over.chain.gz --chain-hg38 .../hg38ToHs1.over.chain.gz \\
        --fai /mnt/sargasso1/CxSV/refs/hs1.fa.fai --genes /mnt/sargasso1/CxSV/resources/genes.final.bed \\
        --version lit-v1.0 --out /mnt/sargasso1/CxSV/resources/literature/lit-v1.0

Inputs: the parsed source files listed in sources.tsv (role 'parsed', in --sources)
and the curated tables in curated/ (facts transcribed from tables, cited row by row).
Every breakpoint is lifted on its own with the item-4 code
(acquisition/catalogues/liftover.lift_points, UCSC liftOver, default -minMatch):
GRCh37 sources hg19 -> hg38 -> hs1 (two chains), GRCh38 sources hg38 -> hs1.
The base next to each breakpoint (pos + 1) is lifted as well: when it lands
before the breakpoint the chain block is reversed, and the side is flipped.

Per junction: lifted (both ends on one hs1 contig for DEL/DUP/INV, any contigs
for BND), failed_unmapped (an end does not lift), failed_split (a DEL/DUP/INV
whose ends land on two contigs); span_changed (> 50% span change) is reported.

Outputs in --out:
  literature_events.tsv                 one row per event: source, citation, case, phenotype,
                                        inheritance, author type and classification, junctions
                                        (source / lifted), precision, hs1 region(s), genes at breakends
  literature_junctions.tsv              every junction: hs1 ends with sides and svtype, the source-build
                                        ends, lift status, precision, evidence
  literature_cn_segments.tsv            curated copy-number segments on hs1 (scoring test only)
  literature_derivatives.tsv            derivative chromosomes as fragment chains (Collins et al.), hs1 + source
  literature_set.<version>.junction_sites.tsv.gz (+ .tbi when tabix is on PATH)
                                        lifted junctions in the shared junction-site schema
                                        (docs/JUNCTIONS.md), source=cxsv_literature
  liftover_report.tsv                   per source: junctions in / lifted / failed / split / span changed
  literature_set.json                   version, inputs with sha256, chains with sha256, counts
  MANIFEST.tsv                          file bytes sha256 for every output
  README.md                             citations and curation rules (copied from this directory)
"""
import argparse
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "catalogues"))

import literature_set as L        # noqa: E402
import junction_sites as js       # noqa: E402
import liftover                   # noqa: E402

JUNCTION_COLUMNS = ["event_id", "junction_id", "svtype", "chrom1", "pos1", "side1", "chrom2", "pos2", "side2",
                    "svlen", "precision", "lift_status", "span_changed", "sides_flipped_by_lift",
                    "src_build", "src_chrom1", "src_pos1", "src_side1", "src_chrom2", "src_pos2", "src_side2",
                    "src_svtype", "hg38_chrom1", "hg38_pos1", "hg38_chrom2", "hg38_pos2", "type_basis", "evidence"]
EVENT_COLUMNS = ["event_id", "source_id", "citation", "doi", "pmid", "case_id", "source_table", "build",
                 "phenotype", "genes_reported", "inheritance", "author_event_type", "author_classification",
                 "known_pathogenic", "carrier_sex", "zygosity", "n_junctions_source", "n_junctions_lifted",
                 "n_breakends_hs1", "n_chromosomes", "precision", "lift_status", "region_hs1", "span_hs1",
                 "genes_at_breakends", "flags", "notes"]
DERIV_COLUMNS = ["event_id", "derivative", "order", "fragment", "chrom", "start", "end", "strand_hs1", "lift_status",
                 "src_build", "src_chrom", "src_start", "src_end", "strand"]
CN_COLUMNS = ["event_id", "chrom", "start", "end", "delta", "lift_status", "src_chrom", "src_start", "src_end",
              "evidence"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def contig_order(fai):
    return {line.split("\t")[0]: i for i, line in enumerate(open(fai)) if line.strip()}


# ─── collection ──────────────────────────────────────────────────────────────
def collect(src_dir, cur_dir, order):
    """All source-build junctions, CN segments and events. Returns (events, junctions, cn, notes)."""
    sources = {}
    for r in L.read_tsv(os.path.join(HERE, "sources.tsv")):
        sources.setdefault(r["source_id"], r)
    events = {r["event_id"]: dict(r) for r in L.read_tsv(os.path.join(cur_dir, "events.tsv"))}
    junctions, cn, notes, derivs = [], [], defaultdict(list), []

    # Nazaryan-Petersen 2018: S2 Table + curated junction types
    s2 = L.read_xlsx(os.path.join(src_dir, "np2018_s009.xlsx"))
    types = {(r["case"], r["junction"]): r for r in L.read_tsv(os.path.join(cur_dir, "nazaryan2018_junction_types.tsv"))}
    jn, skipped = L.parse_nazaryan2018(next(iter(s2.values())), types, order)
    junctions += jn
    for s in skipped:
        notes[s["event_id"]].append(f"J{s['junction']} not representable ({s['reason']})")

    # Collins 2017: derivative chains
    for case, rows in L.read_xlsx(os.path.join(src_dir, "13059_2017_1158_MOESM3_ESM.xlsx")).items():
        e = f"COL17_{case}"
        ders, cnvs = L.parse_collins_derivatives(rows)
        ders, moved = L.fix_derivative_boundaries(ders)
        for m in moved:
            notes[e].append(f"fragment {m} (derivative label offset in the table)")
        k = 0
        for der, frags in ders.items():
            for i, (lab, c, st, en, strand) in enumerate(frags, 1):
                derivs.append(dict(event_id=e, derivative=der, order=i, fragment=lab, src_chrom=c, src_start=st,
                                   src_end=en, strand=strand, src_build="GRCh37"))
            for la, lb, c1, p1, s1, c2, p2, s2 in L.derivative_junctions(frags):
                k += 1
                junctions.append(L.make_junction(e, k, c1, p1, s1, c2, p2, s2, "approximate", order,
                                                 type_basis="derivative", evidence=f"{der} {la}|{lb}"))
        for x in cnvs:
            if "copies" in x:
                delta = round(x["copies"] - 2, 2)
            else:
                delta = -1 if x["class"].lower().startswith("loss") else 1
            cn.append(dict(event_id=e, chrom=x["chrom"], start=x["start"], end=x["end"], delta=delta,
                           evidence=f"Additional file 3 {x['label']} {x['class']}"))

    # Grochowski 2024: DUP-TRP/INV-DUP segments
    gj, gc, gf = L.parse_grochowski2024(L.read_tsv(os.path.join(cur_dir, "grochowski2024_table_s1.tsv")), order)
    for j in gj:
        j["type_basis"] = "model"
    junctions += gj
    cn += [dict(x, evidence="Table S1") for x in gc]
    for f in gf:
        notes[f["event_id"]].append(f["flag"])

    # Hand-transcribed junctions (Eisfeldt 2024, Bahrambeigi 2019)
    for r in L.read_tsv(os.path.join(cur_dir, "junctions_manual.tsv")):
        junctions.append(L.make_junction(r["event_id"], r["junction"], L.ucsc(r["chrom1"]), int(r["pos1"]), r["side1"],
                                         L.ucsc(r["chrom2"]), int(r["pos2"]), r["side2"], r["precision"], order,
                                         type_basis="table", evidence=r["evidence"]))
    for r in L.read_tsv(os.path.join(cur_dir, "cn_segments.tsv")):
        cn.append(dict(event_id=r["event_id"], chrom=L.ucsc(r["chrom"]), start=int(r["start"]), end=int(r["end"]),
                       delta=float(r["delta"]), evidence=r["evidence"]))

    # ClinVar (dbVar nstd102), region level
    cv_events, cv_j, cv_counts = L.parse_clinvar_nstd102(os.path.join(src_dir, "nstd102.GRCh38.variant_call.tsv.gz"),
                                                        order)
    for e in cv_events:
        events[e["event_id"]] = e
    for j in cv_j:
        j["type_basis"] = "region"
    junctions += cv_j

    build_of = {e: ev["build"] for e, ev in events.items()}
    for j in junctions:
        j["src_build"] = build_of[j["event_id"]]
    for x in cn:
        x["src_build"] = build_of[x["event_id"]]
    for e, ev in events.items():
        s = sources[ev["source_id"]]
        ev.update(citation=s["citation"], doi=s["doi"], pmid=s["pmid"])
        cls = ev.get("author_classification", "")
        ev["known_pathogenic"] = "yes" if cls.lower().startswith(("pathogenic", "likely pathogenic")) else "no"
    return events, junctions, cn, notes, cv_counts, derivs


# ─── liftover ────────────────────────────────────────────────────────────────
def lift_chain(points_by_build, chains, workdir):
    """{build: set((chrom, pos))} -> {build: {(chrom, pos): (hs1 chrom, pos) | None}}, plus hg38 positions.

    GRCh37 points are lifted hg19 -> hg38, then the hg38 result -> hs1."""
    out, hg38 = {}, {}
    for build, pts in points_by_build.items():
        pts = set(pts)
        if build == "GRCh37":
            m19 = liftover.lift_points(pts, chains["hg19"], workdir=workdir)
            mid = {p: m19.get(p) for p in pts}
        elif build == "GRCh38":
            mid = {p: p for p in pts}
        else:
            raise ValueError(f"no chain for build {build}")
        m38 = liftover.lift_points({v for v in mid.values() if v}, chains["hg38"], workdir=workdir)
        out[build] = {p: (m38.get(mid[p]) if mid[p] else None) for p in pts}
        hg38[build] = mid
    return out, hg38


def lifted_end(build, chrom, pos, mapped):
    """hs1 (chrom, pos, reversed?) of one breakpoint, or None."""
    a = mapped[build].get((chrom, pos))
    b = mapped[build].get((chrom, pos + 1))
    if not a or js.normalise_contig(a[0]) is None:
        return None
    rev = bool(b and b[0] == a[0] and b[1] < a[1])
    return js.normalise_contig(a[0]), a[1], rev


def lift_all(events, junctions, cn, chains, order, workdir, derivs=()):
    pts = defaultdict(set)
    for f in derivs:
        for p in (f["src_start"], f["src_end"]):
            if p not in L.TELOMERES:
                pts[f["src_build"]].add((f["src_chrom"], int(p)))
    for j in junctions:
        for c, p in ((j["chrom1"], j["pos1"]), (j["chrom2"], j["pos2"])):
            pts[j["src_build"]].update({(c, p), (c, p + 1)})
    for x in cn:
        pts[x["src_build"]].update({(x["chrom"], x["start"]), (x["chrom"], x["end"])})
    for e in events.values():
        if e.get("region_chrom"):
            pts[e["build"]].update({(e["region_chrom"], int(e["region_start"])), (e["region_chrom"], int(e["region_end"]))})
    mapped, hg38 = lift_chain(pts, chains, workdir)

    lifted = []
    for j in junctions:
        b = j["src_build"]
        out = {k: j.get(k) for k in ("event_id", "junction_id", "precision", "type_basis", "evidence")}
        out.update(src_build=b, src_chrom1=j["chrom1"], src_pos1=j["pos1"], src_side1=j["side1"],
                   src_chrom2=j["chrom2"], src_pos2=j["pos2"], src_side2=j["side2"], src_svtype=j["svtype"])
        h1, h2 = hg38[b].get((j["chrom1"], j["pos1"])), hg38[b].get((j["chrom2"], j["pos2"]))
        out.update(hg38_chrom1=h1[0] if h1 else None, hg38_pos1=h1[1] if h1 else None,
                   hg38_chrom2=h2[0] if h2 else None, hg38_pos2=h2[1] if h2 else None)
        e1 = lifted_end(b, j["chrom1"], j["pos1"], mapped)
        e2 = lifted_end(b, j["chrom2"], j["pos2"], mapped)
        if not e1 or not e2:
            out["lift_status"] = "failed_unmapped"
            lifted.append(out)
            continue
        s1 = L.flip(j["side1"]) if e1[2] else j["side1"]
        s2 = L.flip(j["side2"]) if e2[2] else j["side2"]
        if j["svtype"] != "BND" and e1[0] != e2[0]:
            out["lift_status"] = "failed_split"
            lifted.append(out)
            continue
        n = L.make_junction(j["event_id"], 0, e1[0], e1[1], s1, e2[0], e2[1], s2, j["precision"], order)
        out.update({k: n[k] for k in ("chrom1", "pos1", "side1", "chrom2", "pos2", "side2", "svtype", "svlen")})
        out["sides_flipped_by_lift"] = "yes" if (e1[2] or e2[2]) else "no"
        src_span = abs(j["pos2"] - j["pos1"])
        out["span_changed"] = "yes" if (n["svtype"] != "BND" and src_span and
                                        abs(abs(n["pos2"] - n["pos1"]) - src_span) > 0.5 * src_span) else "no"
        out["lift_status"] = "lifted"
        lifted.append(out)

    cn_out = []
    for x in cn:
        b = x["src_build"]
        a, z = mapped[b].get((x["chrom"], x["start"])), mapped[b].get((x["chrom"], x["end"]))
        row = dict(event_id=x["event_id"], delta=x["delta"], src_chrom=x["chrom"], src_start=x["start"],
                   src_end=x["end"], evidence=x["evidence"])
        if a and z and a[0] == z[0] and js.normalise_contig(a[0]):
            s, e = sorted((a[1], z[1]))
            ratio = (e - s) / max(1, x["end"] - x["start"])
            row.update(chrom=a[0], start=s, end=e, lift_status="lifted" if 0.5 <= ratio <= 2 else "span_changed")
        else:
            row.update(lift_status="failed")
        cn_out.append(row)

    for e in events.values():
        if e.get("region_chrom"):
            a = mapped[e["build"]].get((e["region_chrom"], int(e["region_start"])))
            z = mapped[e["build"]].get((e["region_chrom"], int(e["region_end"])))
            e["_region"] = (a[0], *sorted((a[1], z[1]))) if (a and z and a[0] == z[0]) else None
    for f in derivs:
        b = f["src_build"]
        ends = []
        for p in (f["src_start"], f["src_end"]):
            if p in L.TELOMERES:
                ends.append((f["src_chrom"], p))
            else:
                m = mapped[b].get((f["src_chrom"], int(p)))
                ends.append((js.normalise_contig(m[0]), m[1]) if m else None)
        if None in ends or ends[0][0] != ends[1][0]:
            f.update(lift_status="failed")
            continue
        (c, s0), (_, e0) = ends
        nums = [x for x in (s0, e0) if x not in L.TELOMERES]
        rev = len(nums) == 2 and s0 > e0
        f.update(chrom=c, start=e0 if rev else s0, end=s0 if rev else e0,
                 strand_hs1=L.flip(f["strand"]) if rev else f["strand"], lift_status="lifted")
    return lifted, cn_out


# ─── annotation ──────────────────────────────────────────────────────────────
def load_genes(path):
    by = defaultdict(list)
    if path:
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                by[f[0]].append((int(f[1]), int(f[2]), f[3]))
    return by


def genes_at(genes, chrom, pos):
    return {g for s, e, g in genes.get(chrom, ()) if s < pos <= e}


def summarise_events(events, lifted, cn_out, genes, notes):
    by = defaultdict(list)
    for j in lifted:
        by[j["event_id"]].append(j)
    for e, ev in events.items():
        js_ = by.get(e, [])
        ok = [j for j in js_ if j["lift_status"] == "lifted"]
        ev["n_junctions_source"] = len(js_)
        ev["n_junctions_lifted"] = len(ok)
        precs = {j["precision"] for j in js_}
        ev["precision"] = "region" if not js_ else ("bp" if precs == {"bp"} else
                                                    "region" if precs == {"region"} else
                                                    "approximate" if "bp" not in precs else "mixed")
        if not js_:
            ev["lift_status"] = "region_lifted" if ev.get("_region") else ("region_failed" if ev.get("region_chrom")
                                                                            else "no_junctions")
        else:
            ev["lift_status"] = "lifted" if len(ok) == len(js_) else ("partial" if ok else "failed")
        span = defaultdict(list)
        for j in ok:
            span[j["chrom1"]].append(j["pos1"])
            span[j["chrom2"]].append(j["pos2"])
        for x in cn_out:
            if x["event_id"] == e and x.get("lift_status") == "lifted":
                span[x["chrom"]] += [x["start"], x["end"]]
        if ev.get("_region"):
            c, s, t = ev["_region"]
            span[c] += [s, t]
        ev["region_hs1"] = ";".join(f"{c}:{min(p)}-{max(p)}" for c, p in sorted(span.items())) or None
        ev["span_hs1"] = sum(max(p) - min(p) for p in span.values()) if span else None
        ev["n_chromosomes"] = len(span) or None
        ev["n_breakends_hs1"] = L.breakends(ok) if ok else 0
        gs = set()
        for j in ok:
            gs |= genes_at(genes, j["chrom1"], j["pos1"]) | genes_at(genes, j["chrom2"], j["pos2"])
        ev["genes_at_breakends"] = ",".join(sorted(gs)) or None
        flags = list(notes.get(e, []))
        if any(j.get("span_changed") == "yes" for j in ok):
            flags.append("span changed > 50% by liftover")
        if any(j.get("sides_flipped_by_lift") == "yes" for j in ok):
            flags.append("chain block reversed at a breakpoint (sides flipped)")
        ev["flags"] = "; ".join(flags) or None


def liftover_report(events, lifted):
    rows = defaultdict(Counter)
    for j in lifted:
        s = events[j["event_id"]]["source_id"]
        rows[s]["junctions"] += 1
        rows[s][j["lift_status"]] += 1
        rows[s]["span_changed"] += j.get("span_changed") == "yes"
        rows[s]["sides_flipped"] += j.get("sides_flipped_by_lift") == "yes"
    for e in events.values():
        rows[e["source_id"]]["events"] += 1
        rows[e["source_id"]][f"events_{e['lift_status']}"] += 1
    cols = ["source_id", "events", "events_lifted", "events_partial", "events_failed", "events_region_lifted",
            "events_region_failed", "events_no_junctions", "junctions", "lifted", "failed_unmapped", "failed_split",
            "span_changed", "sides_flipped"]
    return [dict({"source_id": s}, **{c: v.get(c, 0) for c in cols[1:]}) for s, v in sorted(rows.items())], cols


def write_sites(path, lifted, version, order_names):
    order = js.ContigOrder(order_names)
    sites, seen = [], set()
    for j in lifted:
        if j["lift_status"] != "lifted":
            continue
        k = (j["chrom1"], j["pos1"], j["chrom2"], j["pos2"], j["svtype"])
        if k in seen:
            continue
        seen.add(k)
        sites.append(js.make_site(j["chrom1"], j["pos1"], j["chrom2"], j["pos2"], j["svtype"],
                                  None if j["svtype"] == "BND" else j["svlen"]))
    tmp = path[:-3]
    with open(tmp, "w", newline="") as fh:
        counts = js.write_sites(fh, sites, "cxsv_literature", version, [], order)
    if shutil.which("bgzip"):
        subprocess.run(["bgzip", "-f", tmp], check=True)
        if shutil.which("tabix"):
            subprocess.run(["tabix", "-f", "-s1", "-b2", "-e2", "-S1", path], check=True)
    else:
        with open(tmp, "rb") as fi, gzip.open(path, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        os.remove(tmp)
    return counts


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sources", required=True, help="directory with the parsed source files (sources.tsv)")
    ap.add_argument("--curated", default=os.path.join(HERE, "curated"))
    ap.add_argument("--chain-hg19", required=True, help="UCSC hg19ToHg38.over.chain.gz")
    ap.add_argument("--chain-hg38", required=True, help="UCSC hg38ToHs1.over.chain.gz (the item-4 chain)")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--genes", default=None, help="genes.final.bed (hs1) for genes_at_breakends")
    ap.add_argument("--version", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--code-commit", default="unknown")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    order = contig_order(a.fai)
    for r in L.read_tsv(os.path.join(HERE, "sources.tsv")):
        if r["role"] == "parsed":
            got = sha256(os.path.join(a.sources, r["file"]))
            if got != r["sha256"]:
                sys.exit(f"{r['file']}: sha256 {got} != sources.tsv {r['sha256']}")
    events, junctions, cn, notes, cv_counts, derivs = collect(a.sources, a.curated, order)
    chains = {"hg19": a.chain_hg19, "hg38": a.chain_hg38}
    lifted, cn_out = lift_all(events, junctions, cn, chains, order, a.out, derivs)
    summarise_events(events, lifted, cn_out, load_genes(a.genes), notes)

    o = lambda n: os.path.join(a.out, n)  # noqa: E731
    L.write_tsv(o("literature_junctions.tsv"), lifted, JUNCTION_COLUMNS,
                f"CxSV literature set {a.version}: junctions (hs1 + source build). docs/LITERATURE_SET.md")
    L.write_tsv(o("literature_events.tsv"), sorted(events.values(), key=lambda e: (e["source_id"], e["event_id"])),
                EVENT_COLUMNS, f"CxSV literature set {a.version}: events. docs/LITERATURE_SET.md")
    L.write_tsv(o("literature_cn_segments.tsv"), cn_out, CN_COLUMNS,
                f"CxSV literature set {a.version}: curated copy-number changes (scoring test only)")
    L.write_tsv(o("literature_derivatives.tsv"), derivs, DERIV_COLUMNS,
                f"CxSV literature set {a.version}: derivative chromosomes (fragment chains) on hs1 and GRCh37")
    rep, cols = liftover_report(events, lifted)
    L.write_tsv(o("liftover_report.tsv"), rep, cols,
                f"Liftover per source; chains hg19ToHg38 sha256 {sha256(a.chain_hg19)[:16]}, "
                f"hg38ToHs1 sha256 {sha256(a.chain_hg38)[:16]}")
    sites_name = f"literature_set.{a.version}.junction_sites.tsv.gz"
    counts = write_sites(o(sites_name), lifted, a.version, list(order))
    shutil.copy(os.path.join(HERE, "README.md"), o("README.md"))
    meta = {"name": "cxsv_literature", "version": a.version,
            "created_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "code_commit": a.code_commit,
            "inputs": [dict(source_id=r["source_id"], file=r["file"], sha256=r["sha256"], role=r["role"])
                       for r in L.read_tsv(os.path.join(HERE, "sources.tsv")) if r["sha256"] != "."],
            "curated_tables": {f: sha256(os.path.join(a.curated, f)) for f in sorted(os.listdir(a.curated))
                               if f.endswith(".tsv")},
            "chains": {"hg19ToHg38": sha256(a.chain_hg19), "hg38ToHs1": sha256(a.chain_hg38)},
            "counts": {"events": len(events), "junctions": len(lifted),
                       "junctions_lifted": sum(j["lift_status"] == "lifted" for j in lifted),
                       "junction_sites": dict(counts), "clinvar_nstd102": cv_counts}}
    with open(o("literature_set.json"), "w") as fh:
        json.dump(meta, fh, indent=2, sort_keys=True)
    files = sorted(f for f in os.listdir(a.out) if f != "MANIFEST.tsv" and os.path.isfile(o(f)))
    with open(o("MANIFEST.tsv"), "w") as fh:
        fh.write("file\tbytes\tsha256\n")
        for f in files:
            fh.write(f"{f}\t{os.path.getsize(o(f))}\t{sha256(o(f))}\n")
    print(json.dumps(meta["counts"], indent=1), file=sys.stderr)


if __name__ == "__main__":
    main()

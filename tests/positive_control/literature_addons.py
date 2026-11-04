#!/usr/bin/env python3
"""Design and build the L-series positive-control add-ons: published structures re-planted (item 11, test B).

Four literature events (docs/LITERATURE_SET.md, "Detection test") are planted
into the positive-control host 3559 with the existing simulator, as add-ons
spliced in at run time exactly like E13 (tests/positive_control/make_spikein_reads.py,
tests/lab/pos_control_full_genome.sh ADDONS=...). The release CRAM and the
E01-E13 truth are not touched.

  L01  BAH19_BAB8964  PLP1 DUP-TRP/INV-DUP (Pelizaeus-Merzbacher; Bahrambeigi 2019)  full mode
  L02  NP18_P2109_123 17p13.3 delINVdel (Nazaryan-Petersen 2018)                      full mode
  L03  COL17_UTR22    balanced germline chromothripsis, chr2/7/8/9 (Collins 2017)     windowed mode
  L04  EIS24_P6       Xq28 deletion with an inverted 9q22.33 insertion (ALD; Eisfeldt 2024) full mode

Full mode writes an events.yaml-format catalogue (events_literature.yaml) whose
segments come from the literature set's hs1 junctions, and builds it with the
unchanged build_haplotypes.py: the whole event span is simulated (derivative +
reference window) and the host reads over span +/- flank_replace are replaced.

Windowed mode (balanced events with Mb-sized fragments, L03): the derivative is
cut inside every fragment longer than 2 x --window into independent contigs, so
only sequence within --window bp of a breakend is simulated: alt.fa = those
derivative contigs, ref.fa = every breakend +/- window (and every deleted
segment), replace.bed = every breakend +/- --replace (and every deleted segment).
In a balanced rearrangement the sequence further than a read length from every
breakend is unchanged, so this gives the same reads near the junctions as
simulating the whole 42 Mb derivative.

Placement check: every L simulation window must lie outside every E01-E13
simulation window (events.yaml segments incl. flanks) and outside the other L
events'. The script stops otherwise.

Usage:
  literature_addons.py --literature LITDIR --fasta hs1.fa --events tests/positive_control/events.yaml \\
      --out BUILD [--window 100000] [--replace 40000]
then per event: build_haplotypes.py (full mode; done here) and make_spikein_reads.py --build BUILD/<L>.
"""
import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]

FLANK_SIM, FLANK_REPLACE = 150_000, 50_000
ADDONS = {
    "L01": dict(literature="BAH19_BAB8964", name="DUP-TRP/INV-DUP (PLP1, Pelizaeus-Merzbacher)", gene="PLP1",
                mode="dup_trp_inv_dup", flank=FLANK_SIM),
    "L02": dict(literature="NP18_P2109_123", name="delINVdel (17p13.3)", gene=None, mode="del_inv_del",
                flank=FLANK_SIM),
    "L03": dict(literature="COL17_UTR22", name="balanced germline chromothripsis (chr2/7/8/9)", gene=None,
                mode="windowed"),
    "L04": dict(literature="EIS24_P6", name="deletion with inverted interchromosomal insertion (Xq28; 9q22.33, ALD)",
                gene="BCAP31", mode="del_ins_inv", flank=100_000),
}


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader([ln for ln in fh if not ln.startswith("#")], delimiter="\t"))


def junctions_of(lit, event):
    js = [j for j in read_tsv(f"{lit}/literature_junctions.tsv") if j["event_id"] == event]
    if not js or any(j["lift_status"] != "lifted" for j in js):
        sys.exit(f"{event}: junctions missing or not all lifted")
    for j in js:
        j["pos1"], j["pos2"] = int(j["pos1"]), int(j["pos2"])
    return sorted(js, key=lambda j: j["junction_id"])


# ─── full-mode templates (segments 0-based half-open, as events.yaml) ────────
def tmpl_dup_trp_inv_dup(js, flank):
    """DUP1 B, TRP C, DUP2 D: der p B C D C- B C D q (E01 structure).
    J (tail-to-tail) = B start / C start; J (head-to-head) = C end / D end."""
    t2t = [j for j in js if j["side1"] == j["side2"] == "-"][0]
    h2h = [j for j in js if j["side1"] == j["side2"] == "+"][0]
    c = t2t["chrom1"]
    a, b = t2t["pos1"] - 1, t2t["pos2"] - 1          # B start (0-based), C start
    cc, d = h2h["pos1"], h2h["pos2"]                 # C end, D end (1-based last base = 0-based end)
    segs = {"p": [c, a - flank, a], "B": [c, a, b], "C": [c, b, cc], "D": [c, cc, d], "q": [c, d, d + flank]}
    return segs, {"der": ["p+", "B+", "C+", "D+", "C-", "B+", "C+", "D+", "q+"]}


def tmpl_del_inv_del(js, flank):
    """delINVdel: J h2h (x1 +, y1 +), J t2t (x2 -, y2 -): der p I- q; D1 = (x1, x2), D2 = (y1, y2) deleted."""
    h2h = [j for j in js if j["side1"] == j["side2"] == "+"][0]
    t2t = [j for j in js if j["side1"] == j["side2"] == "-"][0]
    c = h2h["chrom1"]
    x1, y1, x2, y2 = h2h["pos1"], h2h["pos2"], t2t["pos1"], t2t["pos2"]
    segs = {"p": [c, x1 - flank, x1], "D1": [c, x1, x2 - 1], "I": [c, x2 - 1, y1], "D2": [c, y1, y2 - 1],
            "q": [c, y2 - 1, y2 - 1 + flank]}
    return segs, {"der": ["p+", "I-", "q+"]}


def tmpl_del_ins_inv(js, flank):
    """Deletion D on the recipient with an inverted copy of donor piece I: der(R) p I- q; the donor
    chromosome of the same haplotype stays intact (path 'donor': p_d I q_d)."""
    hh = [j for j in js if j["side1"] == j["side2"] == "+"][0]
    tt = [j for j in js if j["side1"] == j["side2"] == "-"][0]
    rc = "chrX" if "chrX" in (hh["chrom1"], hh["chrom2"]) else hh["chrom1"]

    def ends(j):
        return ((j["chrom1"], j["pos1"]), (j["chrom2"], j["pos2"])) if j["chrom1"] == rc else \
               ((j["chrom2"], j["pos2"]), (j["chrom1"], j["pos1"]))
    (_, r1), (dc, d_end) = ends(hh)                  # recipient end retained to r1; donor piece end
    (_, r2), (_, d_start) = ends(tt)                 # donor piece start; recipient resumes at r2
    segs = {f"p_{rc}": [rc, r1 - flank, r1], "D": [rc, r1, r2 - 1], f"q_{rc}": [rc, r2 - 1, r2 - 1 + flank],
            f"p_{dc}": [dc, d_start - 1 - flank, d_start - 1], "I": [dc, d_start - 1, d_end],
            f"q_{dc}": [dc, d_end, d_end + flank]}
    return segs, {f"der({rc[3:]})": [f"p_{rc}+", "I-", f"q_{rc}+"],
                  f"normal({dc[3:]})": [f"p_{dc}+", "I+", f"q_{dc}+"]}


TEMPLATES = {"dup_trp_inv_dup": tmpl_dup_trp_inv_dup, "del_inv_del": tmpl_del_inv_del,
             "del_ins_inv": tmpl_del_ins_inv}


# ─── windowed mode ───────────────────────────────────────────────────────────
def derivative_contigs(frags, window):
    """Cut a derivative (hs1 fragments, 5'->3') into contigs holding all sequence within `window`
    of a junction. Returns [[(chrom, start0, end0, strand), ...], ...] (0-based half-open pieces)."""
    contigs, cur = [], []
    n = len(frags)
    for i, f in enumerate(frags):
        c, s, e, x = f["chrom"], f["start0"], f["end0"], f["strand"]
        long_ = (e - s) > 2 * window
        if not long_:
            cur.append((c, s, e, x))
            continue
        head = (c, s, s + window, "+") if x == "+" else (c, e - window, e, "-")
        tail = (c, e - window, e, "+") if x == "+" else (c, s, s + window, "-")
        if i > 0:
            cur.append(head)
            contigs.append(cur)
        cur = [tail] if i < n - 1 else []
    if len(cur) > 1:
        contigs.append(cur)
    return [k for k in contigs if len(k) > 1]


def piece_seq(fa, p):
    c, s, e, x = p
    q = fa.fetch(c, s, e).upper()
    return q if x == "+" else revcomp(q)


def merge(iv):
    out = []
    for c, s, e in sorted(iv):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return out


def build_windowed(lid, spec, lit, fa, out, window, replace):
    ev = spec["literature"]
    rows = [r for r in read_tsv(f"{lit}/literature_derivatives.tsv") if r["event_id"] == ev]
    if any(r["lift_status"] != "lifted" for r in rows):
        sys.exit(f"{ev}: derivative fragments not all lifted")
    ders = defaultdict(list)
    for r in sorted(rows, key=lambda r: (r["derivative"], int(r["order"]))):
        c = r["chrom"]
        ln = fa.get_reference_length(c)
        s = 1 if r["start"] == "pter" else int(r["start"])
        e = ln if r["end"] == "qter" else int(r["end"])
        ders[r["derivative"]].append(dict(chrom=c, start0=s - 1, end0=e, strand=r["strand_hs1"], label=r["fragment"]))
    d = f"{out}/{lid}"
    os.makedirs(d, exist_ok=True)
    truth, bk = [], []
    with open(f"{d}/alt.fa", "w") as fh:
        for der, frags in ders.items():
            for x, y in zip(frags, frags[1:]):
                lpos, lside = (x["end0"], "end") if x["strand"] == "+" else (x["start0"] + 1, "start")
                rpos, rside = (y["start0"] + 1, "start") if y["strand"] == "+" else (y["end0"], "end")
                truth.append(dict(event=lid, derivative=der, left=f"{x['label']}{x['strand']}",
                                  right=f"{y['label']}{y['strand']}", chrom1=x["chrom"], pos1=lpos, side1=lside,
                                  chrom2=y["chrom"], pos2=rpos, side2=rside))
                bk += [(x["chrom"], lpos), (y["chrom"], rpos)]
            for k, contig in enumerate(derivative_contigs(frags, window), 1):
                seq = "".join(piece_seq(fa, p) for p in contig)
                fh.write(f">{lid}_{der.replace('(', '_').replace(')', '')}_c{k} "
                         + " ".join(f"{c}:{s}-{e}{x}" for c, s, e, x in contig) + "\n")
                for i in range(0, len(seq), 80):
                    fh.write(seq[i:i + 80] + "\n")
    dels = [(r["chrom"], int(r["start"]) - 1, int(r["end"])) for r in read_tsv(f"{lit}/literature_cn_segments.tsv")
            if r["event_id"] == ev and r["lift_status"] == "lifted" and float(r["delta"]) < 0]
    refw = merge([(c, max(0, p - window), p + window) for c, p in bk] + [(c, max(0, s - window), e + window)
                                                                         for c, s, e in dels])
    repw = merge([(c, max(0, p - replace), p + replace) for c, p in bk] + [(c, s, e) for c, s, e in dels])
    with open(f"{d}/ref.fa", "w") as fh:
        for c, s, e in refw:
            seq = fa.fetch(c, s, e).upper()
            fh.write(f">{lid}_ref_{c}_{s} {c}:{s}-{e}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    with open(f"{d}/replace.bed", "w") as fh:
        for c, s, e in repw:
            fh.write(f"{c}\t{s}\t{e}\t{lid}\n")
    spans = defaultdict(list)
    for c, p in bk:
        spans[c].append(p)
    ev_rows = [f"{c}\t{min(p)}\t{max(p)}\t{lid}\t{spec['name'].replace(' ', '_')}\n" for c, p in sorted(spans.items())]
    return truth, ev_rows, [(c, s, e) for c, s, e in refw], dels


def write_truth(tdir, truth, ev_rows, fa, cn_rows=()):
    from build_haplotypes import bnd_alt     # needs pysam; imported only when building
    os.makedirs(tdir, exist_ok=True)
    cols = ["event", "derivative", "left", "right", "chrom1", "pos1", "side1", "chrom2", "pos2", "side2"]
    with open(f"{tdir}/truth_junctions.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for j in truth:
            fh.write("\t".join(str(j[c]) for c in cols) + "\n")
    with open(f"{tdir}/truth_events.bed", "w") as fh:
        fh.writelines(sorted(ev_rows))
    with open(f"{tdir}/copy_number.tsv", "w") as fh:
        fh.write("event\tsegment\tchrom\tstart\tend\tcopies_on_derivative\n")
        fh.writelines(cn_rows)
    with open(f"{tdir}/truth.vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=CxSV literature add-ons (tests/positive_control/literature_addons.py)\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpos_control_3559\n")
        for i, j in enumerate(sorted(truth, key=lambda j: (j["chrom1"], int(j["pos1"]))), 1):
            rb = fa.fetch(j["chrom1"], int(j["pos1"]) - 1, int(j["pos1"])).upper()
            fh.write(f"{j['chrom1']}\t{j['pos1']}\t{j['event']}_J{i}\t{rb}\t{bnd_alt(j, rb)}\t.\tPASS\t"
                     f"SVTYPE=BND;EVENT={j['event']}\tGT\t0/1\n")


def sim_windows(ev):
    by = {}
    for n, (c, s, e) in ev["segments"].items():
        lo, hi = by.get(c, (s, e))
        by[c] = (min(lo, s), max(hi, e))
    return [(c, s, e) for c, (s, e) in by.items()]


def overlaps(a, b):
    return [(x, y) for x in a for y in b if x[0] == y[0] and x[1] < y[2] and y[1] < x[2]]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--literature", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--events", required=True, help="tests/positive_control/events.yaml (E01-E13 placement)")
    ap.add_argument("--e13", default=None, help="add-on events_E13.yaml (placement check)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=100_000)
    ap.add_argument("--replace", type=int, default=40_000)
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    import pysam
    fa = pysam.FastaFile(a.fasta)
    os.makedirs(a.out, exist_ok=True)
    taken = [w for e in yaml.safe_load(open(a.events))["events"] for w in sim_windows(e)]
    if a.e13:
        taken += [w for e in yaml.safe_load(open(a.e13))["events"] for w in sim_windows(e)]
    lit_events = {r["event_id"]: r for r in read_tsv(f"{a.literature}/literature_events.tsv")}
    full, placed = [], {}
    for lid, spec in ADDONS.items():
        if a.only and lid not in a.only.split(","):
            continue
        le = lit_events[spec["literature"]]
        meta = dict(literature_event=spec["literature"], source=le["citation"], doi=le["doi"],
                    case=le["case_id"], phenotype=le["phenotype"], author_classification=le["author_classification"])
        if spec["mode"] == "windowed":
            truth, ev_rows, wins, dels = build_windowed(lid, spec, a.literature, fa, a.out, a.window, a.replace)
            cn = [f"{lid}\tdel{i}\t{c}\t{s}\t{e}\t0\n" for i, (c, s, e) in enumerate(dels, 1)]
            write_truth(f"{a.out}/{lid}/truth", truth, ev_rows, fa, cn)
            with open(f"{a.out}/{lid}/addon.yaml", "w") as fh:
                yaml.safe_dump(dict(id=lid, name=spec["name"], mode="windowed", window=a.window, replace=a.replace,
                                    junctions=len(truth), literature=meta,
                                    expect=dict(cxsv=True, panel_status="private", score="top")), fh, sort_keys=False)
            placed[lid] = wins
            print(f"{lid}: windowed, {len(truth)} junctions, {len(wins)} reference windows "
                  f"({sum(e - s for _, s, e in wins):,} bp)", file=sys.stderr)
            continue
        segs, paths = TEMPLATES[spec["mode"]](junctions_of(a.literature, spec["literature"]), spec["flank"])
        ev = dict(id=lid, name=spec["name"], gene=spec["gene"], literature=meta,
                  note=f"re-planted from the literature set ({spec['literature']}); heterozygous in the female host",
                  expect=dict(cxsv=True, panel_status="private", score="top"), zygosity="heterozygous",
                  segments=segs, paths=paths)
        full.append(ev)
        placed[lid] = sim_windows(ev)
    for lid, wins in placed.items():
        clash = overlaps(wins, taken) + [x for o, w in placed.items() if o != lid for x in overlaps(wins, w)]
        if clash:
            sys.exit(f"{lid}: simulation window overlaps an existing event: {clash[:3]}")
    if full:
        doc = dict(version=1, reference="T2T-CHM13v2 (UCSC hs1 names); 0-based half-open coordinates",
                   host="sample 3559 (exempt), realigned with minimap2 2.24-r1122 k17", seed=3559,
                   flank_sim=FLANK_SIM, flank_replace=FLANK_REPLACE,
                   notation="as tests/positive_control/events.yaml; segments from the literature set (hs1)",
                   events=full)
        path = f"{a.out}/events_literature.yaml"
        with open(path, "w") as fh:
            yaml.safe_dump(doc, fh, sort_keys=False)
        subprocess.run([sys.executable, f"{HERE}/build_haplotypes.py", path, a.fasta, f"{a.out}/full", "--check"],
                       check=True)
        for ev in full:
            lid = ev["id"]
            os.makedirs(f"{a.out}/{lid}/truth", exist_ok=True)
            for f in ("alt.fa", "ref.fa", "replace.bed"):
                os.replace(f"{a.out}/full/{lid}/{f}", f"{a.out}/{lid}/{f}")
            for f in ("truth_junctions.tsv", "truth_events.bed", "copy_number.tsv"):
                lines = open(f"{a.out}/full/truth/{f}").read().splitlines(True)
                keep = [ln for i, ln in enumerate(lines) if (i == 0 and f.endswith(".tsv")) or lid in ln.split("\t")]
                open(f"{a.out}/{lid}/truth/{f}", "w").writelines(keep)
            with open(f"{a.out}/{lid}/addon.yaml", "w") as fh:
                yaml.safe_dump(dict(ev, mode="full"), fh, sort_keys=False)
    print("placement check: no L simulation window overlaps E01-E13 or another L event", file=sys.stderr)


if __name__ == "__main__":
    main()

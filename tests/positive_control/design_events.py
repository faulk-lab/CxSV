#!/usr/bin/env python3
"""Place the positive-control CxSV catalogue on hs1 and write events.yaml.

Each event is a TEMPLATE: an ordered list of contiguous reference segments with
nominal lengths, and one or more derivative paths (segment order + strand) that
describe the rearranged haplotype. This script anchors every template at its
locus and moves each breakpoint (segment boundary) to the nearest position, within
+/- MAX_SHIFT bp, that lies in unambiguous sequence:

  * no segmental duplication within SEGDUP_PAD bp, and
  * at most MAX_REPEAT_FRAC of the +/- REPEAT_PAD bp window soft-masked
    (lower case) in the hs1 FASTA, i.e. RepeatMasker/TRF-masked.

The soft mask is used rather than resources/repeats.final.bed: that BED holds
UCSC "joined" RepeatMasker elements whose start/end span whole fragmented
repeats including the gaps between fragments, so it reports e.g. most of MECP2
as LINE although only ~2% of it is masked.

Gene events additionally must disrupt at least one exon of their gene (an exon
overlaps a deleted, inverted or copy-changed segment); the anchor is slid until
they do. The result, events.yaml, is the single source of truth: haplotypes,
truth set and figure are all generated from it, never from this script.

Usage (on a machine with the annotation BEDs):
  design_events.py --resources /mnt/sargasso1/CxSV/resources --fasta hs1.fa \
                   --out tests/positive_control/events.yaml [--seed 3559]
"""
import argparse
import bisect
import random
import sys
from collections import defaultdict

import pysam
import yaml

MAX_SHIFT, SEGDUP_PAD, REPEAT_PAD, MAX_REPEAT_FRAC = 5000, 10_000, 250, 0.30
FLANK_SIM, FLANK_REPLACE = 150_000, 50_000

# ─── catalogue: (id, name, locus, segments [(name, length, role)], paths) ─────
# role: keep | del | inv | dup | src (template source) | lost — documentation only;
# the derivative paths define the structure. "p"/"q" are the flanks.
def catalogue(rng):
    shatter = [(f"S{i}", rng.randint(26_000, 40_000), "frag") for i in range(1, 13)]
    return [
        dict(id="E01", name="DUP-TRP/INV-DUP", gene="MECP2", chrom="chrX",
             note="duplicated B and D flank a triplicated C whose middle copy is inverted "
                  "(replication-based, inverted-repeat mediated; MECP2/PLP1 type)",
             segments=[("B", 6000, "dup"), ("C", 12000, "trp"), ("D", 8000, "dup")],
             paths={"der": ["p+", "B+", "C+", "D+", "C-", "B+", "C+", "D+", "q+"]},
             expect=dict(cxsv=True, ideal_class="DUP+INV")),
        dict(id="E02", name="delINVdel", gene="NF1", chrom="chr17",
             note="inversion with a deletion at each breakpoint",
             segments=[("D1", 4000, "del"), ("I", 18000, "inv"), ("D2", 6000, "del")],
             paths={"der": ["p+", "I-", "q+"]},
             expect=dict(cxsv=True, ideal_class="DEL+INV")),
        dict(id="E03", name="dupINVdup", gene="PMP22", chrom="chr17",
             note="inversion of B with inverted duplications of the flanking segments A and C",
             segments=[("A", 5000, "dup"), ("B", 15000, "inv"), ("C", 6000, "dup")],
             paths={"der": ["p+", "A+", "C-", "B-", "A-", "C+", "q+"]},
             expect=dict(cxsv=True, ideal_class="DUP+INV")),
        dict(id="E04", name="DEL+INV+DUP", gene="SHANK3", chrom="chr22",
             note="deletion of B, inversion of C and tandem duplication of E within 24 kb",
             segments=[("A", 3000, "keep"), ("B", 4000, "del"), ("C", 9000, "inv"),
                       ("D", 3000, "keep"), ("E", 5000, "dup")],
             paths={"der": ["p+", "A+", "C-", "D+", "E+", "E+", "q+"]},
             expect=dict(cxsv=True, ideal_class="DEL+DUP+INV")),
        dict(id="E05", name="nested inversion with inner deletion", gene="NRXN1", chrom="chr2",
             note="35 kb inversion containing a second, inner inversion (C, back in "
                  "forward orientation) and a deletion (D)",
             segments=[("A", 8000, "inv"), ("B", 7000, "inv"), ("C", 8000, "inv2"),
                       ("D", 3000, "del"), ("E", 9000, "inv")],
             paths={"der": ["p+", "E-", "C+", "B-", "A-", "q+"]},
             expect=dict(cxsv=True, ideal_class="DEL+INV")),
        dict(id="E06", name="templated insertions at a deletion (FoSTeS/MMBIR)", gene="AUTS2", chrom="chr7",
             note="20 kb deletion X whose junction carries a forward copy of T1 (from ~14 kb "
                  "upstream) and an inverted copy of T2 (from ~30 kb downstream)",
             segments=[("T1", 1200, "src"), ("U1", 13800, "keep"), ("X", 20000, "del"),
                       ("U2", 30000, "keep"), ("T2", 900, "src")],
             paths={"der": ["p+", "T1+", "U1+", "T1+", "T2-", "U2+", "T2+", "q+"]},
             expect=dict(cxsv=True, ideal_class="INS+DEL")),
        dict(id="E07", name="reciprocal translocation with junction deletions",
             gene="SATB2/FOXP2", chrom=["chr2", "chr7"],
             note="t(2;7) between SATB2 and FOXP2 introns; 2 kb lost on chr2 and 1 kb on chr7",
             segments={"chr2": [("d2", 2000, "del")], "chr7": [("d7", 1000, "del")]},
             paths={"der(2)": ["p_chr2+", "q_chr7+"], "der(7)": ["p_chr7+", "q_chr2+"]},
             expect=dict(cxsv=True, ideal_class="BND_cluster")),
        dict(id="E08", name="chromoplexy (balanced 3-chromosome chain)",
             gene="NSD2/EHMT1/ZIC2", chrom=["chr4", "chr9", "chr13"],
             note="chained translocations chr4->chr9->chr13->chr4 with 500 bp losses at each junction",
             segments={"chr4": [("d4", 500, "del")], "chr9": [("d9", 500, "del")],
                       "chr13": [("d13", 500, "del")]},
             paths={"der(4)": ["p_chr4+", "q_chr9+"], "der(9)": ["p_chr9+", "q_chr13+"],
                    "der(13)": ["p_chr13+", "q_chr4+"]},
             expect=dict(cxsv=True, ideal_class="BND_cluster",
                         caveat="each chromosome carries only two breakends, so per-locus "
                                "criteria C1-C5 may not fire")),
        dict(id="E09", name="localized chromothripsis", gene="CHD2", chrom="chr15",
             note="~400 kb shattered into 12 fragments; S4, S8 and S12 lost; the rest "
                  "re-joined in random order and orientation",
             segments=shatter,
             paths={"der": ["p+", "S1+", "S5-", "S3+", "S9-", "S7+", "S2+", "S11-", "S6+", "S10-", "q+"]},
             expect=dict(cxsv=True, ideal_class="Chromothripsis_like")),
        dict(id="E10", name="multi-breakpoint cluster", gene=None, chrom="chr5",
             note="two inversions (C, E) and two deletions (B, F) within 39 kb, intergenic",
             segments=[("A", 6000, "keep"), ("B", 3000, "del"), ("C", 7000, "inv"),
                       ("D", 5000, "keep"), ("E", 8000, "inv"), ("F", 4000, "del"), ("G", 6000, "keep")],
             paths={"der": ["p+", "A+", "C-", "D+", "E-", "G+", "q+"]},
             expect=dict(cxsv=True, ideal_class="Multi_breakpoint_cluster")),
        dict(id="E11", name="simple deletion (specificity control)", gene=None, chrom="chr10",
             note="a single 20 kb deletion: must NOT be classified a CxSV",
             segments=[("X", 20000, "del")], paths={"der": ["p+", "q+"]},
             expect=dict(cxsv=False, ideal_class=None)),
        dict(id="E12", name="simple inversion (specificity control)", gene=None, chrom="chr12",
             note="a single 60 kb balanced inversion: must NOT be classified a CxSV",
             segments=[("I", 60000, "inv")], paths={"der": ["p+", "I-", "q+"]},
             expect=dict(cxsv=False, ideal_class=None)),
    ]


# ─── annotation ──────────────────────────────────────────────────────────────
def read_bed(path, chroms=None, names=False):
    out = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if chroms and f[0] not in chroms:
                continue
            out[f[0]].append((int(f[1]), int(f[2]), f[3] if names and len(f) > 3 else ""))
    for c in out:
        out[c].sort()
    return out


class Track:
    """Sorted intervals with fast 'any overlap' and 'covered bp' queries."""
    def __init__(self, ivs):
        merged = []
        for s, e, *_ in sorted(ivs):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        self.s = [m[0] for m in merged]
        self.e = [m[1] for m in merged]

    def covered(self, a, b):
        i = max(0, bisect.bisect_right(self.s, a) - 1)
        tot = 0
        while i < len(self.s) and self.s[i] < b:
            tot += max(0, min(b, self.e[i]) - max(a, self.s[i]))
            i += 1
        return tot


class SoftMask:
    """Fraction of soft-masked (lower-case) bases in a window of the FASTA."""
    def __init__(self, fasta):
        self.fa = pysam.FastaFile(fasta)

    def covered(self, chrom, a, b):
        seq = self.fa.fetch(chrom, max(0, a), b)
        return sum(ch.islower() or ch in "Nn" for ch in seq)


def good_breakpoint(chrom, x, segdup, rep):
    if segdup[chrom].covered(x - SEGDUP_PAD, x + SEGDUP_PAD) > 0:
        return False
    return rep.covered(chrom, x - REPEAT_PAD, x + REPEAT_PAD) <= MAX_REPEAT_FRAC * 2 * REPEAT_PAD


def adjust(chrom, x, segdup, rep):
    for d in range(0, MAX_SHIFT + 1, 25):
        for y in (x - d, x + d):
            if good_breakpoint(chrom, y, segdup, rep):
                return y
    return None


def place_linear(chrom, anchor, segs, segdup, rep):
    """Contiguous segments from anchor; every boundary moved to good sequence."""
    bounds = [anchor]
    for _, length, _ in segs:
        bounds.append(bounds[-1] + length)
    placed = []
    for b in bounds:
        y = adjust(chrom, b, segdup, rep)
        if y is None or (placed and y - placed[-1] < 0.6 * (b - bounds[len(placed) - 1])):
            return None
        placed.append(y)
    return {name: [chrom, placed[i], placed[i + 1]] for i, (name, _, _) in enumerate(segs)}


def disrupts_exon(seg_coords, segs, exons, gene):
    for (name, _, role) in segs:
        if role in ("keep", "src"):
            continue
        c, s, e = seg_coords[name]
        for es, ee, g in exons.get(c, []):
            if ee > s and es < e and gene in g.split(","):
                return True
    return False


def gene_span(genes, gene):
    for c, ivs in genes.items():
        hits = [(s, e) for s, e, g in ivs if gene in g.split(",")]
        if hits:
            return c, min(h[0] for h in hits), max(h[1] for h in hits)
    raise SystemExit(f"gene {gene} not found")


def gene_free(genes, chrom, s, e, pad=200_000):
    return not any(ge > s - pad and gs < e + pad for gs, ge, _ in genes.get(chrom, []))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--resources", required=True)
    ap.add_argument("--fasta", required=True, help="soft-masked hs1 FASTA (+ .fai)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=3559)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    cat = catalogue(rng)
    chroms = set()
    for ev in cat:
        chroms.update(ev["chrom"] if isinstance(ev["chrom"], list) else [ev["chrom"]])
    clen = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fasta + ".fai")}
    genes = read_bed(f"{a.resources}/genes.final.bed", chroms, names=True)
    exons = read_bed(f"{a.resources}/exons.filtered.bed", chroms, names=True)
    segdup = {c: Track(v) for c, v in read_bed(f"{a.resources}/segdup.final.bed", chroms).items()}
    rep = SoftMask(a.fasta)
    for c in chroms:
        segdup.setdefault(c, Track([]))

    events, used = [], defaultdict(list)

    def clash(chrom, s, e):
        return any(e + 2 * FLANK_SIM > us and s - 2 * FLANK_SIM < ue for us, ue in used[chrom])

    for ev in cat:
        out = {k: ev[k] for k in ("id", "name", "gene", "note", "expect")}
        out["zygosity"] = "heterozygous"
        if isinstance(ev["chrom"], list):                # inter-chromosomal: one junction locus per chromosome
            segs_out, genes_list = {}, ev["gene"].split("/")
            for chrom, g in zip(ev["chrom"], genes_list):
                c, gs, ge = gene_span(genes, g)
                assert c == chrom, (g, c, chrom)
                placed = None
                for k in range(0, 400):
                    anchor = (gs + ge) // 2 + (k // 2) * 1500 * (1 if k % 2 else -1)
                    placed = place_linear(chrom, anchor, ev["segments"][chrom], segdup, rep)
                    if placed and not clash(chrom, *placed[ev["segments"][chrom][0][0]][1:]):
                        break
                    placed = None
                if not placed:
                    sys.exit(f"{ev['id']}: no placement on {chrom}")
                segs_out.update(placed)
                d = placed[ev["segments"][chrom][0][0]]
                segs_out[f"p_{chrom}"] = [chrom, max(0, d[1] - FLANK_SIM), d[1]]
                segs_out[f"q_{chrom}"] = [chrom, d[2], min(clen[chrom], d[2] + FLANK_SIM)]
                used[chrom].append((d[1], d[2]))
            out["segments"] = segs_out
        else:
            chrom, segs = ev["chrom"], ev["segments"]
            total = sum(l for _, l, _ in segs)
            if ev["gene"]:
                c, gs, ge = gene_span(genes, ev["gene"])
                assert c == chrom
                centre, lo, hi = (gs + ge) // 2, gs - total // 2, ge - total // 2
            else:
                centre, lo, hi = {"chr5": 25_000_000, "chr10": 30_000_000, "chr12": 75_000_000}[chrom], 0, clen[chrom]
            placed = None
            for k in range(0, 4000):
                anchor = centre - total // 2 + (k // 2) * 2000 * (1 if k % 2 else -1)
                if anchor < lo or anchor > hi:
                    continue
                cand = place_linear(chrom, anchor, segs, segdup, rep)
                if not cand:
                    continue
                s, e = cand[segs[0][0]][1], cand[segs[-1][0]][2]
                if clash(chrom, s, e):
                    continue
                if ev["gene"] and not disrupts_exon(cand, segs, exons, ev["gene"]):
                    continue
                if not ev["gene"] and not gene_free(genes, chrom, s, e):
                    continue
                placed = cand
                break
            if not placed:
                sys.exit(f"{ev['id']}: no placement found")
            s, e = placed[segs[0][0]][1], placed[segs[-1][0]][2]
            placed["p"] = [chrom, max(0, s - FLANK_SIM), s]                 # flanks clamped to the
            placed["q"] = [chrom, e, min(clen[chrom], e + FLANK_SIM)]       # chromosome (e.g. SHANK3)
            used[chrom].append((s, e))
            out["segments"] = placed
        out["paths"] = ev["paths"]
        events.append(out)
        print(f"{ev['id']} {ev['name']}: " + ", ".join(
            f"{k}={v[0]}:{v[1]}-{v[2]}" for k, v in out["segments"].items() if not k.startswith(("p", "q"))),
            file=sys.stderr)

    doc = {
        "version": 1,
        "reference": "T2T-CHM13v2 (UCSC hs1 names); 0-based half-open coordinates",
        "host": "sample 3559 (exempt), realigned with minimap2 2.24-r1122 k17",
        "seed": a.seed,
        "design": dict(max_shift=MAX_SHIFT, segdup_pad=SEGDUP_PAD, repeat_pad=REPEAT_PAD,
                       max_repeat_frac=MAX_REPEAT_FRAC, repeat_source="hs1 soft mask"),
        "flank_sim": FLANK_SIM, "flank_replace": FLANK_REPLACE,
        "notation": "segments are reference intervals; a derivative path lists segments 5'->3' "
                    "on the rearranged haplotype with strand (+ reference orientation, - reverse "
                    "complement). p/q are flanks. Segments absent from every path are deleted.",
        "events": events,
    }
    with open(a.out, "w") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, width=100)
    print(f"wrote {a.out} ({len(events)} events)", file=sys.stderr)


if __name__ == "__main__":
    main()

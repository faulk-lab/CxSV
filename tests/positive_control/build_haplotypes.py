#!/usr/bin/env python3
"""Build haplotype sequences and the truth set from events.yaml.

For every event:
  <work>/<ID>/alt.fa      derivative sequence(s): each path's segments concatenated
                          5'->3' (reverse-complemented for '-')
  <work>/<ID>/ref.fa      reference window(s): p.start..q.end on each chromosome
  <work>/<ID>/replace.bed host reads overlapping these windows are replaced
                          (event span +/- flank_replace on each chromosome)
Truth set (in <work>/truth/):
  truth_junctions.tsv     every novel adjacency on the derivative haplotype(s)
  truth.vcf               the same junctions as VCF 4.2 breakends (one record per junction)
  truth_events.bed        per event and chromosome: span of rearranged sequence + id + type
  copy_number.tsv         per segment: copies on the derivative haplotype (0 = deleted)

Every derivative is re-checked against the reference after construction
(--check), so the sequences used for simulation provably match events.yaml.

Usage: build_haplotypes.py events.yaml hs1.fa <work> [--check]
"""
import argparse
import os
import sys
from collections import Counter

import pysam
import yaml

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def parse(el):
    return el[:-1], el[-1]


def seg_seq(fa, seg, strand):
    c, s, e = seg
    x = fa.fetch(c, s, e).upper()
    return x if strand == "+" else revcomp(x)


def junctions(ev):
    """Novel adjacencies along each derivative path."""
    segs = ev["segments"]
    out = []
    for der, path in ev["paths"].items():
        for (xa, sa), (ya, sb) in zip(map(parse, path), map(parse, path[1:])):
            X, Y = segs[xa], segs[ya]
            if sa == sb == "+" and X[0] == Y[0] and X[2] == Y[1]:
                continue                         # reference adjacency
            if sa == sb == "-" and X[0] == Y[0] and X[1] == Y[2]:
                continue
            # left breakend
            if sa == "+":
                lpos, lside = X[2], "end"        # last retained base = X.end (1-based)
            else:
                lpos, lside = X[1] + 1, "start"  # first base of X (1-based)
            if sb == "+":
                rpos, rside = Y[1] + 1, "start"
            else:
                rpos, rside = Y[2], "end"
            out.append(dict(event=ev["id"], derivative=der, left=f"{xa}{sa}", right=f"{ya}{sb}",
                            chrom1=X[0], pos1=lpos, side1=lside, chrom2=Y[0], pos2=rpos, side2=rside))
    return out


def bnd_alt(j, ref_base):
    p = f"{j['chrom2']}:{j['pos2']}"
    if j["side1"] == "end":                      # retained sequence is to the left of pos1
        return f"{ref_base}[{p}[" if j["side2"] == "start" else f"{ref_base}]{p}]"
    return f"[{p}[{ref_base}" if j["side2"] == "start" else f"]{p}]{ref_base}"


def chrom_windows(ev):
    """Per chromosome: (sim_start, sim_end, span_start, span_end)."""
    by = {}
    for name, (c, s, e) in ev["segments"].items():
        flank = name == "p" or name == "q" or name.startswith(("p_", "q_"))
        lo, hi, ss, se = by.get(c, (s, e, None, None))
        lo, hi = min(lo, s), max(hi, e)
        if not flank:
            ss = s if ss is None else min(ss, s)
            se = e if se is None else max(se, e)
        by[c] = (lo, hi, ss, se)
    for c, (lo, hi, ss, se) in by.items():
        if ss is None:                           # inter-chromosomal: span = the deleted junction piece
            ss, se = hi, lo
            for name, (cc, s, e) in ev["segments"].items():
                if cc == c and not name.startswith(("p_", "q_")):
                    ss, se = min(ss, s), max(se, e)
        by[c] = (lo, hi, ss, se)
    return by


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("events")
    ap.add_argument("fasta")
    ap.add_argument("work")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    doc = yaml.safe_load(open(a.events))
    fa = pysam.FastaFile(a.fasta)
    fr = doc["flank_replace"]
    os.makedirs(f"{a.work}/truth", exist_ok=True)
    all_j, cn_rows, ev_rows = [], [], []
    for ev in doc["events"]:
        d = f"{a.work}/{ev['id']}"
        os.makedirs(d, exist_ok=True)
        segs = ev["segments"]
        with open(f"{d}/alt.fa", "w") as fh:
            for der, path in ev["paths"].items():
                seq = "".join(seg_seq(fa, segs[n], s) for n, s in map(parse, path))
                if a.check:
                    pos = 0
                    for n, s in map(parse, path):
                        c, x, y = segs[n]
                        piece = seq[pos:pos + (y - x)]
                        ref = fa.fetch(c, x, y).upper()
                        assert piece == (ref if s == "+" else revcomp(ref)), (ev["id"], der, n)
                        pos += y - x
                    assert pos == len(seq)
                tag = der.replace("(", "_").replace(")", "")
                fh.write(f">{ev['id']}_{tag} {' '.join(path)}\n")
                for k in range(0, len(seq), 80):
                    fh.write(seq[k:k + 80] + "\n")
        wins = chrom_windows(ev)
        with open(f"{d}/ref.fa", "w") as fh, open(f"{d}/replace.bed", "w") as fb:
            for c, (lo, hi, ss, se) in sorted(wins.items()):
                seq = fa.fetch(c, lo, hi).upper()
                fh.write(f">{ev['id']}_ref_{c} {c}:{lo}-{hi}\n")
                for k in range(0, len(seq), 80):
                    fh.write(seq[k:k + 80] + "\n")
                fb.write(f"{c}\t{max(0, ss - fr)}\t{se + fr}\t{ev['id']}\n")
                ev_rows.append(f"{c}\t{ss}\t{se}\t{ev['id']}\t{ev['name'].replace(' ', '_')}\n")
        copies = Counter(parse(el)[0] for path in ev["paths"].values() for el in path)
        for n, (c, s, e) in segs.items():
            if n == "p" or n == "q" or n.startswith(("p_", "q_")):
                continue
            cn_rows.append(f"{ev['id']}\t{n}\t{c}\t{s}\t{e}\t{copies.get(n, 0)}\n")
        js = junctions(ev)
        all_j += js
        print(f"{ev['id']}: {len(js)} junctions, derivatives "
              + ", ".join(f"{der}={sum(segs[parse(x)[0]][2] - segs[parse(x)[0]][1] for x in p):,} bp"
                          for der, p in ev["paths"].items()), file=sys.stderr)

    t = f"{a.work}/truth"
    cols = ["event", "derivative", "left", "right", "chrom1", "pos1", "side1", "chrom2", "pos2", "side2"]
    with open(f"{t}/truth_junctions.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for j in all_j:
            fh.write("\t".join(str(j[c]) for c in cols) + "\n")
    with open(f"{t}/truth_events.bed", "w") as fh:
        fh.writelines(sorted(ev_rows))
    with open(f"{t}/copy_number.tsv", "w") as fh:
        fh.write("event\tsegment\tchrom\tstart\tend\tcopies_on_derivative\n")
        fh.writelines(cn_rows)
    contigs = {c: fa.get_reference_length(c) for c in fa.references}
    with open(f"{t}/truth.vcf", "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=CxSV positive control (tests/positive_control/build_haplotypes.py)\n")
        fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        fh.write('##INFO=<ID=EVENT,Number=1,Type=String,Description="Planted CxSV event id">\n')
        fh.write('##INFO=<ID=DERIV,Number=1,Type=String,Description="Derivative haplotype path">\n')
        fh.write('##INFO=<ID=JOIN,Number=1,Type=String,Description="Segments joined (left|right)">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        for c in sorted({j["chrom1"] for j in all_j} | {j["chrom2"] for j in all_j}):
            fh.write(f"##contig=<ID={c},length={contigs[c]}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpos_control_3559\n")
        for i, j in enumerate(sorted(all_j, key=lambda j: (j["chrom1"], j["pos1"])), 1):
            rb = fa.fetch(j["chrom1"], j["pos1"] - 1, j["pos1"]).upper()
            fh.write(f"{j['chrom1']}\t{j['pos1']}\t{j['event']}_J{i}\t{rb}\t{bnd_alt(j, rb)}\t.\tPASS\t"
                     f"SVTYPE=BND;EVENT={j['event']};DERIV={j['derivative'].replace('(', '_').replace(')', '')};"
                     f"JOIN={j['left']}|{j['right']}\tGT\t0/1\n")
    print(f"truth: {len(all_j)} junctions over {len(doc['events'])} events -> {t}", file=sys.stderr)


if __name__ == "__main__":
    main()

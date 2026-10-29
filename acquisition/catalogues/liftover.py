"""Lift junction sites between builds, each breakpoint independently.

Used for GRCh38-only catalogues (gnomAD-SV) with the UCSC chain
hg38ToHs1.over.chain.gz. Every distinct breakpoint (chrom, pos) is written
as a 1-bp BED interval and lifted with UCSC ``liftOver`` (default
-minMatch=0.95, which for a 1-bp feature means the base must sit in an
aligned chain block). Then, per site:

* INS: only pos1 is lifted; pos2 = pos1.
* any end unmapped, or mapped to a contig outside hs1 -> lost (unmapped_end).
* DEL/DUP/INV whose two ends land on different contigs -> lost
  (split_contig). BND ends may land anywhere.
* ends are re-oriented afterwards (a lift into an inverted block can swap
  them); svlen keeps the source value.

``span_changed_gt50pct`` counts kept DEL/DUP/INV sites whose lifted span
differs from the source span by more than 50% (reported, not filtered).
"""
import os
import subprocess
import tempfile
from collections import Counter, defaultdict

from junction_sites import normalise_contig


def lift_points(points, chain, liftover_bin="liftOver", workdir=None):
    """{(chrom, pos)} -> {(chrom, pos): (chrom', pos')} for mapped points."""
    points = sorted(points)
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        bed, out, unm = (os.path.join(tmp, n) for n in ("in.bed", "out.bed", "unmapped.bed"))
        with open(bed, "w") as fh:
            for i, (c, p) in enumerate(points):
                fh.write(f"{c}\t{p - 1}\t{p}\tp{i}\n")
        subprocess.run([liftover_bin, bed, chain, out, unm], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        mapped = {}
        with open(out) as fh:
            for line in fh:
                c, s, e, name = line.rstrip("\n").split("\t")[:4]
                idx = int(name[1:])
                if int(e) - int(s) != 1:
                    continue
                mapped[points[idx]] = (c, int(e))
    return mapped


def apply_lift(sites, mapped):
    """Apply a point map to sites. Returns (kept_sites, stats Counter-dict)."""
    stats = defaultdict(Counter)
    kept = []
    for s in sites:
        t = s["svtype"]
        stats[t]["n_in"] += 1
        m1 = mapped.get((s["chrom1"], s["pos1"]))
        m2 = m1 if t == "INS" else mapped.get((s["chrom2"], s["pos2"]))
        c1 = normalise_contig(m1[0]) if m1 else None
        c2 = normalise_contig(m2[0]) if m2 else None
        if c1 is None or c2 is None:
            stats[t]["lost_unmapped_end"] += 1
            continue
        if t != "BND" and c1 != c2:
            stats[t]["lost_split_contig"] += 1
            continue
        n = dict(s)
        n["chrom1"], n["pos1"] = c1, m1[1]
        n["chrom2"], n["pos2"] = c2, m2[1]
        if t in ("DEL", "DUP", "INV"):
            src = abs(s["pos2"] - s["pos1"])
            new = abs(n["pos2"] - n["pos1"])
            if src > 0 and abs(new - src) > 0.5 * src:
                stats[t]["span_changed_gt50pct"] += 1
        stats[t]["n_out"] += 1
        kept.append(n)
    return kept, stats


def lift_sites(sites, chain, liftover_bin="liftOver", workdir=None):
    pts = set()
    for s in sites:
        pts.add((s["chrom1"], s["pos1"]))
        if s["svtype"] != "INS":
            pts.add((s["chrom2"], s["pos2"]))
    mapped = lift_points(pts, chain, liftover_bin, workdir)
    return apply_lift(sites, mapped)

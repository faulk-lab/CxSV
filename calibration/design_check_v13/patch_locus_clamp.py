#!/usr/bin/env python3
"""Design-check fix (panel-v1.3 item g): clamp a CN locus span to the contig end.

CN segments are built on 1 kb bins, so a segment reaching the last (partial) bin ends past the
contig (HG03667: chr2:242559540-242697000, chr2 length 242,696,752); check_locus then raised and
cn_linkage.py stopped. locus_span gets an optional contig length and clamps the end to it.
  patch_locus_clamp.py REPO"""
import sys

repo = sys.argv[1]
lib = f"{repo}/workflow/lib/cxsv_cnlink.py"
s = open(lib).read()
old_sig = 'def locus_span(c, chrom, rows=()):'
new_sig = 'def locus_span(c, chrom, rows=(), length=None):'
old_doc = '''    (BOC_5883, 17p13.3).\"\"\"'''
new_doc = '''    (BOC_5883, 17p13.3). With `length` (the contig length) the end is clamped to it: segments
    are built on 1 kb bins, so one that reaches the last, partial bin ends past the contig
    (HG03667 chr2, design check panel-v1.3).\"\"\"'''
old_ret = '''    hi = max([s["end"] for s in segs] + [int(r[1]) for r in rows])
    return max(0, lo), hi'''
new_ret = '''    hi = max([s["end"] for s in segs] + [int(r[1]) for r in rows])
    if length is not None:
        hi = min(hi, int(length))
    return max(0, lo), hi'''
for a, b in ((old_sig, new_sig), (old_doc, new_doc), (old_ret, new_ret)):
    assert s.count(a) == 1, a
    s = s.replace(a, b)
open(lib, "w").write(s)

scr = f"{repo}/workflow/scripts/cn_linkage.py"
s = open(scr).read()
a = 's0, e0 = cl.check_locus(chrom, *cl.locus_span(k, chrom, rows_of[r["ClusterID"]]), k, lengths.get(chrom))'
b = ('s0, e0 = cl.check_locus(chrom, *cl.locus_span(k, chrom, rows_of[r["ClusterID"]], lengths.get(chrom)),\n'
     '                                     k, lengths.get(chrom))')
assert s.count(a) == 1
open(scr, "w").write(s.replace(a, b))

t = f"{repo}/tests/unit/test_cn_locus_coords.py"
s = open(t).read()
anchor = "def test_boc5883_loss_gets_segment_anchored_locus(tmp_path):"
assert s.count(anchor) == 1
test = '''def test_locus_span_clamped_to_contig_end():
    # HG03667 (panel-v1.2-r2): a gain reaching the last, partial 1 kb bin of chr2
    gain = _cand([dict(chrom="chr2", start=242_600_000, end=242_697_000, state=3, ploidy=2)])
    rows = [(242_559_540, 242_559_542)]
    assert cl.locus_span(gain, "chr2", rows) == (242_559_540, 242_697_000)
    s, e = cl.locus_span(gain, "chr2", rows, 242_696_752)
    assert (s, e) == (242_559_540, 242_696_752)
    assert cl.check_locus("chr2", s, e, gain, 242_696_752) == (s, e)
    with pytest.raises(ValueError):
        cl.check_locus("chr2", 242_559_540, 242_697_000, gain, 242_696_752)


'''
open(t, "w").write(s.replace(anchor, test + anchor))
print("patched")

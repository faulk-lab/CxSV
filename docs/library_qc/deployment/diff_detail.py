#!/usr/bin/env python3
"""For files compare_qc_masked.py calls DIFFERENT: VCF bodies (non-header lines) and header lines that
differ; PDFs after removing CreationDate/ModDate; .snf reported only by size. Usage: diff_detail.py OLD NEW COMPARE_TSV"""
import gzip, re, sys, os
old, new, tsv = sys.argv[1:4]
def lines(p):
    op = gzip.open if p.endswith('.gz') else open
    with op(p, 'rt', errors='replace') as fh:
        return fh.read().splitlines()
for row in open(tsv).read().splitlines()[1:]:
    rel, st = row.split('\t')
    if st != 'DIFFERENT':
        continue
    po, pn = os.path.join(old, rel), os.path.join(new, rel)
    if '.vcf' in rel:
        lo, ln = lines(po), lines(pn)
        bo = [l for l in lo if not l.startswith('#')]; bn = [l for l in ln if not l.startswith('#')]
        ho = [l for l in lo if l.startswith('#')]; hn = [l for l in ln if l.startswith('#')]
        dh = [(a, b) for a, b in zip(ho, hn) if a != b]
        print(f"{rel}: body records {len(bo)} vs {len(bn)}, body identical {bo == bn}; header lines {len(ho)} vs {len(hn)}, differing {len(dh) + abs(len(ho)-len(hn))}")
        for a, b in dh[:3]:
            print("   old:", a[:160]); print("   new:", b[:160])
    elif rel.endswith('.pdf'):
        strip = lambda b: re.sub(rb'/(CreationDate|ModDate)\s*\([^)]*\)', b'', b)
        print(f"{rel}: identical after removing CreationDate/ModDate: {strip(open(po,'rb').read()) == strip(open(pn,'rb').read())}")
    else:
        print(f"{rel}: sizes {os.path.getsize(po)} vs {os.path.getsize(pn)}")

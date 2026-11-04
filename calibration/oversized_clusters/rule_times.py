#!/usr/bin/env python3
"""Wall time of every job in a Snakemake log: '[timestamp]' + 'rule NAME:' ... 'Finished jobid: N (Rule: NAME)'.
Usage: rule_times.py SNAKEMAKE.log > rule_times.tsv (rule, jobs, first_start, last_end, max_job_h, sum_job_h)"""
import re, sys
from datetime import datetime
from collections import defaultdict
ts, start, rule_of, cur_rule = None, {}, {}, None
jobs = defaultdict(list)
for line in open(sys.argv[1], errors="replace"):
    m = re.match(r"^\[(\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})\]", line)
    if m:
        ts = datetime.strptime(" ".join(m.group(1).split()), "%a %b %d %H:%M:%S %Y"); continue
    m = re.match(r"^(?:local)?rule (\S+):", line)
    if m:
        cur_rule = m.group(1); continue
    m = re.match(r"^\s+jobid: (\d+)", line)
    if m and cur_rule:
        start[m.group(1)] = ts; rule_of[m.group(1)] = cur_rule; cur_rule = None; continue
    m = re.match(r"^Finished job(?:id)?:? (\d+)", line.strip())
    if m and m.group(1) in start:
        jobs[rule_of[m.group(1)]].append((start[m.group(1)], ts))
print("rule\tjobs\tfirst_start\tlast_end\tmax_job_h\tsum_job_h")
for r, v in sorted(jobs.items(), key=lambda kv: min(x[0] for x in kv[1])):
    d = [(e - s).total_seconds() / 3600 for s, e in v]
    print(f"{r}\t{len(v)}\t{min(s for s, _ in v)}\t{max(e for _, e in v)}\t{max(d):.2f}\t{sum(d):.2f}")

#!/usr/bin/env python3
"""Simulate heterozygous reads for every planted event and splice them into the host.

Inputs come from the earlier stages:
  design_events.py     -> events.yaml
  build_haplotypes.py  -> <build>/<EVENT>/{alt.fa, ref.fa, replace.bed}
  host realignment     -> host CRAM (minimap2 2.24, -x map-ont -a --eqx -K 4g, hs1)
  badread training     -> <model>/{host_error_model, host_qscore_model, host_read_stats.tsv}
                          (identity_max_sim caps the beta maximum; the host's p99.5
                          identity is 100%, which would make badread emit perfect reads)

For each event:
  1. Local depth = median host depth over the event's replace windows.
  2. Badread simulates the derivative haplotype(s) (alt.fa) and the matching
     reference windows (ref.fa) at local_depth / ploidy each, so the locus stays
     at its original depth with one derivative and one reference copy. On a male
     host chrX/chrY are hemizygous: the derivative is simulated alone at full depth.
  3. Reads are renamed sim_<EVENT>_<alt|ref>_<n>; their origin is kept in
     sim_read_origins.tsv (badread's header comment).

Then all simulated reads are aligned exactly like the panel, and the host is spliced:
  - every host read with ANY alignment overlapping a replace window is dropped
    (all of its records, genome-wide),
  - every simulated read with ANY alignment overlapping a replace window is kept
    (all of its records), everything else simulated is discarded,
  - host-minus-dropped and kept-simulated are merged into one coordinate-sorted CRAM.
Simulated reads carry RG ID:sim (same SM as the host) and have no MM/ML tags.

Usage:
  simulate_and_splice.py --events events.yaml --build BUILD --host HOST.cram \
      --ref hs1.fa --index hs1.mmi --model MODEL_DIR --out OUT_DIR \
      [--sex auto|male|female] [--threads 16] [--seed 3559]
"""
import argparse
import os
import statistics
import subprocess
import sys

import yaml

SM = "pos_control_3559"


def sh(cmd, **kw):
    print("+", cmd, flush=True)
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", **kw)


def sh_out(cmd):
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash",
                          capture_output=True, text=True).stdout


def read_bed(path):
    rows = []
    for line in open(path):
        if line.strip():
            c, s, e = line.split("\t")[:3]
            rows.append((c, int(s), int(e)))
    return rows


def median_depth(cram, ref, chrom, start, end):
    out = sh_out(f"samtools depth -a --reference {ref} -r {chrom}:{start + 1}-{end} {cram} | cut -f3")
    vals = [int(v) for v in out.split()]
    return statistics.median(vals) if vals else 0.0


def host_sex(cram, ref):
    """Male if chrX mean depth is < 0.7 x the autosomal mean depth.

    Uses `samtools coverage` columns (#rname ... meandepth); reads the host's
    precomputed <cram-stem>.../coverage.tsv when present, else runs it.
    """
    cov = os.path.join(os.path.dirname(cram), "coverage.tsv")
    text = open(cov).read() if os.path.exists(cov) else sh_out(f"samtools coverage --reference {ref} {cram}")
    depth, length = {}, {}
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        depth[f[0]], length[f[0]] = float(f[6]), int(f[2]) - int(f[1]) + 1
    auto = [f"chr{i}" for i in range(1, 23)]
    a = sum(depth[c] * length[c] for c in auto) / sum(length[c] for c in auto)
    ratio = depth["chrX"] / a
    return ("male" if ratio < 0.7 else "female"), ratio


def read_stats(path):
    d = {}
    for line in open(path):
        k, v = line.rstrip("\n").split("\t")
        d[k] = v
    return d


def simulate(fasta, depth, model, stats, seed, prefix, hap, out_fq, origins):
    """Run badread on one FASTA at the given depth; rename reads; record origins."""
    if depth <= 0:
        return 0
    tmp = out_fq + ".raw.fq"
    sh(f"badread simulate --reference {fasta} --quantity {depth:.3f}x "
       f"--error_model {model}/host_error_model --qscore_model {model}/host_qscore_model "
       f"--length {stats['mean_len']},{stats['sd_len']} "
       f"--identity {stats['identity_mean']},{stats.get('identity_max_sim', stats['identity_p99.5'])},{stats['identity_sd']} "
       f"--junk_reads 0 --random_reads 0 --chimeras 0 --glitches 0,0,0 "
       f"--start_adapter 0,0 --end_adapter 0,0 --seed {seed} > {tmp} 2> {tmp}.log")
    n = 0
    with open(tmp) as fh, open(out_fq, "a") as out:
        while True:
            h = fh.readline()
            if not h:
                break
            seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
            n += 1
            name = f"sim_{prefix}_{hap}_{n:06d}"
            origins.write(f"{name}\t{h[1:].strip()}\n")
            out.write(f"@{name}\n{seq}+\n{qual}")
    os.remove(tmp)
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--events", required=True)
    ap.add_argument("--build", required=True)
    ap.add_argument("--host", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sex", default="auto", choices=["auto", "male", "female"])
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--seed", type=int, default=3559)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    events = yaml.safe_load(open(a.events))["events"]
    stats = read_stats(os.path.join(a.model, "host_read_stats.tsv"))
    sex, ratio = host_sex(a.host, a.ref)
    if a.sex != "auto":
        sex = a.sex
    print(f"host sex: {sex} (chrX/autosome mean-depth ratio {ratio:.2f})", flush=True)

    sim_fq = os.path.join(a.out, "sim_reads.fq")
    if os.path.exists(sim_fq):
        os.remove(sim_fq)
    all_bed = os.path.join(a.out, "replace_all.bed")
    plan = open(os.path.join(a.out, "simulation_plan.tsv"), "w")
    plan.write("event\tchroms\tlocal_depth\tploidy\talt_depth\tref_depth\talt_reads\tref_reads\n")
    with open(os.path.join(a.out, "sim_read_origins.tsv"), "w") as origins, open(all_bed, "w") as bedout:
        for i, ev in enumerate(events):
            eid = ev["id"]
            windows = read_bed(os.path.join(a.build, eid, "replace.bed"))
            for c, s, e in windows:
                bedout.write(f"{c}\t{s}\t{e}\t{eid}\n")
            depth = statistics.median(median_depth(a.host, a.ref, c, s, e) for c, s, e in windows)
            chroms = sorted({c for c, _, _ in windows})
            hemizygous = sex == "male" and all(c in ("chrX", "chrY") for c in chroms)
            ploidy = 1 if hemizygous else 2
            alt_depth = depth / ploidy
            ref_depth = depth - alt_depth
            seed = a.seed * 100 + i
            n_alt = simulate(os.path.join(a.build, eid, "alt.fa"), alt_depth, a.model, stats,
                             seed, eid, "alt", sim_fq, origins)
            n_ref = simulate(os.path.join(a.build, eid, "ref.fa"), ref_depth, a.model, stats,
                             seed + 50, eid, "ref", sim_fq, origins)
            plan.write(f"{eid}\t{','.join(chroms)}\t{depth:.1f}\t{ploidy}\t{alt_depth:.2f}\t"
                       f"{ref_depth:.2f}\t{n_alt}\t{n_ref}\n")
            plan.flush()
    plan.close()

    sim_bam = os.path.join(a.out, "sim_reads.hs1.minimap2.bam")
    sh(f"minimap2 -t {a.threads} -x map-ont -a --eqx -K 4g "
       f"-R '@RG\\tID:sim\\tSM:{SM}\\tDS:badread-simulated-positive-control-reads' "
       f"{a.index} {sim_fq} 2> {a.out}/sim_minimap2.stderr "
       f"| samtools sort -@ 4 -m 2G -o {sim_bam} - && samtools index {sim_bam}")

    drop = os.path.join(a.out, "host_reads_dropped.txt")
    keep = os.path.join(a.out, "sim_reads_kept.txt")
    sh(f"samtools view --reference {a.ref} -L {all_bed} {a.host} | cut -f1 | sort -u > {drop}")
    sh(f"samtools view -L {all_bed} {sim_bam} | cut -f1 | sort -u > {keep}")
    sim_keep = os.path.join(a.out, "sim_reads_kept.bam")
    sh(f"samtools view -b -N {keep} -o {sim_keep} {sim_bam} && samtools index {sim_keep}")

    final = os.path.join(a.out, "pos_control_3559.hs1.minimap2.cram")
    sh(f"samtools view -@ 4 -u --reference {a.ref} -N ^{drop} {a.host} "
       f"| samtools merge -@ {a.threads} --reference {a.ref} -O cram -o {final} - {sim_keep}")
    sh(f"samtools index {final}")
    print(f"dropped host reads: {sum(1 for _ in open(drop))}; kept simulated reads: "
          f"{sum(1 for _ in open(keep))}", flush=True)


if __name__ == "__main__":
    main()

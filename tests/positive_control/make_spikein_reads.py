#!/usr/bin/env python3
"""Simulate the heterozygous reads of a spike-in event (E13) for several host genomes.

A spike-in event is planted into the positive control AND into some background
genomes at regression time (docs/POSITIVE_CONTROL.md, "E13"). This script makes
the fixed, checksummed read sets that the regression splices in, one per host,
with exactly the method of simulate_and_splice.py:

  1. local depth = median depth of the host CRAM over the event's replace windows;
  2. badread simulates the derivative (alt.fa) and the reference windows (ref.fa)
     at depth/2 each (autosomal heterozygous; male chrX/chrY hemizygous as in
     simulate_and_splice.py), seed = --seed * 100 + host index (+50 for ref);
  3. reads are aligned like the panel (minimap2 2.24, -x map-ont -a --eqx -K 4g,
     hs1 k17 index) with read group ID:sim_<EVENT>, SM = the host's own SM;
  4. only reads with an alignment overlapping a replace window are kept.

Outputs in --out:
  <host>.<EVENT>.sim.bam(.bai)   reads to splice into that host
  simulation_plan.tsv            host, cram, sm, local depth, ploidy, alt/ref depth, reads kept
  sim_read_origins.<host>.tsv    badread header of every simulated read

Splicing (tests/lab/pos_control_regression.sh, splice_spikein): every host read
with an alignment overlapping a replace window is dropped, and the kept
simulated reads are merged in.

Usage:
  make_spikein_reads.py --build BUILD/E13 --event E13 --hosts hosts.tsv --ref hs1.fa \\
      --index hs1.map-ont.k17.mmi --model MODEL_DIR --out OUT [--threads 16] [--seed 3559]
hosts.tsv: name<TAB>cram<TAB>sex (male|female)
"""
import argparse
import os
import statistics
import subprocess

from simulate_and_splice import median_depth, read_bed, read_stats, sh, sh_out, simulate


def host_sm(cram, ref):
    for line in sh_out(f"samtools view -H --reference {ref} {cram}").splitlines():
        if line.startswith("@RG"):
            for f in line.split("\t"):
                if f.startswith("SM:"):
                    return f[3:]
    return os.path.basename(cram).split(".")[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", required=True, help="build_haplotypes.py output for the event (alt.fa, ref.fa, replace.bed)")
    ap.add_argument("--event", required=True)
    ap.add_argument("--hosts", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--seed", type=int, default=3559)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    stats = read_stats(os.path.join(a.model, "host_read_stats.tsv"))
    windows = read_bed(os.path.join(a.build, "replace.bed"))
    bed = os.path.join(a.out, f"{a.event}.replace.bed")
    with open(bed, "w") as fh:
        for c, s, e in windows:
            fh.write(f"{c}\t{s}\t{e}\t{a.event}\n")
    hosts = [l.rstrip("\n").split("\t") for l in open(a.hosts) if l.strip() and not l.startswith("#")]
    plan = open(os.path.join(a.out, "simulation_plan.tsv"), "w")
    plan.write("host\tcram\tsm\tlocal_depth\tploidy\talt_depth\tref_depth\talt_reads\tref_reads\treads_kept\n")
    for k, (name, cram, sex) in enumerate(hosts):
        sm = host_sm(cram, a.ref)
        depth = statistics.median(median_depth(cram, a.ref, c, s, e) for c, s, e in windows)
        chroms = {c for c, _, _ in windows}
        ploidy = 1 if (sex == "male" and chroms <= {"chrX", "chrY"}) else 2
        alt_depth = depth / ploidy
        ref_depth = depth - alt_depth
        fq = os.path.join(a.out, f"{name}.fq")
        if os.path.exists(fq):
            os.remove(fq)
        seed = a.seed * 100 + k
        with open(os.path.join(a.out, f"sim_read_origins.{name}.tsv"), "w") as origins:
            n_alt = simulate(os.path.join(a.build, "alt.fa"), alt_depth, a.model, stats, seed,
                             f"{a.event}_{name}", "alt", fq, origins)
            n_ref = simulate(os.path.join(a.build, "ref.fa"), ref_depth, a.model, stats, seed + 50,
                             f"{a.event}_{name}", "ref", fq, origins)
        allb = os.path.join(a.out, f"{name}.all.bam")
        sh(f"minimap2 -t {a.threads} -x map-ont -a --eqx -K 4g "
           f"-R '@RG\\tID:sim_{a.event}\\tSM:{sm}\\tDS:badread-simulated-spike-in-{a.event}' "
           f"{a.index} {fq} 2> {a.out}/{name}.minimap2.stderr "
           f"| samtools sort -@ 4 -m 2G -o {allb} - && samtools index {allb}")
        keep = os.path.join(a.out, f"{name}.kept.txt")
        sh(f"samtools view -L {bed} {allb} | cut -f1 | sort -u > {keep}")
        out_bam = os.path.join(a.out, f"{name}.{a.event}.sim.bam")
        sh(f"samtools view -b -N {keep} -o {out_bam} {allb} && samtools index {out_bam}")
        n_keep = sum(1 for _ in open(keep))
        plan.write(f"{name}\t{cram}\t{sm}\t{depth:.1f}\t{ploidy}\t{alt_depth:.2f}\t{ref_depth:.2f}\t"
                   f"{n_alt}\t{n_ref}\t{n_keep}\n")
        plan.flush()
        for f in (fq, allb, allb + ".bai", keep):
            os.remove(f)
    plan.close()
    subprocess.run(f"cd {a.out} && md5sum *.sim.bam *.sim.bam.bai {a.event}.replace.bed simulation_plan.tsv "
                   f"> MD5SUMS", shell=True, check=True, executable="/bin/bash")


if __name__ == "__main__":
    main()

#!/usr/bin/env bash
# Aggregate-only rerun of a panel build WITHOUT the flagged libraries (scratch, test
# release; docs/LIBRARY_QC.md). Per-sample outputs and held aggregate inputs
# (Sniffles2 population, callability masks) are hard links of OLD's files.
#   exclusion_rerun.sh OLD_BUILD EXC CODE FLAGGED_TSV SIF CORES [--dry-run]
# FLAGGED_TSV: cohort, sample_id (header). Writes only under EXC.
set -euo pipefail
OLD=$1 EXC=$2 CODE=$3 FLAG=$4 SIF=$5 CORES=$6 DRY=${7:-}
T=$CODE/scripts/panel_build
FORCE="panel_sample_manifest cohort_cxsv_tables build_panel build_panel_sites panel_callability_masks panel_site_callability calling_fingerprint panel_release_id panel_release_metadata panel_release"
mkdir -p "$EXC/sheets" "$EXC/provenance" "$EXC/wd"
python3 - "$OLD" "$EXC" "$FLAG" <<'PY'
import csv, json, os, sys
old, exc, flag = sys.argv[1:4]
drop = {(r["cohort"], r["sample_id"]) for r in csv.DictReader(open(flag), delimiter="\t")}
cfg = json.load(open(f"{old}/aggregate.yaml"))
n = 0
for c, spec in cfg["panel_cohorts"].items():
    rows = list(csv.DictReader(open(spec["samples"]), delimiter="\t"))
    keep = [r for r in rows if (c, r["sample_id"]) not in drop]
    n += len(rows) - len(keep)
    out = f"{exc}/sheets/{c}.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(keep)
    spec["samples"] = out
assert n == len(drop), (n, len(drop))
cfg["release"]["test_release"] = True
cfg["release"]["description"] = (f"SCRATCH test release, not for use: {os.path.basename(old)} aggregate rerun "
                                  f"without {n} library-QC-flagged libraries (docs/LIBRARY_QC.md)")
for name, root in (("aggregate.dry.yaml", f"{old}/results/panel"), ("aggregate.yaml", f"{exc}/results/panel")):
    cfg["panel_results_root"] = root
    json.dump(cfg, open(f"{exc}/{name}", "w"), indent=2)
print(f"dropped {n} libraries")
PY
X="singularity exec --containall --cleanenv -B /mnt/sargasso1 -B $EXC/wd:/tmp $SIF"
SM="snakemake -s $CODE/workflow/panel/Snakefile --cores $CORES --nolock"
if [ ! -s "$EXC/provenance/reused_upstream.tsv" ]; then
  $X $SM -n --directory $EXC/wd/dry --configfile $EXC/aggregate.dry.yaml --forcerun $FORCE \
      --allowed-rules $FORCE -- $OLD/results/panel/release/release_manifest.tsv > $EXC/provenance/dryrun.txt 2>&1 || { tail -30 $EXC/provenance/dryrun.txt; exit 1; }
  awk '/^Job stats/{f=1} f{print} f && /^$/{exit}' $EXC/provenance/dryrun.txt
  python3 $T/link_upstream.py --old $OLD/results/panel --new $EXC/results/panel --dryrun $EXC/provenance/dryrun.txt \
      --inventory $EXC/provenance/reused_upstream.tsv --include-log svim.log --apply
fi
[ "$DRY" = --dry-run ] && exit 0
$X $SM --directory $EXC/wd/run --configfile $EXC/aggregate.yaml --forcerun $FORCE --allowed-rules $FORCE \
    --keep-going --printshellcmds -- $EXC/results/panel/release/release_manifest.tsv > $EXC/provenance/run.log 2>&1
python3 $T/link_upstream.py --verify $EXC/provenance/reused_upstream.tsv

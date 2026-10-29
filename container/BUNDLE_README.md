# CxSV offline bundle @BUNDLE_ID@

Everything needed to run the CxSV patient scorer (and, if needed, the panel
builder) on a machine with **no network access**. Nothing outside this
directory is used. Built by `scripts/make_bundle.sh` from commit
`@COMMIT@` of https://github.com/faulk-lab/CxSV; full documentation is in
`repo/docs/RELEASE_AND_BUNDLE.md`.

| | |
|---|---|
| code | `repo/` (git archive of @COMMIT@); also `source/*.tar.gz` and `source/*.bundle` (git bundle) |
| container | `container/cxsv-pipeline.sif` (Singularity/Apptainer), `container/cxsv-pipeline.docker.tar.gz` (Docker); image `@IMAGE@`, id `@IMAGE_ID@`, SIF sha256 `@SIF_SHA256@` |
| reference | `reference/hs1.fa` (T2T-CHM13v2, UCSC names) + `.fai` + `ref_cache/` (htslib MD5 cache for CRAM decoding) |
| annotation | `resources/*.bed` |
| panel release | `panel_release/@RELEASE_ID@/` (test release: @TEST_RELEASE@) |
| catalogues | `catalogues/<name>/<version>/` junction-site tables with their manifests |
| integrity | `MANIFEST.sha256`, `bundle_manifest.tsv`, `bundle.json` |

**If "test release" above says `true`, the panel is a TEST panel (a few
genomes, sliced to small regions). Its public/private labels must not be used
for any real patient.**

## 1. Check the bundle

    ./run_cxsv.sh verify-bundle       # sha256 of every file
    ./run_cxsv.sh verify-release      # the panel release against its own manifest
    ./run_cxsv.sh netcheck            # the container sees loopback only

## 2. Run a patient

Write a run config, e.g. `family_042.yaml`:

    patient_runs:
      family_042:
        samples: /data/family_042/sheet.tsv   # sample_id, alignment[, index][, individual_id][, platform]
    threads: {samtools: 4, sv_callers: 16}

Alignments must be against hs1 (CRAM or BAM, indexed). Then

    ./run_cxsv.sh patient --workdir /scratch/family_042 --config family_042.yaml \
        --bind /data/family_042 --cores 16

`--bind` mounts an input directory read-only at the same path; the work
directory receives `results/patient/`. The container runtime is found
automatically (singularity, apptainer, docker; or `--runtime NAME`). The run
stops if the panel release does not verify against its manifest or if the
calling fingerprint differs from the panel's.

Main outputs per sample, under `results/patient/samples/<run>/<sample>/private/`:
`cxsv_panel_status.tsv` (every CxSV locus: Panel_Status, Panel_Carriers,
Panel_Callable = panel individuals callable at all its junctions,
Panel_Individuals), `private_cxsv_master_table.tsv`, `private_cxsv.vcf`,
`junction_panel_status.tsv`.

## 3. Without the run script

    singularity exec --containall --net --network none \
        --bind "$PWD:/bundle:ro" --bind /scratch/family_042:/work --bind /data/family_042:/data/family_042:ro \
        container/cxsv-pipeline.sif \
        snakemake -s /bundle/repo/workflow/patient/Snakefile --directory /work --cores 16 \
            --configfile /bundle/config/bundle_site.yaml /work/config/family_042.yaml

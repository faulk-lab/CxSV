# =============================================================================
# CxSV pipeline image: ONE image with every tool both workflows need, installed
# package-for-package from envs/cxsv-pipeline.explicit.txt (conda package URLs
# + MD5s of the env that built the panel). Runs under Docker, Singularity CE
# and Apptainer (no entrypoint, no daemon, any uid; tools on PATH).
#
# The image holds the environment only. The workflow code comes from the
# repository at a fixed commit (in the offline bundle: repo/, a git archive),
# so one image serves every commit that keeps the same lock.
#
# Build (context = repository root; scripts/make_bundle.sh does this from a
# git archive of the commit):
#   docker build -f container/Dockerfile \
#       --build-arg ENV_LOCK_SHA256=$(sha256sum envs/cxsv-pipeline.explicit.txt | cut -c1-64) \
#       -t cxsv-pipeline:<tag> .
# SIF:  singularity build cxsv-pipeline.sif docker-archive://cxsv-pipeline.docker.tar
# Run:  docs/RELEASE_AND_BUNDLE.md ("Running the workflows in the container").
# =============================================================================

# micromamba 2.10.0 on Debian 12 slim, pinned by digest (linux/amd64).
FROM mambaorg/micromamba:2.10.0-debian12-slim@sha256:bd7f1077166c1c71e9375a0e54835f151d5f216b96301b75f888c53c46d8e78b

ARG ENV_LOCK_SHA256=unknown
ARG IMAGE_TAG=cxsv-pipeline

USER root
COPY envs/cxsv-pipeline.explicit.txt /opt/cxsv-env/cxsv-pipeline.explicit.txt

# Explicit install: exactly the listed package files (checked against their
# MD5s), no solver. Then drop the package cache.
RUN micromamba create -y -p /opt/cxsv --file /opt/cxsv-env/cxsv-pipeline.explicit.txt \
 && micromamba clean -a -y \
 && rm -rf /opt/conda/pkgs \
 && echo "$ENV_LOCK_SHA256  /opt/cxsv-env/cxsv-pipeline.explicit.txt" > /opt/cxsv-env/SHA256 \
 && chmod -R a+rX /opt/cxsv /opt/cxsv-env

# Tools on PATH without activation (Singularity/Apptainer ignore entrypoints).
# C.UTF-8: one fixed locale for sort/awk. No bytecode writes into read-only
# mounts. CXSV_* are recorded in panel release.json / patient run records.
ENV PATH=/opt/cxsv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    CONDA_PREFIX=/opt/cxsv \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONNOUSERSITE=1 \
    XDG_CACHE_HOME=/tmp/.cache \
    CXSV_IMAGE=${IMAGE_TAG} \
    CXSV_ENV_LOCK_SHA256=${ENV_LOCK_SHA256}

LABEL org.opencontainers.image.title="cxsv-pipeline" \
      org.opencontainers.image.description="CxSV panel builder + patient scorer environment (envs/cxsv-pipeline.explicit.txt)" \
      org.opencontainers.image.source="https://github.com/faulk-lab/CxSV" \
      cxsv.env_lock_sha256="${ENV_LOCK_SHA256}"

# Smoke test at build time: every fingerprinted tool and the workflow engine.
RUN snakemake --version && sniffles --version && cuteSV --version && svim --version \
 && samtools --version | head -1 && bcftools --version | head -1 && bedtools --version \
 && mosdepth --version && (SURVIVOR 2>&1 | grep -i version) && Rscript -e 'library(ggplot2)' >/dev/null \
 && rm -rf /tmp/.cache /tmp/* && chmod 1777 /tmp

WORKDIR /work
ENTRYPOINT []
CMD ["bash"]

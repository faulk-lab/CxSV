#!/usr/bin/env Rscript
# plot_chromosome_density.R
# ──────────────────────────────────────────────────────────────────────────────
# Chromosome-level SV density plot (SVs per Mb).
# Uses CHM13v2.0 chromosome lengths; adjust chr_lengths if using hg38.
#
# Usage: Rscript plot_chromosome_density.R <long_read_filtered.vcf> <output.pdf>

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

this_file <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(this_file)), "plot_utils.R"))

args    <- commandArgs(trailingOnly = TRUE)
vcf_in  <- args[1]
outfile <- args[2]

# ── CHM13v2.0 / hs1 chromosome lengths (bp) ──────────────────────────────────
# Source: https://github.com/marbl/CHM13
chr_lengths <- tibble::tribble(
  ~Chr,     ~Length,
  "chr1",  248387328,
  "chr2",  242696752,
  "chr3",  201105948,
  "chr4",  193574945,
  "chr5",  182045439,
  "chr6",  172126628,
  "chr7",  160567428,
  "chr8",  146259331,
  "chr9",  150617247,
  "chr10", 134758134,
  "chr11", 135127769,
  "chr12", 133324548,
  "chr13", 113566686,
  "chr14", 101161492,
  "chr15",  99753195,
  "chr16",  96330374,
  "chr17",  84276897,
  "chr18",  80542538,
  "chr19",  61707364,
  "chr20",  66210255,
  "chr21",  45090682,
  "chr22",  51324926,
  "chrX",  154259566,
  "chrY",   62460029
)

# ── Parse VCF: extract CHROM and SVTYPE ───────────────────────────────────────
df <- parse_vcf_sv(vcf_in)
if (is.null(df)) stop("No variant records in VCF.")

# Keep only canonical autosomes + sex chromosomes
df <- df %>%
  filter(Chr %in% chr_lengths$Chr, SVTYPE %in% CANONICAL_SV_TYPES)

sv_per_chr <- df %>%
  count(Chr, name = "N_SV") %>%
  left_join(chr_lengths, by = "Chr") %>%
  mutate(
    SV_per_Mb = N_SV / (Length / 1e6),
    Chr        = factor(Chr, levels = chr_lengths$Chr)  # natural chr order
  ) %>%
  filter(!is.na(Length))

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(sv_per_chr, aes(x = Chr, y = SV_per_Mb, fill = SV_per_Mb)) +
  geom_col(color = "white", linewidth = 0.3, width = 0.8) +
  geom_text(aes(label = sprintf("%.1f", SV_per_Mb)),
            vjust = -0.4, size = 2.8, color = "grey30") +
  scale_fill_gradient(low = "#dce9f5", high = "#1d3557",
                      name = "SVs / Mb") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title    = "Chromosome-Level SV Density — Long Read",
    subtitle = "Normalized by chromosome length (CHM13v2.0 / hs1)",
    x        = "Chromosome",
    y        = "SVs per Megabase"
  ) +
  theme_cxsv(base_size = 12) +
  theme(
    axis.text.x     = element_text(angle = 50, hjust = 1, size = 8),
    legend.position = "right"
  )

ggsave(outfile, plot = p, width = 11, height = 5, device = "pdf")
message("Saved: ", outfile)

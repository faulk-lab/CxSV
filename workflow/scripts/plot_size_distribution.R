#!/usr/bin/env Rscript
# plot_size_distribution.R
# ──────────────────────────────────────────────────────────────────────────────
# Violin + boxplot of SV size on log10 scale, faceted by SV type.
# Usage: Rscript plot_size_distribution.R <merged.vcf> <output.pdf>

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

# ── Parse VCF for SVTYPE and SVLEN ───────────────────────────────────────────
df <- parse_vcf_sv(vcf_in)
if (is.null(df)) stop("No variant records found in VCF.")

df <- df %>%
  filter(!is.na(SVLEN), SVLEN > 0, SVTYPE %in% CANONICAL_SV_TYPES) %>%
  mutate(SVTYPE = factor(SVTYPE, levels = CANONICAL_SV_TYPES))

if (nrow(df) == 0) stop("No parseable SVLEN values found.")

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(x = SVTYPE, y = SVLEN, fill = SVTYPE)) +
  geom_violin(alpha = 0.6, color = "white", linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = 21, outlier.size = 1.2,
               outlier.alpha = 0.4, fill = "white", color = "grey30") +
  scale_y_log10(
    labels = label_number(scale_cut = cut_short_scale()),
    breaks = 10^(1:8),
    minor_breaks = NULL
  ) +
  scale_fill_manual(values = SV_TYPE_COLORS) +
  annotation_logticks(sides = "l", short = unit(0.1,"cm"),
                      mid   = unit(0.15,"cm"), long = unit(0.2,"cm")) +
  labs(
    title    = "SV Size Distribution by Type — Long Read",
    subtitle = "Log₁₀ scale; box shows IQR, whiskers = 1.5×IQR",
    x        = "SV Type",
    y        = "SV Length (bp)",
    fill     = "SV Type"
  ) +
  theme_cxsv(base_size = 13) +
  theme(legend.position = "none")

ggsave(outfile, plot = p, width = 8, height = 5.5, device = "pdf")
message("Saved: ", outfile)

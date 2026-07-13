#!/usr/bin/env Rscript
# plot_sv_count_by_type.R
# ──────────────────────────────────────────────────────────────────────────────
# Bar plot: SV count by type for long-read calls.
# Usage: Rscript plot_sv_count_by_type.R <sv_counts.tsv> <output.pdf>

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

this_file <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(this_file)), "plot_utils.R"))

args    <- commandArgs(trailingOnly = TRUE)
infile  <- args[1]
outfile <- args[2]

# ── Load data ─────────────────────────────────────────────────────────────────
df <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(df) <- c("Count", "SVTYPE")

# Normalize type names
df$SVTYPE <- toupper(trimws(df$SVTYPE))

# Keep canonical SV types; lump rare ones into "OTHER"
df$SVTYPE <- ifelse(df$SVTYPE %in% CANONICAL_SV_TYPES, df$SVTYPE, "OTHER")
df <- df %>%
  group_by(SVTYPE) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  arrange(desc(Count))

# Ordered factor for x-axis
df$SVTYPE <- factor(df$SVTYPE, levels = df$SVTYPE)

fill_vals <- SV_TYPE_COLORS[as.character(df$SVTYPE)]
fill_vals[is.na(fill_vals)] <- "#AAAAAA"

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(x = SVTYPE, y = Count, fill = SVTYPE)) +
  geom_col(color = "white", linewidth = 0.4, width = 0.7) +
  geom_text(aes(label = comma(Count)),
            vjust = -0.5, size = 3.4, fontface = "bold", color = "grey25") +
  scale_fill_manual(values = fill_vals) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Structural Variant Count by Type — Long Read",
    subtitle = "Ensemble call set (Sniffles2 + CuteSV + SVIM), SURVIVOR merged",
    x        = "SV Type",
    y        = "Number of SVs",
    fill     = "SV Type"
  ) +
  theme_cxsv(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(face = "bold", size = 12)
  )

ggsave(outfile, plot = p, width = 7, height = 5, device = "pdf")
message("Saved: ", outfile)

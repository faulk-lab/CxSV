#!/usr/bin/env Rscript
# plot_genomic_context.R
# ──────────────────────────────────────────────────────────────────────────────
# Bar plot of CxSV genomic context distribution.
# Reads the Genomic_Context column from cxsv_master_table.tsv, which is
# pre-computed by annotate_cxsv.py and is annotation-only — no loci are
# filtered by context.
#
# Expected context labels (from assign_genomic_context in annotate_cxsv.py):
#   Exonic  Promoter  Intronic  Repeat_intergenic  SegDup_intergenic  Intergenic
#
# The plot intentionally shows that most CxSVs are in intergenic / repeat /
# segdup regions — the grey zone where chromoanagenesis preferentially acts.
#
# Usage: Rscript plot_genomic_context.R <cxsv_master_table.tsv> <output.pdf>

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

df <- tryCatch(
  read.table(infile, header = TRUE, sep = "\t",
             stringsAsFactors = FALSE, quote = ""),
  error = function(e) stop(paste("Cannot read master table:", e$message))
)

if (nrow(df) == 0) {
  pdf(outfile, width = 7, height = 5)
  plot.new()
  text(0.5, 0.5, "No CxSV loci in master table.", cex = 1.4, col = "grey40")
  dev.off()
  quit(save = "no", status = 0)
}

# ── Validate column ───────────────────────────────────────────────────────────
if (!"Genomic_Context" %in% colnames(df)) {
  stop("'Genomic_Context' column not found in master table.\n",
       "Re-run annotate_cxsv.py to regenerate the master table.")
}

# ── Ordered factor — coding regions first, grey zone last ────────────────────
context_order <- c(
  "Exonic",
  "Promoter",
  "Intronic",
  "Repeat_intergenic",
  "SegDup_intergenic",
  "Intergenic"
)
# Any unexpected labels get lumped as "Other"
df$Genomic_Context[!df$Genomic_Context %in% context_order] <- "Other"
if ("Other" %in% df$Genomic_Context) {
  context_order <- c(context_order, "Other")
}
df$Genomic_Context <- factor(df$Genomic_Context, levels = context_order)

counts <- df %>%
  count(Genomic_Context, .drop = FALSE) %>%
  mutate(Pct = n / sum(n) * 100)

# ── Palette ───────────────────────────────────────────────────────────────────
# Warm = genic (coding risk)   Cool/grey = grey-zone (expected majority)
ctx_colors <- c(
  Exonic             = "#E63946",
  Promoter           = "#F4A261",
  Intronic           = "#457B9D",
  Repeat_intergenic  = "#6D6875",
  SegDup_intergenic  = "#2A9D8F",
  Intergenic         = "#A8DADC",
  Other              = "#CCCCCC"
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(counts, aes(x = Genomic_Context, y = n, fill = Genomic_Context)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = sprintf("%d\n(%.1f%%)", n, Pct)),
    vjust = -0.35, size = 3.3, fontface = "bold", color = "grey25",
    lineheight = 0.9
  ) +
  scale_fill_manual(values = ctx_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title    = "CxSV Genomic Context",
    subtitle = paste0(
      "n = ", nrow(df), " loci  |  ",
      "Annotation-only — no loci filtered by context\n",
      "Grey-zone enrichment (Repeat/SegDup/Intergenic) is expected for chromoanagenesis loci"
    ),
    x = "Genomic Context",
    y = "Number of CxSV Loci"
  ) +
  theme_cxsv(base_size = 13) +
  theme(
    plot.subtitle   = element_text(color = "grey35", size = 8.5, lineheight = 1.3),
    legend.position = "none",
    axis.text.x     = element_text(angle = 30, hjust = 1, size = 10, face = "bold")
  )

ggsave(outfile, plot = p, width = 8, height = 5.5, device = "pdf")
message("Saved: ", outfile)

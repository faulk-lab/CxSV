#!/usr/bin/env Rscript
# plot_complexity_class.R
# ──────────────────────────────────────────────────────────────────────────────
# Three-panel figure for the CxSV master table:
#   Panel A — Complexity class frequency (bar plot)
#   Panel B — SV type co-occurrence per complexity class (heatmap)
#   Panel C — CxSV criteria frequency (C1-C5 bar plot)
#
# Usage: Rscript plot_complexity_class.R <cxsv_master_table.tsv> <output.pdf>

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
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
  pdf(outfile, width = 9, height = 10)
  plot.new()
  text(0.5, 0.5, "No CxSV loci in master table.", cex = 1.4, col = "grey40")
  dev.off()
  quit(save = "no", status = 0)
}

# ── Colour palette ────────────────────────────────────────────────────────────
class_pal <- c(
  "DEL+INV"                  = "#E63946",
  "DUP+INV"                  = "#457B9D",
  "DEL+DUP+INV"              = "#9B5DE5",
  "INS+DEL"                  = "#E9C46A",
  "BND_cluster"              = "#F4A261",
  "Multi_breakpoint_cluster" = "#2A9D8F",
  "Chromothripsis_like"      = "#1D1D1D",
  "Nested_SV"                = "#6D6875",
  "Overlapping_types"        = "#A8DADC",
  "CN_plus_orientation_change" = "#264653",
  "Complex_other"            = "#AAAAAA"
)

# ── PANEL A: complexity class bar ─────────────────────────────────────────────
class_counts <- df %>%
  count(Cmplx_Class) %>%
  arrange(desc(n)) %>%
  mutate(Cmplx_Class = factor(Cmplx_Class, levels = Cmplx_Class))

pA <- ggplot(class_counts, aes(x = Cmplx_Class, y = n, fill = Cmplx_Class)) +
  geom_col(color = "white", linewidth = 0.35, width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.2, fontface = "bold",
            color = "grey25") +
  scale_fill_manual(values = class_pal, na.value = "#CCCCCC") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "A  CxSV Complexity Class Frequency",
    x = NULL, y = "Count"
  ) +
  theme_cxsv(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 38, hjust = 1, size = 8.5)
  )

# ── PANEL B: SV type co-occurrence heatmap ───────────────────────────────────
sv_types_exp <- df %>%
  select(Cmplx_Class, SV_Types) %>%
  mutate(SV_Types = strsplit(SV_Types, ",")) %>%
  unnest(SV_Types) %>%
  mutate(SV_Types = toupper(trimws(SV_Types))) %>%
  filter(SV_Types %in% CANONICAL_SV_TYPES)

heat_df <- sv_types_exp %>%
  count(Cmplx_Class, SV_Types) %>%
  group_by(Cmplx_Class) %>%
  mutate(Prop = n / sum(n)) %>%
  ungroup()

pB <- ggplot(heat_df, aes(x = SV_Types, y = Cmplx_Class, fill = Prop)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%", Prop * 100)),
            size = 2.8, color = "white", fontface = "bold") +
  scale_fill_gradient(low = "#dce9f5", high = "#1d3557",
                      labels = percent_format(accuracy = 1)) +
  labs(
    title = "B  SV Type Composition by Complexity Class",
    x = "SV Type", y = NULL, fill = "% within\nclass"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", size = 12),
    axis.text.y  = element_text(size = 8.5),
    panel.grid   = element_blank(),
    legend.key.height = unit(0.6, "cm")
  )

# ── PANEL C: CxSV criteria frequency ─────────────────────────────────────────
criteria_pal <- c(
  "C1_min_breakpoints"         = "#E63946",
  "C2_nested_SV"               = "#457B9D",
  "C3_overlapping_types"       = "#2A9D8F",
  "C4_CN_plus_orient"          = "#F4A261",
  "C5_chromoanagenesis_cluster"= "#9B5DE5"
)
criteria_labels <- c(
  "C1_min_breakpoints"         = "C1: >=3 breakpoints",
  "C2_nested_SV"               = "C2: Nested SV",
  "C3_overlapping_types"       = "C3: Overlapping types",
  "C4_CN_plus_orient"          = "C4: CN + orientation",
  "C5_chromoanagenesis_cluster"= "C5: Chromoanagenesis"
)

if ("CxSV_Criteria" %in% colnames(df)) {
  crit_df <- df %>%
    select(CxSV_Criteria) %>%
    filter(!is.na(CxSV_Criteria), CxSV_Criteria != "NONE", CxSV_Criteria != "") %>%
    mutate(Criteria = strsplit(CxSV_Criteria, ",")) %>%
    unnest(Criteria) %>%
    filter(Criteria %in% names(criteria_labels)) %>%
    count(Criteria) %>%
    mutate(
      Label    = criteria_labels[Criteria],
      Criteria = factor(Criteria, levels = names(criteria_labels))
    )
  
  pC <- ggplot(crit_df, aes(x = Criteria, y = n, fill = Criteria)) +
    geom_col(color = "white", linewidth = 0.35, width = 0.7) +
    geom_text(aes(label = n), vjust = -0.4, size = 3.2, fontface = "bold",
              color = "grey25") +
    scale_fill_manual(values = criteria_pal) +
    scale_x_discrete(labels = criteria_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title   = "C  CxSV Criteria Frequency",
      subtitle = "Loci can satisfy multiple criteria",
      x = NULL, y = "Loci matching criterion"
    ) +
    theme_cxsv(base_size = 11) +
    theme(
      plot.subtitle = element_text(color = "grey40", size = 8.5),
      legend.position = "none",
      axis.text.x   = element_text(angle = 30, hjust = 1, size = 9)
    )
} else {
  # Fallback if CxSV_Criteria column absent (old table)
  pC <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, size = 4, color = "grey50",
             label = "CxSV_Criteria column not found.\nRe-run classify_cxsv.py.") +
    theme_void() +
    labs(title = "C  CxSV Criteria Frequency")
}

# ── Combine ───────────────────────────────────────────────────────────────────
combined <- (pA / pB / pC) +
  plot_layout(heights = c(1, 1.1, 0.9)) +
  plot_annotation(
    title    = "CxSV Complexity & Criteria Summary",
    subtitle = paste0(
      "Long-read ensemble call set (Sniffles2 + CuteSV + SVIM), n = ",
      nrow(df), " CxSV loci\n",
      "Definition: >=3 bp (C1) | nested (C2) | overlapping types (C3) | ",
      "CN+orient (C4) | 10-50kb cluster (C5)"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey35", size = 9, lineheight = 1.3)
    )
  )

ggsave(outfile, plot = combined, width = 9, height = 13, device = "pdf")
message("Saved: ", outfile)

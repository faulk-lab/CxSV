#!/usr/bin/env Rscript
# plot_cxsv_venn.R
# ──────────────────────────────────────────────────────────────────────────────
# Three-way Venn diagram of CxSV caller concordance.
#
# Uses SURVIVOR SUPP_VEC from per-sample merged VCFs — NOT individual caller
# VCFs. SURVIVOR records which callers supported each merged call in the
# SUPP_VEC INFO field during the merge_sv_callers step.
#
# Bit order matches the filelist passed to SURVIVOR (merge_sv_callers rule):
#   bit 0 = CuteSV
#   bit 1 = SVIM
#   bit 2 = Sniffles2
# e.g. SUPP_VEC=101 means CuteSV + Sniffles2 called it, SVIM did not.
#
# Only calls that overlap a CxSV locus (IsCxSV==YES in cxsv_summary.tsv)
# are counted. Calls are pooled across all samples for a population-level view.
#
# Usage:
#   Rscript plot_cxsv_venn.R  "vcf1,vcf2,...,vcfN"  summary.tsv  output.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
    if (!requireNamespace("VennDiagram", quietly = TRUE)) {
      stop("Install ggVennDiagram:  install.packages('ggVennDiagram')")
    }
  }
})

this_file <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(this_file)), "plot_utils.R"))

args        <- commandArgs(trailingOnly = TRUE)
vcf_list    <- strsplit(args[1], ",")[[1]]
summary_tsv <- args[2]
outfile     <- args[3]

# ── Load CxSV loci from cxsv_summary.tsv ─────────────────────────────────────
cxsv_df <- tryCatch({
  df <- read.table(summary_tsv, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "")
  df[df$IsCxSV == "YES", ]
}, error = function(e) {
  stop(paste("Cannot read cxsv_summary.tsv:", e$message))
})

if (nrow(cxsv_df) == 0) {
  pdf(outfile, width = 7, height = 5)
  plot.new()
  text(0.5, 0.5, "No CxSV loci detected.\nVenn diagram cannot be produced.",
       cex = 1.4, col = "grey40", adj = c(0.5, 0.5))
  dev.off()
  message("No CxSV loci -- empty plot written to ", outfile)
  quit(save = "no", status = 0)
}

# ── Parse SUPP_VEC from each per-sample SURVIVOR merged VCF ──────────────────
# For each call overlapping a CxSV locus, assign to caller sets by bit position.
# Each call gets a unique ID = "sample::variantID" to avoid cross-sample collision.

cutesv_ids   <- character(0)
svim_ids     <- character(0)
sniffles_ids <- character(0)

n_vcfs_read  <- 0
n_calls_seen <- 0
n_cxsv_hits  <- 0

for (vcf_path in vcf_list) {
  vcf_path <- trimws(vcf_path)
  if (!file.exists(vcf_path)) {
    message("WARN: VCF not found, skipping: ", vcf_path)
    next
  }
  n_vcfs_read <- n_vcfs_read + 1
  sample_name <- sub("_merged_sv\\.vcf$", "", basename(vcf_path))

  for (ln in read_vcf_lines(vcf_path)) {
    parts <- strsplit(ln, "\t")[[1]]
    if (length(parts) < 8) next
    n_calls_seen <- n_calls_seen + 1

    chrom <- parts[1]
    pos   <- suppressWarnings(as.integer(parts[2]))
    vid   <- parts[3]
    info  <- parts[8]
    if (is.na(pos)) next

    # Parse END for overlap check
    end_field <- parse_info_field(info, "END")
    end       <- if (!is.na(end_field)) as.integer(end_field) else pos

    # Check if this call overlaps any CxSV locus
    is_cxsv_hit <- any(
      cxsv_df$Chr   == chrom &
        cxsv_df$Start <= end   &
        cxsv_df$End   >= pos
    )
    if (!is_cxsv_hit) next
    n_cxsv_hits <- n_cxsv_hits + 1

    # Parse SUPP_VEC
    vec <- parse_info_field(info, "SUPP_VEC")
    if (is.na(vec)) next

    call_id <- paste0(sample_name, "::", vid)
    
    # Bit 0 = CuteSV
    if (nchar(vec) >= 1 && substr(vec, 1, 1) == "1")
      cutesv_ids <- c(cutesv_ids, call_id)
    # Bit 1 = SVIM
    if (nchar(vec) >= 2 && substr(vec, 2, 2) == "1")
      svim_ids   <- c(svim_ids,   call_id)
    # Bit 2 = Sniffles2
    if (nchar(vec) >= 3 && substr(vec, 3, 3) == "1")
      sniffles_ids <- c(sniffles_ids, call_id)
  }
}

cutesv_ids   <- unique(cutesv_ids)
svim_ids     <- unique(svim_ids)
sniffles_ids <- unique(sniffles_ids)

message(sprintf(
  "Venn input: %d VCFs | %d calls total | %d CxSV-overlapping",
  n_vcfs_read, n_calls_seen, n_cxsv_hits
))
message(sprintf(
  "Caller sets: CuteSV=%d  SVIM=%d  Sniffles2=%d",
  length(cutesv_ids), length(svim_ids), length(sniffles_ids)
))

venn_list <- list(
  CuteSV    = cutesv_ids,
  SVIM      = svim_ids,
  Sniffles2 = sniffles_ids
)

# ── Build subtitle and caption ────────────────────────────────────────────────
subtitle_txt <- sprintf(
  "CuteSV n=%d  |  SVIM n=%d  |  Sniffles2 n=%d  |  %d samples pooled",
  length(cutesv_ids), length(svim_ids), length(sniffles_ids), n_vcfs_read
)
caption_txt <- paste0(
  "Source: SURVIVOR SUPP_VEC field in per-sample merged VCFs.\n",
  "Only calls overlapping CxSV loci (cxsv_summary.tsv, IsCxSV=YES) are counted.\n",
  "Bit order: CuteSV=bit0, SVIM=bit1, Sniffles2=bit2.\n",
  "CxSV definition: >=3 bp (C1) | nested (C2) | overlapping types (C3) | ",
  "CN+orient (C4) | 10-50kb cluster (C5).\n",
  "Intergenic/repeat/segdup loci are included -- most CxSVs are in grey-zone regions."
)

# ── Plot ──────────────────────────────────────────────────────────────────────
pdf(outfile, width = 8.5, height = 7.5)

if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
  library(ggVennDiagram)
  p <- ggVennDiagram(
    venn_list,
    label_alpha = 0,
    label       = "count"
  ) +
    scale_fill_gradient(low = "#dce9f5", high = "#1d3557") +
    scale_color_manual(values = rep("grey20", 3)) +
    labs(
      title    = "CxSV Caller Concordance — SURVIVOR SUPP_VEC",
      subtitle = subtitle_txt,
      caption  = caption_txt
    ) +
    theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(color = "grey30", size = 9,  hjust = 0.5),
      plot.caption  = element_text(color = "grey50", size = 7,  hjust = 0,
                                   lineheight = 1.3),
      legend.position = "none",
      plot.margin   = margin(10, 15, 10, 15)
    )
  print(p)
  
} else {
  library(VennDiagram)
  grid.newpage()
  v <- venn.diagram(
    x            = venn_list,
    filename     = NULL,
    fill         = unname(CALLER_COLORS[c("CuteSV", "SVIM", "Sniffles2")]),
    alpha        = 0.45,
    cex          = 1.5,
    cat.cex      = 1.3,
    cat.fontface = "bold",
    main         = "CxSV Caller Concordance (SURVIVOR SUPP_VEC)",
    main.cex     = 1.3,
    main.fontface = "bold",
    sub          = subtitle_txt,
    sub.cex      = 0.82
  )
  grid.draw(v)
}

dev.off()
message("Saved: ", outfile)

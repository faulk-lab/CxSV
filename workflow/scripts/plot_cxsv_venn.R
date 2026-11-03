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

# Vectorised per file: the earlier per-line loop compared every record with
# every CxSV locus row and grew the id vectors with c(), which does not finish
# on the full panel. Overlap is now a per-chromosome lookup over the loci
# sorted by Start with a running maximum of End: a record [pos, end] overlaps
# some locus iff the loci with Start <= end include one with End >= pos.
loci <- cxsv_df[!is.na(cxsv_df$Start) & !is.na(cxsv_df$End), c("Chr", "Start", "End")]
loci_index <- lapply(split(loci, loci$Chr), function(d) {
  o <- order(d$Start)
  list(start = d$Start[o], max_end = cummax(d$End[o]))
})
overlaps_cxsv <- function(chrom, pos, end) {
  hit <- logical(length(chrom))
  for (ch in intersect(unique(chrom), names(loci_index))) {
    sel <- which(chrom == ch)
    ix  <- loci_index[[ch]]
    k   <- findInterval(end[sel], ix$start)   # number of loci with Start <= end
    ok  <- k > 0L
    hit[sel[ok]] <- ix$max_end[k[ok]] >= pos[sel[ok]]
  }
  hit
}

cutesv_list   <- vector("list", length(vcf_list))
svim_list     <- vector("list", length(vcf_list))
sniffles_list <- vector("list", length(vcf_list))

n_vcfs_read  <- 0
n_calls_seen <- 0
n_cxsv_hits  <- 0
n_end_na     <- 0

for (i in seq_along(vcf_list)) {
  vcf_path <- trimws(vcf_list[i])
  if (!file.exists(vcf_path)) {
    message("WARN: VCF not found, skipping: ", vcf_path)
    next
  }
  n_vcfs_read <- n_vcfs_read + 1
  # unchanged from the per-line version (see the follow-up note in the commit):
  sample_name <- sub("_merged_sv\\.vcf$", "", basename(vcf_path))

  parts <- strsplit(read_vcf_lines(vcf_path), "\t", fixed = TRUE)
  parts <- parts[lengths(parts) >= 8]
  n_calls_seen <- n_calls_seen + length(parts)
  if (length(parts) == 0) next

  chrom <- vapply(parts, `[`, character(1), 1L)
  pos   <- suppressWarnings(as.integer(vapply(parts, `[`, character(1), 2L)))
  vid   <- vapply(parts, `[`, character(1), 3L)
  info  <- vapply(parts, `[`, character(1), 8L)
  keep  <- !is.na(pos)
  chrom <- chrom[keep]; pos <- pos[keep]; vid <- vid[keep]; info <- info[keep]

  # END for the overlap check (first END= match, as parse_info_field), else POS
  end <- pos
  m   <- regexpr("(?<=END=)[^;]+", info, perl = TRUE)
  end[m != -1L] <- suppressWarnings(as.integer(regmatches(info, m)))
  n_end_na <- n_end_na + sum(is.na(end))
  hit <- overlaps_cxsv(chrom, pos, end) & !is.na(end)
  n_cxsv_hits <- n_cxsv_hits + sum(hit)
  if (!any(hit)) next

  info <- info[hit]; vid <- vid[hit]
  vec  <- rep(NA_character_, length(info))
  m    <- regexpr("(?<=SUPP_VEC=)[^;]+", info, perl = TRUE)
  vec[m != -1L] <- regmatches(info, m)
  ok   <- !is.na(vec)
  vec  <- vec[ok]
  call_id <- paste0(sample_name, "::", vid[ok])

  # Bit 0 = CuteSV, bit 1 = SVIM, bit 2 = Sniffles2
  cutesv_list[[i]]   <- call_id[nchar(vec) >= 1 & substr(vec, 1, 1) == "1"]
  svim_list[[i]]     <- call_id[nchar(vec) >= 2 & substr(vec, 2, 2) == "1"]
  sniffles_list[[i]] <- call_id[nchar(vec) >= 3 & substr(vec, 3, 3) == "1"]
}
if (n_end_na > 0) message("WARN: ", n_end_na, " records with a non-integer END were not counted as CxSV hits")

cutesv_ids   <- unlist(cutesv_list,   use.names = FALSE)
svim_ids     <- unlist(svim_list,     use.names = FALSE)
sniffles_ids <- unlist(sniffles_list, use.names = FALSE)
if (is.null(cutesv_ids))   cutesv_ids   <- character(0)
if (is.null(svim_ids))     svim_ids     <- character(0)
if (is.null(sniffles_ids)) sniffles_ids <- character(0)

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

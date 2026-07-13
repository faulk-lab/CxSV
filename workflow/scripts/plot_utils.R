#!/usr/bin/env Rscript
# plot_utils.R
# ──────────────────────────────────────────────────────────────────────────────
# Shared helpers sourced by the plot_*.R scripts: VCF parsing, the SV-type /
# caller color palettes, and the common ggplot theme. Keeps this boilerplate
# in one place instead of copy-pasted across every plotting script.
#
# Usage in a plot script:
#   source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))), "plot_utils.R"))

suppressPackageStartupMessages({
  library(ggplot2)
})

# ── SV type / caller palettes ────────────────────────────────────────────────
CANONICAL_SV_TYPES <- c("DEL", "DUP", "INV", "INS", "BND", "TRA")

SV_TYPE_COLORS <- c(
  DEL   = "#E63946",  # bold red
  DUP   = "#457B9D",  # steel blue
  INV   = "#2A9D8F",  # teal
  INS   = "#E9C46A",  # amber
  BND   = "#9B5DE5",  # purple
  TRA   = "#F4A261",  # orange
  OTHER = "#A8DADC"   # light teal
)

CALLER_COLORS <- c(Sniffles2 = "#E63946", CuteSV = "#457B9D", SVIM = "#2A9D8F")

# ── VCF parsing ───────────────────────────────────────────────────────────────

#' Read a VCF and return only the data lines (header stripped).
read_vcf_lines <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines[!grepl("^#", lines)]
}

#' Extract one INFO subfield (e.g. "SVTYPE", "SVLEN", "END") from a VCF INFO
#' string. Returns NA if the key is absent — safer than sub()-based extraction,
#' which silently returns the whole string unchanged on a non-match.
parse_info_field <- function(info, key) {
  m <- regmatches(info, regexpr(paste0("(?<=", key, "=)[^;]+"), info, perl = TRUE))
  if (length(m) == 0) NA_character_ else m
}

#' Parse a VCF file into a data.frame with Chr, Pos, SVTYPE, SVLEN, End columns.
#' Returns NULL if the file has no variant records.
parse_vcf_sv <- function(path) {
  data_lines <- read_vcf_lines(path)
  if (length(data_lines) == 0) return(NULL)

  records <- lapply(data_lines, function(ln) {
    parts <- strsplit(ln, "\t")[[1]]
    info  <- parts[8]
    list(
      Chr    = parts[1],
      Pos    = suppressWarnings(as.integer(parts[2])),
      SVTYPE = parse_info_field(info, "SVTYPE"),
      SVLEN  = suppressWarnings(abs(as.integer(parse_info_field(info, "SVLEN")))),
      End    = suppressWarnings(as.integer(parse_info_field(info, "END")))
    )
  })

  data.frame(
    Chr    = vapply(records, `[[`, character(1), "Chr"),
    Pos    = vapply(records, `[[`, integer(1),   "Pos"),
    SVTYPE = vapply(records, `[[`, character(1), "SVTYPE"),
    SVLEN  = vapply(records, `[[`, integer(1),   "SVLEN"),
    End    = vapply(records, `[[`, integer(1),   "End"),
    stringsAsFactors = FALSE
  )
}

# ── Shared theme ──────────────────────────────────────────────────────────────

#' Common ggplot theme used across all CxSV plots. Callers can layer further
#' theme() tweaks on top (e.g. legend.position, axis.text.x rotation).
theme_cxsv <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title          = element_text(face = "bold", size = base_size + 1),
      plot.subtitle       = element_text(color = "grey40", size = base_size - 3),
      panel.grid.major.y  = element_line(color = "grey92", linewidth = 0.4)
    )
}

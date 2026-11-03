#!/usr/bin/env Rscript
# plot_breakpoint_resolution.R
# ──────────────────────────────────────────────────────────────────────────────
# Density plot of breakpoint confidence interval (CI) widths per caller.
# CI width = CIPOS_right - CIPOS_left  (from CIPOS/CIEND INFO fields).
# If CI is absent, width is recorded as 1 (exact breakpoint resolution).
#
# Usage: Rscript plot_breakpoint_resolution.R \
#            <sniffles1.vcf,sniffles2.vcf,...> \
#            <cutesv1.vcf,cutesv2.vcf,...> \
#            <svim1.vcf,svim2.vcf,...> \
#            <output.pdf>

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

this_file <- sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
source(file.path(dirname(normalizePath(this_file)), "plot_utils.R"))

args    <- commandArgs(trailingOnly = TRUE)
sniffles_files <- strsplit(args[1], ",")[[1]]
cutesv_files   <- strsplit(args[2], ",")[[1]]
svim_files     <- strsplit(args[3], ",")[[1]]
outfile        <- args[4]

# ── Helper: extract CI widths from a list of VCF files ───────────────────────
# Vectorised per file (the earlier per-line loop grew `widths` with c() once per
# record, which is quadratic: ~20M records for the 1,017-sample ONT cohort).
# Same rule as before: a record with >= 8 tab-separated fields contributes
# |lo| + |hi| from its first CIPOS=lo,hi match, or 1 when CIPOS is absent.
ci_widths <- function(vcf_paths, caller_name) {
  per_file <- vector("list", length(vcf_paths))
  for (i in seq_along(vcf_paths)) {
    path <- vcf_paths[i]
    if (!file.exists(path)) next
    lines <- read_vcf_lines(path)
    if (length(lines) == 0) next
    parts <- strsplit(lines, "\t", fixed = TRUE)
    parts <- parts[lengths(parts) >= 8]
    if (length(parts) == 0) next
    info <- vapply(parts, `[`, character(1), 8L)
    w    <- rep(1, length(info))   # no CI = exact
    m    <- regexpr("(?<=CIPOS=)-?[0-9]+,-?[0-9]+", info, perl = TRUE)
    has  <- m != -1L
    if (any(has)) {
      lohi <- strsplit(regmatches(info, m), ",", fixed = TRUE)
      lo   <- as.integer(vapply(lohi, `[`, character(1), 1L))
      hi   <- as.integer(vapply(lohi, `[`, character(1), 2L))
      w[has] <- abs(hi) + abs(lo)   # total CI width
    }
    per_file[[i]] <- w
  }
  widths <- as.numeric(unlist(per_file, use.names = FALSE))
  if (length(widths) == 0) return(NULL)
  message(sprintf("%s: %d records from %d VCFs", caller_name, length(widths), length(vcf_paths)))
  data.frame(CI_Width = widths, Caller = caller_name, stringsAsFactors = FALSE)
}

df_list <- list(
  ci_widths(sniffles_files, "Sniffles2"),
  ci_widths(cutesv_files,   "CuteSV"),
  ci_widths(svim_files,     "SVIM")
)
df <- bind_rows(Filter(Negate(is.null), df_list))

if (nrow(df) == 0) stop("No CI data extracted. Check VCF paths.")

df$Caller <- factor(df$Caller, levels = names(CALLER_COLORS))

# ── Compute medians for annotation ───────────────────────────────────────────
medians <- df %>%
  group_by(Caller) %>%
  summarise(med = median(CI_Width), .groups = "drop")

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(df, aes(x = CI_Width + 1, fill = Caller, color = Caller)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  geom_vline(data = medians, aes(xintercept = med + 1, color = Caller),
             linetype = "dashed", linewidth = 0.8) +
  geom_text(data = medians,
            aes(x = med + 1, y = Inf, label = paste0("med=", round(med)), color = Caller),
            vjust = 2, hjust = -0.1, size = 3.2, fontface = "italic", show.legend = FALSE) +
  scale_x_log10(labels = comma, breaks = c(1, 10, 100, 1000)) +
  scale_fill_manual(values = CALLER_COLORS) +
  scale_color_manual(values = CALLER_COLORS) +
  labs(
    title    = "Breakpoint Resolution Comparison — Long Read Callers",
    subtitle = "Distribution of CIPOS interval widths (bp); dashed = median",
    x        = "CI Width + 1 (bp, log₁₀ scale)",
    y        = "Density",
    fill     = "Caller",
    color    = "Caller"
  ) +
  theme_cxsv(base_size = 13) +
  theme(
    legend.position = c(0.85, 0.75),
    legend.background = element_rect(fill = "white", color = "grey80")
  )

ggsave(outfile, plot = p, width = 8, height = 5, device = "pdf")
message("Saved: ", outfile)

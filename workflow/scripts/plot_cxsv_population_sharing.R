#!/usr/bin/env Rscript
# plot_cxsv_population_sharing.R
# ──────────────────────────────────────────────────────────────────────────────
# STANDALONE — not part of `rule all` in the Snakefile. Run manually once you
# have a per-locus sample-sharing BED (Chr/Start/End/N_Samples), e.g. from
# `bedtools intersect -c` between cxsv_population.bed and per-sample CxSV BEDs.
#
# Three plots:
#   1. Histogram — frequency distribution of sample sharing
#   2. Genome map — where shared loci are located, colored by frequency
#   3. Stacked bar — sharing categories per chromosome
#
# Usage:
#   Rscript workflow/scripts/plot_cxsv_population_sharing.R \
#       cxsv_sample_counts.bed \
#       results/plots/cxsv_population_sharing.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

args    <- commandArgs(trailingOnly=TRUE)
infile  <- args[1]
outfile <- args[2]

# ── Load data ─────────────────────────────────────────────────────────────────
df <- read.table(infile, header=FALSE, sep="\t",
                 col.names=c("Chr","Start","End","N_Samples"),
                 stringsAsFactors=FALSE)

# Keep standard chromosomes in order
chr_order <- paste0("chr", c(1:22,"X","Y"))
df <- df %>%
  filter(Chr %in% chr_order) %>%
  mutate(
    Chr      = factor(Chr, levels=chr_order),
    Mid      = (Start + End) / 2,
    Size     = End - Start,
    # Sharing category
    Category = case_when(
      N_Samples >= 100 ~ "Fixed (100 samples)",
      N_Samples >= 80  ~ "Very common (80-99)",
      N_Samples >= 50  ~ "Common (50-79)",
      N_Samples >= 20  ~ "Low frequency (20-49)",
      TRUE             ~ "Rare (<20)"
    ),
    Category = factor(Category, levels=c(
      "Fixed (100 samples)", "Very common (80-99)",
      "Common (50-79)", "Low frequency (20-49)", "Rare (<20)"
    ))
  )

cat_colors <- c(
  "Fixed (100 samples)" = "#1d3557",
  "Very common (80-99)" = "#457b9d",
  "Common (50-79)"      = "#a8dadc",
  "Low frequency (20-49)"= "#f4a261",
  "Rare (<20)"          = "#e63946"
)

n_total  <- nrow(df)
n_fixed  <- sum(df$N_Samples >= 100)
n_common <- sum(df$N_Samples >= 50)
n_rare   <- sum(df$N_Samples < 10)

# ── CHM13 chromosome lengths (bp) ─────────────────────────────────────────────
chr_lengths <- c(
  chr1=248387328, chr2=242696752, chr3=201105948, chr4=193574945,
  chr5=182045439, chr6=172126628, chr7=160567428, chr8=146259331,
  chr9=150617247, chr10=134758134, chr11=135127769, chr12=133324548,
  chr13=113566686, chr14=101161492, chr15=99753195,  chr16=96330374,
  chr17=84276897,  chr18=80542538,  chr19=61707364,  chr20=66210255,
  chr21=45090682,  chr22=51324926,  chrX=154259566,  chrY=62460029
)

# ── PLOT 1: Histogram ─────────────────────────────────────────────────────────
p1 <- ggplot(df, aes(x=N_Samples, fill=Category)) +
  geom_histogram(binwidth=2, color="white", linewidth=0.2) +
  geom_vline(xintercept=c(50,80,100), linetype="dashed",
             color="grey30", linewidth=0.5) +
  annotate("text", x=51, y=Inf, label="50%", vjust=1.5, hjust=-0.1,
           size=3, color="grey30") +
  annotate("text", x=81, y=Inf, label="80%", vjust=1.5, hjust=-0.1,
           size=3, color="grey30") +
  annotate("text", x=101, y=Inf, label="100%", vjust=1.5, hjust=-0.1,
           size=3, color="grey30") +
  scale_fill_manual(values=cat_colors) +
  scale_x_continuous(breaks=seq(0,100,10),
                     labels=paste0(seq(0,100,10))) +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(
    title    = "A   CxSV Population Sharing — Frequency Distribution",
    subtitle = sprintf(
      "n=%s loci total  |  Fixed (100 samples): %s (%s%%)  |  Common (>=50): %s (%s%%)  |  Rare (<10): %s (%s%%)",
      formatC(n_total, format="d", big.mark=","),
      formatC(n_fixed, format="d", big.mark=","),
      round(100*n_fixed/n_total,1),
      formatC(n_common, format="d", big.mark=","),
      round(100*n_common/n_total,1),
      formatC(n_rare, format="d", big.mark=","),
      round(100*n_rare/n_total,1)
    ),
    x    = "Number of samples sharing locus (out of 100)",
    y    = "Number of CxSV loci",
    fill = "Sharing category"
  ) +
  theme_classic(base_size=12) +
  theme(
    plot.title    = element_text(face="bold", size=13),
    plot.subtitle = element_text(color="grey35", size=8.5, lineheight=1.3),
    legend.position = "bottom",
    legend.title  = element_text(face="bold", size=9),
    legend.text   = element_text(size=8),
    panel.grid.major.y = element_line(color="grey92")
  ) +
  guides(fill=guide_legend(nrow=2))

# ── PLOT 2: Genome map ────────────────────────────────────────────────────────
# Build chromosome backbone
chr_df <- data.frame(
  Chr    = names(chr_lengths),
  Length = as.numeric(chr_lengths)
) %>%
  filter(Chr %in% chr_order) %>%
  mutate(Chr = factor(Chr, levels=chr_order))

p2 <- ggplot() +
  # Chromosome backbones
  geom_segment(data=chr_df,
               aes(x=0, xend=Length/1e6, y=Chr, yend=Chr),
               color="grey85", linewidth=2.5) +
  # CxSV loci colored by sharing category
  geom_point(data=df %>% arrange(Category),
             aes(x=Mid/1e6, y=Chr, color=Category, size=Category),
             alpha=0.7, shape="|") +
  scale_color_manual(values=cat_colors) +
  scale_size_manual(values=c(
    "Fixed (100 samples)"   = 3.5,
    "Very common (80-99)"   = 2.5,
    "Common (50-79)"        = 2,
    "Low frequency (20-49)" = 1.5,
    "Rare (<20)"            = 2
  )) +
  scale_y_discrete(limits=rev(chr_order)) +
  scale_x_continuous(labels=function(x) paste0(x," Mb")) +
  labs(
    title  = "B   Genomic Distribution of CxSV Loci by Population Frequency",
    x      = "Genomic position",
    y      = NULL,
    color  = "Sharing category",
    size   = "Sharing category"
  ) +
  theme_classic(base_size=11) +
  theme(
    plot.title      = element_text(face="bold", size=13),
    legend.position = "bottom",
    legend.title    = element_text(face="bold", size=9),
    legend.text     = element_text(size=8),
    axis.text.y     = element_text(size=8),
    panel.grid.major.x = element_line(color="grey94")
  ) +
  guides(
    color = guide_legend(nrow=2, override.aes=list(size=3, shape=15)),
    size  = "none"
  )

# ── PLOT 3: Per-chromosome stacked bar ────────────────────────────────────────
chr_cat <- df %>%
  count(Chr, Category) %>%
  group_by(Chr) %>%
  mutate(Pct = n / sum(n) * 100) %>%
  ungroup()

p3 <- ggplot(chr_cat, aes(x=Chr, y=Pct, fill=Category)) +
  geom_col(width=0.75, color="white", linewidth=0.2) +
  geom_hline(yintercept=50, linetype="dashed", color="grey50", linewidth=0.4) +
  scale_fill_manual(values=cat_colors) +
  scale_y_continuous(labels=function(x) paste0(x,"%"),
                     expand=expansion(mult=c(0,0.02))) +
  labs(
    title  = "C   Sharing Category Composition per Chromosome",
    x      = NULL,
    y      = "% of CxSV loci",
    fill   = "Sharing category"
  ) +
  theme_classic(base_size=11) +
  theme(
    plot.title      = element_text(face="bold", size=13),
    axis.text.x     = element_text(angle=45, hjust=1, size=8),
    legend.position = "bottom",
    legend.title    = element_text(face="bold", size=9),
    legend.text     = element_text(size=8),
    panel.grid.major.y = element_line(color="grey92")
  ) +
  guides(fill=guide_legend(nrow=2))

# ── Combine and save ──────────────────────────────────────────────────────────
combined <- (p1 / p2 / p3) +
  plot_layout(heights=c(1, 1.4, 0.9)) +
  plot_annotation(
    title    = "CxSV Population Sharing Analysis",
    subtitle = sprintf(
      "1000G-ONT cohort | n=100 samples | n=%s CxSV loci | CHM13v2.0 (hs1)\n%s%% of loci shared by >=50 samples | %s%% fixed in all 100 samples",
      formatC(n_total, format="d", big.mark=","),
      round(100*n_common/n_total,1),
      round(100*n_fixed/n_total,1)
    ),
    theme=theme(
      plot.title    = element_text(face="bold", size=15),
      plot.subtitle = element_text(color="grey35", size=9, lineheight=1.4)
    )
  )

ggsave(outfile, plot=combined, width=13, height=18, device="pdf")
cat("Saved:", outfile, "\n")

# ── Print summary table ───────────────────────────────────────────────────────
cat("\n=== Population Sharing Summary ===\n")
df %>%
  count(Category) %>%
  mutate(Pct=round(100*n/sum(n),1)) %>%
  arrange(desc(n)) %>%
  { cat(sprintf("  %-30s %6s  %5s%%\n", "Category", "Loci", "")); . } %>%
  apply(1, function(r) cat(sprintf("  %-30s %6s  %5s%%\n", r[1], r[2], r[3])))

cat(sprintf("\n  Overlap method: bedtools intersect (any coordinate overlap)\n"))
cat(sprintf("  This means calls within the same CxSV locus boundaries count\n"))
cat(sprintf("  as 'shared', with up to ~%d bp wiggle room from SURVIVOR merging\n",
            300))

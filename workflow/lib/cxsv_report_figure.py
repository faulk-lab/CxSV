"""Tier 1 points graphic of the clinician report (roadmap item 16l; manuscript Figure 3c/3d style).

  tier1_points_figure(rows, tiers, out_png, max_rows=15) -> dict | None

rows: cxsv_scores.tsv rows (dicts, score_cxsv.py) of the sample's Tier 1 loci.
Left panel: total points per locus (dot on a stem from 0), dashed lines at the tier
thresholds, the complexity class beside the dot; dot colour = panel status.
Right panel: points per feature category stacked per locus (positive to the right,
negative to the left), black diamond at the total, dashed line at the Tier 1 threshold.
Category colours are those of manuscript Figure 3d. Rows are ordered by total points;
more than max_rows Tier 1 loci are capped to the max_rows highest (the caller states
the cap). With no Tier 1 locus nothing is drawn and None is returned (the report
prints a note instead of an empty plot; no_tier1_note() writes that note as the image).

matplotlib only (Agg). Font-size ladder 8/7/6 pt, open frame, frameless legends,
Type-42 fonts, as the lab's figure rules (figure-style) prescribe.
"""
from __future__ import annotations

import re

# feature -> category (manuscript fig3_scores.py FGROUP)
FGROUP = {"panel_privacy": "privacy", "catalogue_frequency": "privacy",
          "gene_dosage": "gene/region dosage", "region_dosage": "gene/region dosage",
          "exons": "exons, gene count", "gene_count": "exons, gene count", "no_genes": "exons, gene count",
          "complexity_class": "class, size", "size": "class, size",
          "caller_agreement": "evidence, artefact", "read_support": "evidence, artefact",
          "artefact_context": "evidence, artefact",
          "ct_interleaving": "chromothripsis criteria", "ct_breakpoint_tier": "chromothripsis criteria",
          "ct_cn_oscillation": "chromothripsis criteria", "ct_clustering": "chromothripsis criteria",
          "ct_bfb_modifier": "chromothripsis criteria",
          "hpo_match": "phenotype, inheritance", "omim": "phenotype, inheritance",
          "inheritance": "phenotype, inheritance"}
GORDER = ["privacy", "gene/region dosage", "exons, gene count", "class, size", "evidence, artefact",
          "chromothripsis criteria", "phenotype, inheritance"]
GCOL = {"privacy": "#8c8c8c", "gene/region dosage": "#2f5d8a", "exons, gene count": "#7fa7cf",
        "class, size": "#c46b1c", "evidence, artefact": "#b5a36b", "chromothripsis criteria": "#6a3d9a",
        "phenotype, inheritance": "#d4a5c4"}
STATUS_COL = {"private": "#2f5d8a", "shared_junctions": "#7fa7cf", "public": "#bdbdbd"}
STATUS_LABEL = {"private": "private to this genome", "shared_junctions": "partly matched in panel",
                "public": "public (in healthy panel)"}

RC = {"font.family": "sans-serif", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
      "legend.fontsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6,
      "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3,
      "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
      "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 300,
      "axes.titlelocation": "left", "axes.titleweight": "normal", "pdf.fonttype": 42, "ps.fonttype": 42}


def category_points(row):
    """{category: points} of one scores row (sum of the features of each category)."""
    out = {g: 0.0 for g in GORDER}
    for f, g in FGROUP.items():
        v = row.get(f + "_points")
        if v not in (None, "", "NA"):
            out[g] += float(v)
    return out


def _gene_label(row):
    """Most dosage-relevant affected gene (ClinGen HI/TS 1-3 first), '+k' for the rest."""
    genes = [g.split(":", 1)[0] for g in str(row.get("Genes_Affected", "")).split(",") if g and g != "none"]
    if not genes:
        return "no gene"
    curated = re.findall(r"([^;\s(]+):\w+\(HI ([^,]+),TS ([^,]+),", str(row.get("gene_dosage_value", "")))
    best = [g for g, hi, ts in curated if hi in ("1", "2", "3") or ts in ("1", "2", "3")]
    lead = best[0] if best else genes[0]
    name = rf"$\mathit{{{lead}}}$" if re.fullmatch(r"[A-Za-z0-9]+", lead) else lead
    return name + (f" +{len(genes) - 1}" if len(genes) > 1 else "")


def locus_label(row):
    s, e = int(row["Start"]), int(row["End"])
    span = f"{s / 1e6:.2f}\u2013{e / 1e6:.2f} Mb" if e - s >= 10_000 else f"{s:,}\u2013{e:,}"
    return f"{row['Chr']}:{span}  {_gene_label(row)}"


def select_rows(rows, max_rows):
    """Tier 1 rows by total points (highest first), capped at max_rows: (shown, n_total)."""
    rows = sorted(rows, key=lambda r: (-float(r["Total_Score"]), r["Locus_ID"]))
    return rows[:max_rows], len(rows)


def tier1_points_figure(rows, tiers, out_png, max_rows=15, out_pdf=None):
    if not rows:
        return None
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    shown, n_total = select_rows(rows, max_rows)
    n = len(shown)
    thr = sorted(((float(t["min"]), t["name"]) for t in tiers if float(t["min"]) >= 0), reverse=True)
    t1 = max(m for m, _ in thr) if thr else 0.99
    with plt.rc_context(RC):
        h = 1.55 + 0.27 * n
        fig, (axl, axr) = plt.subplots(1, 2, figsize=(7.0, h), sharey=True,
                                       gridspec_kw={"width_ratios": [1, 1.15], "wspace": 0.08})
        y = list(range(n))[::-1]
        totals = [float(r["Total_Score"]) for r in shown]
        cats = [category_points(r) for r in shown]
        # left: total points, stem from 0, tier thresholds
        xmax = max(totals + [t1]) + 0.25
        statuses = []
        for yi, r, tot in zip(y, shown, totals):
            st = r.get("Panel_Status", "private")
            statuses.append(st)
            axl.plot([0, tot], [yi, yi], color="#bbbbbb", lw=0.8, zorder=1)
            axl.scatter(tot, yi, s=18, color=STATUS_COL.get(st, "#2f5d8a"), zorder=3,
                        edgecolor="black" if st != "private" else "none", linewidth=0.4)
            axl.text(tot + 0.05, yi, r.get("Cmplx_Class", ""), va="center", fontsize=6)
        for m, name in thr:
            axl.axvline(m, color="black", lw=0.6, ls="--", zorder=0)
            lab = "Tier 1" if m == t1 else re.sub(r"^Tier(\d+).*", r"T\1", name)
            axl.text(m + 0.02, n - 0.35, lab, fontsize=6, ha="left", va="bottom")
        longest = max(len(r.get("Cmplx_Class", "")) for r in shown)
        axl.set_xlim(-0.08, xmax + 0.045 * longest)
        axl.set_ylim(-0.7, n - 0.1)
        axl.set_yticks(y)
        axl.set_yticklabels([locus_label(r) for r in shown])
        axl.set_xlabel("total points")
        axl.set_title("Total points per Tier 1 locus")
        present = [s for s in ("private", "shared_junctions", "public") if s in statuses]
        status_handles = [Line2D([], [], marker="o", ls="", color=STATUS_COL[st], markersize=4,
                                 markeredgecolor="black" if st != "private" else STATUS_COL[st],
                                 markeredgewidth=0.4, label=STATUS_LABEL[st]) for st in present]
        # right: stacked feature-category points, diamond = total
        used = [g for g in GORDER if any(abs(c[g]) > 1e-9 for c in cats)]
        lo = hi = 0.0
        for yi, c, tot in zip(y, cats, totals):
            pos = neg = 0.0
            for g in GORDER:
                v = c[g]
                if abs(v) < 1e-9:
                    continue
                if v > 0:
                    axr.barh(yi, v, left=pos, color=GCOL[g], height=0.7, lw=0)
                    pos += v
                else:
                    axr.barh(yi, v, left=neg, color=GCOL[g], height=0.7, lw=0)
                    neg += v
            axr.scatter(tot, yi, marker="D", color="black", s=12, zorder=4)
            lo, hi = min(lo, neg), max(hi, pos, tot)
        axr.axvline(0, color="black", lw=0.6)
        axr.axvline(t1, color="black", lw=0.6, ls="--")
        axr.set_xlim(lo - 0.12, hi + 0.15)
        axr.set_xlabel("points (diamond = total)")
        axr.set_title("Points by feature category")
        axr.tick_params(axis="y", length=0)
        handles = status_handles + [Patch(color=GCOL[g], label=g) for g in used] + \
            [Line2D([], [], marker="D", ls="", color="black", markersize=4, label="total points")]
        # one legend under both panels: dot colour (left) and categories (right)
        fig.subplots_adjust(bottom=0.62 / h, top=1 - 0.30 / h)     # fixed margins in inches: x labels fit
        fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=4,
                   handlelength=1, columnspacing=1.0, handletextpad=0.4)
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        if out_pdf:
            fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
    from PIL import Image
    with Image.open(out_png) as im:
        w, hh = im.size
    return {"path": out_png, "n_total": n_total, "n_shown": n, "width_px": w, "height_px": hh}


def no_tier1_note(out_png, text="No Tier 1 locus in this sample: no Tier 1 points graphic."):
    """The figure output of a sample without Tier 1 loci: a one-line note, not an empty plot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with plt.rc_context(RC):
        fig = plt.figure(figsize=(4.0, 0.35))
        fig.text(0.0, 0.5, text, va="center", ha="left", fontsize=8)
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)
    return out_png

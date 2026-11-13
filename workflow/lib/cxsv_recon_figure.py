"""Drawing of CxSV structures: reference row(s) and derivative row(s) (positive-control style).

Shared by the positive-control event figure (tests/positive_control/plot_events.py) and the
structural-reconstruction panel of the clinician report (workflow/scripts/clinician_report.py,
reconstruction from workflow/lib/cxsv_recon.py).

Glyphs: lettered coloured segments (a segment keeps its letter and colour on every row),
grey flanking sequence drawn shortened, hatched blocks = deleted on the derivative, arrow
blocks pointing right = reference orientation and left = inverted (reverse complement), red
triangles = novel junctions. Segments are to scale in kb; inter-chromosomal structures are
drawn schematically (coloured by chromosome, "not to scale").

Report panels add, when the drawing dict carries them: junction ids over the triangles
(J1, J2 ... = the junction table of the report), hollow triangles where no junction was called
(an incomplete walk joined across a gap), the read-depth copy number and the genes under each
reference segment (gene names in italics), and an 'ambiguous' stamp.

Not fingerprinted: reporting only.
"""
from __future__ import annotations

import textwrap
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.patches import Polygon, Rectangle     # noqa: E402

SEG_COLOURS = ["#4E79A7", "#F28E2B", "#59A14F", "#B07AA1", "#76B7B2", "#EDC948",
               "#9C755F", "#FF9DA7", "#86BCB6", "#D37295", "#A0CBE8", "#8CD17D"]
CHROM_COLOURS = ["#4E79A7", "#F28E2B", "#59A14F"]
FLANK = "#BABABA"
JUNCTION = "#D62728"
AMBIGUOUS = "#B35900"
NOTE_GREY = "#555555"
AX_WIDTH_IN = 1.75   # approximate panel width; used only to decide label placement

STYLE = {                                             # publication rules (figure-style ladder 7/6/5)
    "font.family": "sans-serif", "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7,
    "legend.fontsize": 6, "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.6,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 2, "xtick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.titlelocation": "left",
}


def apply_style():
    """The figure-style rcParams used by report panels (mechanics only: sizes, ticks, fonts)."""
    matplotlib.rcParams.update(STYLE)


def parse(el):
    return el[:-1], el[-1]


def is_flank(name):
    return name in ("p", "q") or name.startswith(("p_", "q_"))


below_labels = []   # (row y, level, x0, x1) of labels drawn under narrow blocks
span_scale = [1.0, 1.0]   # [arrow-tip scale, data units per inch]


def block(ax, x, w, y, h, colour, strand=None, hatch=None, label=None, alpha=1.0):
    if strand is None or w < 1e-9:
        ax.add_patch(Rectangle((x, y), w, h, facecolor=colour if not hatch else "white",
                               edgecolor=colour, hatch=hatch, linewidth=0.8, alpha=alpha))
    else:
        tip = min(w * 0.35, h * 0.9 * span_scale[0])
        if strand == "+":
            pts = [(x, y), (x + w - tip, y), (x + w, y + h / 2), (x + w - tip, y + h), (x, y + h)]
        else:
            pts = [(x + w, y), (x + tip, y), (x, y + h / 2), (x + tip, y + h), (x + w, y + h)]
        ax.add_patch(Polygon(pts, closed=True, facecolor=colour, edgecolor="white", linewidth=0.6))
    if label and w > 0:
        need = (0.05 * len(label) + 0.03) * span_scale[1]   # data units the text needs
        if w >= need:
            kw = dict(color="white")
            if hatch:
                kw = dict(color=colour, bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))
            ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=6,
                    fontweight="bold", **kw)
        else:
            # narrow block: label underneath, staggered when it would collide with a neighbour
            half = (0.04 * len(label) + 0.02) * span_scale[1] / 2
            cx, level = x + w / 2, 0
            while any(abs(y - py) < 1e-6 and lv == level and cx - half < px1 and cx + half > px0
                      for py, lv, px0, px1 in below_labels):
                level += 1
            below_labels.append((y, level, cx - half, cx + half))
            ax.text(cx, y - 0.04 - 0.3 * level, label, ha="center", va="top", fontsize=5,
                    fontweight="bold", color=colour)


def draw_event(ax, ev):
    """One structure: reference row(s) then derivative row(s). ev = events.yaml shape
    {id, name, gene, segments: {name: [chrom, start, end]}, paths: {der: ['p+', ...]}} plus the
    optional report keys junction_labels, gap_marks, segment_notes, status, title, subtitle."""
    segs, paths = ev["segments"], ev["paths"]
    inter = any(k.startswith("p_") for k in segs)
    internal = sorted((n for n in segs if not is_flank(n)), key=lambda n: (segs[n][0], segs[n][1]))
    used = Counter(parse(e)[0] for p in paths.values() for e in p)
    notes = ev.get("segment_notes") or {}
    jlabels = ev.get("junction_labels") or {}
    gaps = ev.get("gap_marks") or {}
    if inter:
        chroms = list(dict.fromkeys(segs[n][0] for n in segs))
        ccol = {c: CHROM_COLOURS[i % 3] for i, c in enumerate(chroms)}
        width = {}
        for n in segs:
            width[n] = 30.0 if is_flank(n) else 6.0
        colour = {n: ccol[segs[n][0]] for n in segs}
        label = {n: (segs[n][0].replace("chr", "") if is_flank(n) else (n if notes else "")) for n in segs}
    else:
        span = (segs[internal[-1]][2] - segs[internal[0]][1]) / 1000 if internal else 1.0
        width = {n: (segs[n][2] - segs[n][1]) / 1000 for n in internal}
        fl = max(span * 0.12, 1.5)
        width["p"] = width["q"] = fl
        colour = {n: SEG_COLOURS[i % len(SEG_COLOURS)] for i, n in enumerate(internal)}
        colour["p"] = colour["q"] = FLANK
        label = {n: (n[1:] if len(internal) > 9 and n[0] == "S" and n[1:].isdigit() else n)
                 for n in internal}
        label["p"] = label["q"] = ""
    h, gap = 0.7, 0.55
    jlab_x = []                                       # (row, x) of junction labels drawn (stagger)
    note_gap = 0.55 if notes and not inter else 0.0
    below_labels.clear()
    rows = []
    if inter:
        for c in chroms:
            names = [n for n in segs if segs[n][0] == c]
            names.sort(key=lambda n: segs[n][1])
            rows.append((f"{c}", names, None, None))
    else:
        rows.append(("reference", ["p"] + internal + ["q"], None, None))
    for der, path in paths.items():
        rows.append((der if der != "der" else "derivative", [parse(e)[0] for e in path],
                     [parse(e)[1] for e in path], der))
    n_ref = len(rows) - len(paths)
    total_w = max(sum(width[n] for n in names) for _, names, _, _ in rows)
    span_scale[0] = total_w / 14
    span_scale[1] = total_w * 1.04 / ev.get("ax_width_in", AX_WIDTH_IN)
    ylabels = []
    for i, (rname, names, strands, der) in enumerate(rows):
        y = -(i * (h + gap)) - (0.35 + note_gap if i >= n_ref else 0)
        x = 0.0
        ref_row = strands is None
        for j, n in enumerate(names):
            w = width[n]
            deleted = ref_row and not is_flank(n) and used.get(n, 0) == 0
            if ref_row:
                block(ax, x, w, y, h, colour[n], hatch="////" if deleted else None, label=label[n])
                if notes and not inter and n in notes:
                    _segment_note(ax, x, w, y, notes[n], colour[n])
            else:
                block(ax, x, w, y, h, colour[n], strand=strands[j], label=label[n])
                if j > 0:
                    a, sa = names[j - 1], strands[j - 1]
                    A, B = segs[a], segs[n]
                    adjacent = (sa == strands[j] == "+" and A[0] == B[0] and A[2] == B[1]) or \
                               (sa == strands[j] == "-" and A[0] == B[0] and A[1] == B[2])
                    is_gap = j in (gaps.get(der) or ())
                    if not adjacent or is_gap:
                        ax.plot([x], [y + h + 0.12], marker="v", markersize=4.5, clip_on=False, zorder=5,
                                color=JUNCTION, markerfacecolor="white" if is_gap else JUNCTION)
                        lab = (jlabels.get(der) or [None] * len(names))[j] if der in jlabels else None
                        if lab:
                            need = (0.045 * len(lab) + 0.02) * span_scale[1]
                            lvl = 0
                            while any(r == i and lv == lvl and abs(px - x) < need for r, px, lv in jlab_x):
                                lvl += 1
                            jlab_x.append((i, x, lvl))
                            ax.text(x, y + h + 0.3 + 0.22 * lvl, lab, ha="center", va="bottom", fontsize=5,
                                    color=JUNCTION)
            x += w
        ylabels.append((y + h / 2, rname))
    ax.set_yticks([yy for yy, _ in ylabels])
    ax.set_yticklabels([t for _, t in ylabels], fontsize=6)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-total_w * 0.02, total_w * 1.02)
    ax.set_ylim(-(len(rows) - 1) * (h + gap) - 1.0 - note_gap, h + (0.75 if jlabels else 0.45))
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    if inter:
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
        ax.text(1.0, -0.02, "not to scale", transform=ax.transAxes, ha="right", va="top",
                fontsize=6, color="#555555")
    else:
        ax.tick_params(axis="x", labelsize=6, length=2)
        ax.set_xlabel("kb (flanks shortened)", fontsize=6, labelpad=1)
    if "subtitle" in ev and ev.get("subtitle") is not None:
        ax.annotate(ev["subtitle"], (0, 1), xycoords="axes fraction", xytext=(0, 3),
                    textcoords="offset points", fontsize=6, color="#444444", va="bottom", ha="left",
                    fontstyle="italic")
    else:
        gene = ev.get("gene")
        loc = segs[internal[0]][0] if not inter else "/".join(c for c in chroms)
        if gene:
            gtxt = "/".join("$\\it{" + g + "}$" for g in gene.split("/"))
        else:
            gtxt = "intergenic"
        spike = ev.get("spike_in")
        if spike:      # add-on event planted in background genomes too (E13): expected public
            gtxt += f"  ({loc}; also in {len(spike['backgrounds'])} panel genomes)"
        else:
            gtxt += "  (" + loc + ")"
        ax.annotate(gtxt, (0, 1), xycoords="axes fraction", xytext=(0, 3),
                    textcoords="offset points", fontsize=6, color="#444444", va="bottom", ha="left")
    title = ev.get("title") or textwrap.fill(f"{ev['id']}  {ev['name']}", 34)
    if ev.get("title_line2"):
        title = title + "\n" + ev["title_line2"]
    ax.annotate(title, (0, 1), xycoords="axes fraction", xytext=(0, 12),
                textcoords="offset points", fontsize=7, va="bottom", ha="left")
    if ev.get("status") == "ambiguous":
        ax.annotate("ambiguous (see text)", (1, 1), xycoords="axes fraction", xytext=(0, 12),
                    textcoords="offset points", fontsize=6, va="bottom", ha="right", color=AMBIGUOUS,
                    fontweight="bold")


def _chrom_key(c):
    s = c.replace("chr", "")
    return (0, int(s), "") if s.isdigit() else (1, 0, s)


def _segment_note(ax, x, w, y, note, colour):
    """Copy number and genes under a reference segment, only when the block is wide enough."""
    lines = [t for t in (note.get("cn"), note.get("genes")) if t]
    if not lines:
        return
    need = max(0.04 * len(t) + 0.02 for t in lines) * span_scale[1]
    if w < need:
        return
    yy = y - 0.14
    if note.get("cn"):
        ax.text(x + w / 2, yy, note["cn"], ha="center", va="top", fontsize=5, color=NOTE_GREY,
                bbox=dict(boxstyle="square,pad=0.05", fc="white", ec="none"))
        yy -= 0.26
    if note.get("genes"):
        ax.text(x + w / 2, yy, note["genes"], ha="center", va="top", fontsize=5, color="#333333",
                fontstyle="italic")


LEGEND_ITEMS = [("hatch", "deleted on the derivative"), ("+", "reference orientation"),
                ("-", "inverted (reverse complement)"), ("flank", "flanking sequence"),
                ("junction", "novel junction")]


def draw_legend(fig, rect=(0.02, 0.0, 0.96, 0.03), items=None):
    """The shared one-line legend under the panels."""
    items = items or LEGEND_ITEMS
    lax = fig.add_axes(rect)
    lax.set_xlim(0, 100); lax.set_ylim(0, 1); lax.axis("off")
    r = fig.canvas.get_renderer()
    units_per_px = 100 / lax.get_window_extent(r).width
    widths = [lax.text(0, 0, t, fontsize=6).get_window_extent(r).width * units_per_px for _, t in items]
    for t in list(lax.texts):
        t.remove()
    gap = (100 - sum(widths) - 4.6 * len(items)) / max(1, len(items) - 1)
    x = 0.0
    for (kind, text), tw in zip(items, widths):
        if kind == "hatch":
            lax.add_patch(Rectangle((x, 0.2), 3, 0.6, facecolor="white", edgecolor="#555555", hatch="////"))
        elif kind == "flank":
            lax.add_patch(Rectangle((x, 0.2), 3, 0.6, facecolor=FLANK))
        elif kind == "+":
            lax.add_patch(Polygon([(x, 0.2), (x + 2, 0.2), (x + 3, 0.5), (x + 2, 0.8), (x, 0.8)], facecolor="#666666"))
        elif kind == "-":
            lax.add_patch(Polygon([(x + 3, 0.2), (x + 1, 0.2), (x, 0.5), (x + 1, 0.8), (x + 3, 0.8)], facecolor="#666666"))
        elif kind == "gap":
            lax.plot([x + 1.5], [0.5], marker="v", color=JUNCTION, markerfacecolor="white", markersize=5)
        else:
            lax.plot([x + 1.5], [0.5], marker="v", color=JUNCTION, markersize=5)
        lax.text(x + 4.0, 0.5, text, fontsize=6, va="center")
        x += 4.6 + tw + gap
    return lax


# ─── report panel ────────────────────────────────────────────────────────────
def recon_figure(ev, out_png, width_in=6.6):
    """One reconstruction panel (+ legend) as PNG. Returns {path, width_px, height_px}."""
    with matplotlib.rc_context(STYLE):
        n_rows = len(ev["paths"]) + (len({v[0] for v in ev["segments"].values()})
                                     if any(k.startswith("p_") for k in ev["segments"]) else 1)
        height = 0.95 + 0.42 * n_rows + (0.25 if ev.get("segment_notes") else 0)
        fig = plt.figure(figsize=(width_in, height))
        ax = fig.add_axes((0.11, 0.30, 0.86, 0.50))
        ev = dict(ev, ax_width_in=width_in * 0.86)
        draw_event(ax, ev)
        items = list(LEGEND_ITEMS)
        if any(ev.get("gap_marks", {}).values()):
            items.append(("gap", "no junction called (gap)"))
        draw_legend(fig, rect=(0.02, 0.01, 0.96, 0.07), items=items)
        fig.savefig(out_png, dpi=300)
        w, h = fig.get_size_inches() * 300
        plt.close(fig)
    return {"path": out_png, "width_px": int(round(w)), "height_px": int(round(h))}


def truth_vs_recon_figure(pairs, out_png):
    """Validation figure: per positive-control Tier 1 locus, truth (left) and reconstruction (right)."""
    import cxsv_recon as R
    pairs = sorted(pairs, key=lambda t: (t[2]["id"], t[0]["Locus_ID"]))
    with matplotlib.rc_context(STYLE):
        n = len(pairs)
        fig, axes = plt.subplots(n, 2, figsize=(7.2, 1.5 * n + 0.7), squeeze=False)
        for k, (row, res, truth, cmp) in enumerate(pairs):
            t = dict(truth, ax_width_in=2.9)
            t["title"] = f"{truth['id']} truth: {textwrap.shorten(truth.get('name', ''), 34, placeholder='...')}"
            draw_event(axes[k][0], t)
            ev = R.layout(res)
            ev.update(segment_notes=None, junction_labels=None, status=None, title_line2=None, ax_width_in=2.9)
            amb = f"; ambiguous, {cmp['n_equivalent']} equal orders" if cmp["n_equivalent"] > 1 else (
                "; ambiguous" if cmp["status"] != "resolved" else "")
            ev["title"] = f"{row['Locus_ID']}: {cmp['outcome'].replace('_', ' ')}{amb}"
            ev["subtitle"] = "reconstructed from the run outputs"
            draw_event(axes[k][1], ev)
        fig.tight_layout(rect=(0, 0.035, 1, 0.995), h_pad=2.4, w_pad=2.0)
        draw_legend(fig, rect=(0.05, 0.016, 0.9, 0.012), items=LEGEND_ITEMS[:3])
        draw_legend(fig, rect=(0.05, 0.002, 0.9, 0.012),
                    items=LEGEND_ITEMS[3:] + [("gap", "no junction called (gap)")])
        fig.savefig(out_png, dpi=300)
        fig.savefig(out_png.rsplit(".", 1)[0] + ".pdf")
        plt.close(fig)
    return out_png

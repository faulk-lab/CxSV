#!/usr/bin/env python3
"""Draw the structure of every planted positive-control CxSV from events.yaml.

One panel per event. Top track: the reference, segments to scale (kb), deleted
segments hatched. Lower track(s): each derivative haplotype, the same segments
in their rearranged order; an arrow pointing right is in reference orientation,
left is reverse-complemented. A segment keeps its letter and colour in every
track. Triangles mark the novel junctions listed in the truth set. Flanks are
drawn as short stubs; inter-chromosomal events (translocation, chromoplexy) are
drawn schematically, coloured by chromosome, not to scale. Spike-in events
(spike_in in events.yaml, e.g. E13) say how many panel genomes also carry them.

Usage: plot_events.py events.yaml out_prefix   (writes out_prefix.png and .pdf)
"""
import sys
import textwrap
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
import yaml

SEG_COLOURS = ["#4E79A7", "#F28E2B", "#59A14F", "#B07AA1", "#76B7B2", "#EDC948",
               "#9C755F", "#FF9DA7", "#86BCB6", "#D37295", "#A0CBE8", "#8CD17D"]
CHROM_COLOURS = ["#4E79A7", "#F28E2B", "#59A14F"]
FLANK = "#BABABA"
JUNCTION = "#D62728"
AX_WIDTH_IN = 1.75   # approximate panel width; used only to decide label placement


def parse(el):
    return el[:-1], el[-1]


def is_flank(name):
    return name in ("p", "q") or name.startswith(("p_", "q_"))


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


below_labels = []   # (row y, level, x0, x1) of labels drawn under narrow blocks
span_scale = [1.0, 1.0]   # [arrow-tip scale, data units per inch]


def draw_event(ax, ev):
    segs, paths = ev["segments"], ev["paths"]
    inter = any(k.startswith("p_") for k in segs)
    internal = sorted((n for n in segs if not is_flank(n)), key=lambda n: (segs[n][0], segs[n][1]))
    used = Counter(parse(e)[0] for p in paths.values() for e in p)
    if inter:
        chroms = list(dict.fromkeys(segs[n][0] for n in segs))
        ccol = {c: CHROM_COLOURS[i % 3] for i, c in enumerate(chroms)}
        width = {}
        for n in segs:
            width[n] = 30.0 if is_flank(n) else 6.0
        colour = {n: ccol[segs[n][0]] for n in segs}
        label = {n: (segs[n][0].replace("chr", "") if is_flank(n) else "") for n in segs}
    else:
        span = (segs[internal[-1]][2] - segs[internal[0]][1]) / 1000
        width = {n: (segs[n][2] - segs[n][1]) / 1000 for n in internal}
        fl = max(span * 0.12, 1.5)
        width["p"] = width["q"] = fl
        colour = {n: SEG_COLOURS[i % len(SEG_COLOURS)] for i, n in enumerate(internal)}
        colour["p"] = colour["q"] = FLANK
        label = {n: (n[1:] if len(internal) > 9 and n[0] == "S" and n[1:].isdigit() else n)
                 for n in internal}
        label["p"] = label["q"] = ""
    h, gap = 0.7, 0.55
    below_labels.clear()
    # reference track(s)
    y = 0.0
    rows = []
    if inter:
        for c in chroms:
            names = [n for n in segs if segs[n][0] == c]
            names.sort(key=lambda n: segs[n][1])
            rows.append((f"{c}", names, None))
    else:
        rows.append(("reference", ["p"] + internal + ["q"], None))
    for der, path in paths.items():
        rows.append((der if der != "der" else "derivative", [parse(e)[0] for e in path], [parse(e)[1] for e in path]))
    n_ref = len(rows) - len(paths)
    total_w = max(sum(width[n] for n in names) for _, names, _ in rows)
    span_scale[0] = total_w / 14
    span_scale[1] = total_w * 1.04 / AX_WIDTH_IN
    ylabels = []
    for i, (rname, names, strands) in enumerate(rows):
        y = -(i * (h + gap)) - (0.35 if i >= n_ref else 0)
        x = 0.0
        ref_row = strands is None
        for j, n in enumerate(names):
            w = width[n]
            deleted = ref_row and not is_flank(n) and used.get(n, 0) == 0
            if ref_row:
                block(ax, x, w, y, h, colour[n], hatch="////" if deleted else None,
                      label=label[n])
            else:
                block(ax, x, w, y, h, colour[n], strand=strands[j], label=label[n])
                if j > 0:
                    a, sa = names[j - 1], strands[j - 1]
                    A, B = segs[a], segs[n]
                    adjacent = (sa == strands[j] == "+" and A[0] == B[0] and A[2] == B[1]) or \
                               (sa == strands[j] == "-" and A[0] == B[0] and A[1] == B[2])
                    if not adjacent:
                        ax.plot([x], [y + h + 0.12], marker="v", color=JUNCTION, markersize=4.5,
                                clip_on=False, zorder=5)
            x += w
        ylabels.append((y + h / 2, rname))
    ax.set_yticks([yy for yy, _ in ylabels])
    ax.set_yticklabels([t for _, t in ylabels], fontsize=6)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-total_w * 0.02, total_w * 1.02)
    ax.set_ylim(-(len(rows) - 1) * (h + gap) - 1.0, h + 0.45)
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
    title = textwrap.fill(f"{ev['id']}  {ev['name']}", 34)
    ax.annotate(title, (0, 1), xycoords="axes fraction", xytext=(0, 12),
                textcoords="offset points", fontsize=7, va="bottom", ha="left")


def main(events_yaml, out):
    doc = yaml.safe_load(open(events_yaml))
    evs = doc["events"]
    ncol = 3
    nrow = (len(evs) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, 1.8 * nrow + 0.5))
    for ax, ev in zip(axes.flat, evs):
        draw_event(ax, ev)
    for ax in axes.flat[len(evs):]:
        ax.axis("off")
    fig.tight_layout(rect=(0, 0.04, 1, 0.985), h_pad=2.6, w_pad=1.2)
    lax = fig.add_axes((0.02, 0.0, 0.96, 0.03))
    lax.set_xlim(0, 100); lax.set_ylim(0, 1); lax.axis("off")
    items = [("hatch", "deleted on the derivative"), ("+", "reference orientation"),
             ("-", "inverted (reverse complement)"), ("flank", "flanking sequence"),
             ("junction", "novel junction")]
    r = fig.canvas.get_renderer()
    units_per_px = 100 / lax.get_window_extent(r).width
    widths = [lax.text(0, 0, t, fontsize=6).get_window_extent(r).width * units_per_px for _, t in items]
    for t in list(lax.texts):
        t.remove()
    gap = (100 - sum(widths) - 4.6 * len(items)) / (len(items) - 1)
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
        else:
            lax.plot([x + 1.5], [0.5], marker="v", color=JUNCTION, markersize=5)
        lax.text(x + 4.0, 0.5, text, fontsize=6, va="center")
        x += 4.6 + tw + gap
    fig.savefig(out + ".png", dpi=300)
    fig.savefig(out + ".pdf")
    return fig


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])

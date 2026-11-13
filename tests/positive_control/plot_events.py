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

The drawing primitives live in workflow/lib/cxsv_recon_figure.py, shared with the
structural-reconstruction panel of the clinician report.

Usage: plot_events.py events.yaml out_prefix   (writes out_prefix.png and .pdf)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import yaml  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "workflow", "lib"))
from cxsv_recon_figure import draw_event, draw_legend  # noqa: E402,F401


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
    draw_legend(fig, rect=(0.02, 0.0, 0.96, 0.03))
    fig.savefig(out + ".png", dpi=300)
    fig.savefig(out + ".pdf")
    return fig


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])

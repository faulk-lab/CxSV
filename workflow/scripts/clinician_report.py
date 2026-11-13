#!/usr/bin/env python3
"""Clinician-facing CxSV research report of one patient sample (.docx + .pdf; roadmap item 16l).

Ported from the accepted test-run generator (make_clinician_reports.py, 2026-11-11) into
the patient workflow (rule clinician_report, workflow/rules/clinician_report.smk). Every
number comes from this run's outputs; the fixed prose (tier meanings, feature table,
limitations, footer) from config/clinician_report.yaml. Structure (2 pages, PI-accepted):

  header        sample, test, reference, bundle, panel release, calling fingerprint + status
  result box    one line: Tier 1 count, or "No Tier 1 ..."; caution when QC flags apply
  1. Sample quality   coverage, read length, library QC, sex check (sex_call, x_ratio,
                      y_ratio, ambiguous flag), depth noise (CN sigma vs threshold)
  2. Findings   every Tier 1 locus in full (genes, ClinGen HI/TS, panel frequency,
                score with every non-zero feature), the Tier 1 points graphic,
                Tier 2 table, Tier 3 / likely-benign counts
  3. What the tiers mean   tier box (thresholds from scoring_params.json) + main features
  4. Limitations
  footer        research use only

No per-cohort or healthy-panel comparator text: a report describes one sample.
Reports never state a diagnosis.

Length: 2 pages without Tier 1 reconstructions; each reconstruction panel may add pages, up
to page_limit() = 2 + ceil(panels / reconstruction_panels_per_page) (config; default 2 per
page). The PDF page count and the limit are printed in the run log; over the limit the
script warns (the report is still written).

  clinician_report.py --sample S --scores scoring/cxsv_scores.tsv --params scoring/scoring_params.json \
      --cn-stats cn/stats.json --cn-only cn/cn_only.tsv --library-qc qc/library_qc.tsv \
      --depth-summary cn/depth.mosdepth.summary.txt --compatibility panel_compatibility.txt \
      --text config/clinician_report.yaml|text.json --out-docx X.docx --out-pdf X.pdf --out-figure X.png \
      [--junctions junctions/junctions.tsv --cxsv-vcf cxsv/cxsv_only.vcf --caller-vcf Sniffles2=... \
       --caller-vcf cuteSV=... --caller-vcf SVIM=... --locus-depth scoring/locus_depth.tsv \
       --cn-segments cn/segments.tsv --genes genes.bed --out-recon-dir DIR --out-recon-tsv X.tsv]

Bundle identity is read from the environment the bundle launcher sets (CXSV_BUNDLE_ID,
CXSV_COMMIT; container/run_cxsv.sh); outside a bundle the report says so.
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_report_figure as rfig  # noqa: E402
import cxsv_scoring as sc          # noqa: E402

TOP = "Tier1_high_priority"
CLASS_PLAIN = {
    "DEL+INV": "deletion with an inverted segment",
    "DUP+INV": "duplication with an inverted segment",
    "DEL+DUP+INV": "deletion, duplication and inversion in one event",
    "Chromothripsis_like": "chromothripsis-like cluster (several clustered breaks rejoined)",
    "Multi_breakpoint_cluster": "cluster of several nearby breakpoints",
    "Nested_SV": "nested structural variant (one SV inside another)",
    "Complex_inversion": "complex inversion",
    "Complex_other": "complex rearrangement, other",
    "BND_cluster": "cluster of breakend (junction-only) calls",
    "CN_plus_orientation_change": "copy-number change with an orientation change",
    "DUP-TRP_INV-DUP": "duplication-triplication with an inverted duplication (DUP-TRP/INV-DUP)",
    "CN_linked_loss": "copy-number loss linked to a junction call",
    "CN_linked_gain": "copy-number gain linked to a junction call",
    "CN_linked_gain_loss": "linked copy-number gain and loss",
    "CN_nested_gain": "nested copy-number gain",
}
HOW = {"broken": "disrupted", "lost": "deleted", "gained": "duplicated"}


# ─── inputs ──────────────────────────────────────────────────────────────────
def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("# ")), delimiter="\t"))


def read_json(path):
    if not path or not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def read_text_config(path):
    """config/clinician_report.yaml (key clinician_report) or the JSON the rule writes."""
    with open(path) as fh:
        if path.endswith(".json"):
            d = json.load(fh)
        else:
            import yaml
            d = yaml.safe_load(fh)
    return d.get("clinician_report", d)


def mean_coverage(path):
    """'total' mean depth of a mosdepth summary (cn/depth.mosdepth.summary.txt), else None."""
    for r in read_tsv(path):
        if r.get("chrom") == "total":
            return float(r["mean"])
    return None


def compatibility(path):
    """(run fingerprint, status) from panel_compatibility.txt."""
    fp, status = None, None
    if path and os.path.exists(path):
        for line in open(path):
            if line.startswith("this run fingerprint:"):
                fp = line.split(":", 1)[1].strip()
            elif line.startswith("status:"):
                status = line.split(":", 1)[1].strip()
    return fp, status


def load(a):
    T = read_text_config(a.text)
    P = read_json(a.params)
    D = {"sample": a.sample, "T": T, "P": P, "scores": read_tsv(a.scores), "cn": read_json(a.cn_stats) or None,
         "cn_only": read_tsv(a.cn_only), "libqc": (read_tsv(a.library_qc) or [{}])[0],
         "coverage": mean_coverage(a.depth_summary), "fp": compatibility(a.compatibility),
         "quarantine": read_tsv(a.quarantine) if getattr(a, "quarantine", None) else None,
         "bundle": os.environ.get("CXSV_BUNDLE_ID") or None, "commit": os.environ.get("CXSV_COMMIT") or None,
         "date": a.date or _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d")}
    return D


# ─── formatting ──────────────────────────────────────────────────────────────
def fmt_int(x):
    try:
        return f"{int(float(x)):,}"
    except (TypeError, ValueError):
        return str(x)


def size_txt(n):
    n = int(float(n))
    return f"{n / 1e6:.2f} Mb" if n >= 1e6 else (f"{n / 1e3:.1f} kb" if n >= 1e3 else f"{n} bp")


def malformed(r):
    return int(r["End"]) - int(r["Start"]) < 10


def coords(r):
    if malformed(r):
        return f"{r['Locus_ID']} (malformed locus coordinates)"
    return f"{r['Chr']}:{fmt_int(r['Start'])}\u2013{fmt_int(r['End'])}"


def dosage_genes(row):
    """Affected genes with ClinGen HI/TS 1-3 (from gene_dosage_value), in table order."""
    return [g for g, hi, ts in re.findall(r"([^;\s(]+):\w+\(HI ([^,]+),TS ([^,]+),", str(row.get("gene_dosage_value", "")))
            if hi in ("1", "2", "3") or ts in ("1", "2", "3")]


def genes_ranked(row):
    """Genes_Affected with the dosage-sensitive genes first (result box, tables)."""
    items = [x for x in str(row.get("Genes_Affected", "")).split(",") if x and x not in ("-", "none")]
    key = set(dosage_genes(row))
    return ",".join(sorted(items, key=lambda x: (x.split(":", 1)[0] not in key)))


def short_value(f, v, r):
    """Compact feature value for the report: dosage-sensitive genes only for gene_dosage,
    else the value cut at 70 characters (the scores table keeps the full value)."""
    if f == "gene_dosage":
        keep = [m.group(0) for m in re.finditer(r"([^;\s(]+):\w+\(HI ([^,]+),TS ([^,]+),LOEUF ([^)]+)\)", v)
                if m.group(1) in set(dosage_genes(r))]
        return "; ".join(keep) or v[:70]
    return v if len(v) <= 70 else v[:67].rstrip() + "\u2026"


def genes_plain(g, max_genes=None):
    """'NF1:broken,RNF135:lost' -> 'NF1 (disrupted), RNF135 (deleted)'."""
    out = []
    for item in str(g).split(","):
        if ":" in item:
            name, how = item.split(":", 1)
            out.append(f"{name} ({HOW.get(how, how)})")
        elif item and item not in ("-", "none"):
            out.append(item)
    if not out:
        return "none"
    if max_genes and len(out) > max_genes:
        k = len(out) - max_genes
        return ", ".join(out[:max_genes]) + f" and {k} more gene{'' if k == 1 else 's'}"
    return ", ".join(out)


def dosage_plain(row):
    """Genes with ClinGen HI/TS 1-3 or LOEUF < 0.35, plus curated ClinGen regions hit."""
    keep = []
    for m in re.finditer(r"([^;\s(]+):(\w+)\(HI ([^,]+),TS ([^,]+),LOEUF ([^)]+)\)", str(row.get("gene_dosage_value"))):
        g, how, hi, ts, lo = m.groups()
        try:
            low = lo != "NA" and float(lo) < 0.35
        except ValueError:
            low = False
        if hi in ("1", "2", "3") or ts in ("1", "2", "3") or low:
            keep.append(f"{g} {HOW.get(how, how)}: ClinGen HI {hi}, TS {ts}, LOEUF {lo}")
    reg = str(row.get("region_dosage_value", "none"))
    if reg not in ("none", "NA", "") and not reg.startswith("off"):
        keep.append(f"ClinGen region: {reg}")
    return "; ".join(keep) if keep else "no gene or region with ClinGen HI/TS 1\u20133"


def population(r):
    carriers = int(float(r.get("Panel_Carriers") or 0))
    n = fmt_int(r.get("Panel_Individuals"))
    st = r.get("Panel_Status")
    if st == "private":
        pan = f"private: 0 of {n} panel individuals (callable in {fmt_int(r.get('Panel_Callable'))})"
    elif st == "public":
        pan = f"public: {fmt_int(carriers)} of {n} presumed-healthy panel individuals carry it"
    else:
        pan = f"partly matched: some junctions seen in the panel, no individual carries the full event"
    afs = [float(x) for x in re.findall(r"AF ([0-9.eE+-]+)", str(r.get("catalogue_frequency_value", "")))]
    if str(r.get("Catalogues_Containing")) not in ("none", "NA", "") and \
            str(r.get("Catalogue_Max_AF")) not in (".", "NA", ""):
        cat = (f"external catalogues: in {str(r['Catalogues_Containing']).replace(',', ', ')}, "
               f"max AF {float(r['Catalogue_Max_AF']):.2g}")
    elif afs:
        cat = f"external catalogues: partial junction matches only (max AF {max(afs):.2g})"
    else:
        cat = "external catalogues: not found"
    return pan, cat


def short_population(r):
    st, n = r.get("Panel_Status"), fmt_int(r.get("Panel_Individuals"))
    if st == "private":
        return f"private (0/{n})"
    if st == "public":
        return f"public ({fmt_int(r.get('Panel_Carriers'))}/{n})"
    return "partly matched (0 carriers)"


def all_points(r):
    """Every non-zero feature of a scores row: [(feature, points, value)] by |points|."""
    out = []
    for f in sc.FEATURE_ORDER:
        p = float(r.get(f + "_points") or 0)
        if abs(p) > 1e-9:
            out.append((f, p, str(r.get(f + "_value", ""))))
    out.sort(key=lambda x: (-abs(x[1]), sc.FEATURE_ORDER.index(x[0])))
    return out


def gate_notes(r):
    """Why the reported tier differs from the tier by points (ceiling / library QC / depth noise)."""
    out = []
    for col, label in (("Tier_Ceiling", None), ("Library_QC_Gate", "library QC conservative mode"),
                       ("Depth_Noise_Gate", "depth-noise gate")):
        v = str(r.get(col, "not_applied"))
        if v in ("not_applied", "NA", ""):
            continue
        out.append(v if label is None else f"{label}: {v}")
    return out


def interpretation(r, D):
    """Automatic, rule-based one-line interpretation (no curation, never a diagnosis)."""
    tier, st = r.get("Tier"), r.get("Panel_Status")
    noisy = sample_noise(D)["state"] == "high"
    if malformed(r):
        return "Locus coordinates are malformed (pipeline defect); not interpretable."
    if r.get("CN_Derived") == "yes" and noisy:
        return "Copy-number-derived call in a sample with a noisy depth profile; unreliable."
    if str(r.get("Tier_Ceiling", "not_applied")) not in ("not_applied", "NA", ""):
        return (f"Reaches Tier 1 by points but is carried by {fmt_int(r.get('Panel_Carriers'))} presumed-healthy "
                "panel individuals: capped at Tier 2; most likely a benign polymorphism.")
    if tier == TOP:
        return ("High-priority candidate for expert review. It must be confirmed in a clinically validated "
                "(accredited) assay before any clinical use; this report does not make a diagnosis.")
    if st == "public":
        return f"Recurrent public locus ({fmt_int(r.get('Panel_Carriers'))} panel carriers); likely benign polymorphism."
    if st == "private":
        return "Uncertain significance; review only with phenotype correlation."
    return "Junctions partly seen in healthy genomes; uncertain, most likely low significance."


# ─── sample-level QC ─────────────────────────────────────────────────────────
def sample_sex(D):
    return sc.resolve_sex(D["cn"], (D["P"].get("patient") or {}).get("sex"))


def sample_noise(D):
    cfg = D["P"].get("depth_noise")
    if cfg is None:
        return {"state": "not_configured", "sigma": (D["cn"] or {}).get("sigma"), "threshold": None,
                "stratum": (D["cn"] or {}).get("stratum", "NA"),
                "text": "depth-noise gate not configured in this run's scoring parameters"}
    return sc.depth_noise(D["cn"], cfg)


def qc_rows(D):
    q, SX, dn = D["libqc"], sample_sex(D), sample_noise(D)
    rows = [["Mean coverage", f"{D['coverage']:.1f}x (mosdepth, whole genome)" if D["coverage"] else "not available"]]
    if q.get("median_read_len"):
        rows.append(["Read length", f"median {fmt_int(q['median_read_len'])} bp "
                                    f"({fmt_int(q.get('reads_sampled'))} reads sampled by library QC)"])
    if q:
        rows.append(["Library QC", f"{q.get('qc_flag', 'NA')} (calibration group {q.get('qc_group', 'NA')}, "
                                   f"from the {q.get('qc_group_source', 'NA')}); {q.get('qc_reason', '')}".rstrip("; ")])
    if SX["ambiguous"]:
        sexv = f"AMBIGUOUS ({SX['sex_call']}) \u2014 {sc.AMBIGUOUS_SEX_FLAG}"
    elif SX["sex_call"] in ("male", "female"):
        sexv = f"{SX['sex_call']} (read-depth sex call)"
    else:
        sexv = f"not called ({SX['source']})"
    sexv += f"; sex_call {SX['sex_call']}, x_ratio {sc.fnum(SX['x_ratio'], 4)}, y_ratio {sc.fnum(SX['y_ratio'], 4)}"
    if SX["given"]:
        sexv += f"; sample sheet: {SX['given']}"
    if SX["discordant"]:
        sexv += " \u2014 DISAGREES with the read-depth call: check sample identity"
    rows.append(["Sex-chromosome check", sexv])
    n_cn_only = len(D["cn_only"])
    if dn["state"] == "high":
        dv = (f"HIGH \u2014 {dn['text']}; {n_cn_only} copy-number-only candidates. Copy-number-only and "
              "copy-number-derived calls in this sample are unreliable.")
    elif dn["state"] == "normal":
        dv = f"normal (CN sigma {sc.fnum(dn['sigma'], 3)}, threshold {sc.fnum(dn['threshold'], 3)} for {dn['stratum']})"
    else:
        dv = dn["text"]
    rows.append(["Depth-profile noise", dv])
    Q = D.get("quarantine")
    if Q is None:
        qv = "not checked (no scoring/quarantined_loci.tsv)"
    else:
        rem = [r for r in Q if r.get("action") == "removed"]
        drp = [r for r in Q if r.get("action") == "junctions_dropped"]
        contigs = sorted({c for r in Q for c in (r.get("nonref_contigs") or "").split(",") if c})
        qv = (f"{len(rem)} locus/loci removed from scoring; {sum(int(r.get('nonref_junctions') or 0) for r in drp)} "
              f"junction(s) dropped from {len(drp)} scored locus/loci")
        if contigs:
            qv += f" (contig(s) absent from the reference: {', '.join(contigs)}; caller artefact, not scored)"
    rows.append(["Quarantined (non-reference contig)", qv])
    fp, status = D["fp"]
    rows.append(["Panel compatibility", f"{status or 'NA'} (calling fingerprint {fp[:12] if fp else 'NA'})"])
    return rows


# ─── content model ───────────────────────────────────────────────────────────
def tier_threshold_text(P):
    tiers = sorted(P.get("tiers") or [], key=lambda t: -float(t["min"]))
    out = {}
    for i, t in enumerate(tiers):
        lo = float(t["min"])
        if i == 0:
            out[t["name"]] = f"\u2265 {lo:.2f}"
        elif lo <= -100:
            out[t["name"]] = f"< {float(tiers[i - 1]['min']):.2f}"
        else:
            out[t["name"]] = f"{lo:.2f} to < {float(tiers[i - 1]['min']):.2f}"
    return out


def tier1_threshold(P):
    return max(float(t["min"]) for t in P.get("tiers") or [{"min": 0.99}])


def result_box(D, t1):
    SX, dn = sample_sex(D), sample_noise(D)
    caveats = []
    if SX["ambiguous"]:
        caveats.append("sex-chromosome check ambiguous")
    if SX["discordant"]:
        caveats.append("sex disagrees with the sample sheet")
    if dn["state"] == "high":
        caveats.append("noisy depth profile: copy-number calls unreliable")
    if D["libqc"].get("qc_flag") == "FLAG":
        caveats.append("library QC flag")
    if not t1:
        text, tone = "No Tier 1 complex structural variant detected.", "neutral"
    else:
        r = t1[0]
        lead = (f"{size_txt(r['Size'])} {r['Cmplx_Class']} at {coords(r)}"
                f" ({genes_plain(genes_ranked(r), max_genes=3)})")
        more = f"; {len(t1) - 1} further Tier 1 locus/loci listed below" if len(t1) > 1 else ""
        text = (f"{len(t1)} Tier 1 finding{'s' if len(t1) > 1 else ''} requiring expert review and orthogonal "
                f"confirmation in an accredited laboratory: {lead}{more}.")
        tone = "alert"
    if caveats:
        text += " Caution: " + "; ".join(caveats) + " (see Sample quality)."
        tone = "caution" if tone == "neutral" else tone
    return text, tone


def sample_blocks(D, figure=None, recons=None):
    """The report as renderer-neutral blocks: (kind, ...)."""
    P, T, rows = D["P"], D["T"], D["scores"]
    t1 = sorted([r for r in rows if r.get("Tier") == TOP], key=lambda r: -float(r["Total_Score"]))
    t2 = sorted([r for r in rows if r.get("Tier") == "Tier2"], key=lambda r: -float(r["Total_Score"]))
    counts = {k: sum(1 for r in rows if r.get("Tier") == k) for k in ("Tier3", "Benign_likely")}
    s = D["sample"]
    B = [("title", f"Complex structural variant (CxSV) research report \u2014 {s}")]
    fp, status = D["fp"]
    n_panel = fmt_int(rows[0]["Panel_Individuals"]) if rows else "NA"
    B.append(("kv", [
        ["Sample ID", s, "Specimen", T.get("specimen", "")],
        ["Test", T.get("test", ""), "Reference genome", T.get("reference", "")],
        ["Pipeline bundle", D["bundle"] or "not run from an offline bundle",
         "Background panel", f"{P.get('panel_release_id', 'NA')}; {n_panel} presumed-healthy individuals"],
        ["Calling fingerprint", f"{fp[:12] if fp else 'NA'} ({status or 'NA'} with panel)",
         "Report date", D["date"] + (f"; code {D['commit'][:12]}" if D["commit"] else "")],
    ] + ([["Rendering", D["rendered_from"], "", ""]] if D.get("rendered_from") else [])))
    text, tone = result_box(D, t1)
    B.append(("box", "Result summary", text, tone))

    B.append(("h", "1. Sample quality"))
    B.append(("kv2", qc_rows(D)))

    B.append(("h", "2. Findings"))
    thr = tier1_threshold(P)
    if t1:
        for r in t1:
            pan, cat = population(r)
            B.append(("h2", f"Tier 1: {coords(r)} (hs1)"))
            kv = [["Size / structure", f"{size_txt(r['Size'])}; {CLASS_PLAIN.get(r['Cmplx_Class'], r['Cmplx_Class'])} "
                                       f"[{r['Cmplx_Class']}; SV types {r.get('SV_Types', 'NA')}; "
                                       f"{fmt_int(r.get('N_Junctions'))} junctions]"],
                  ["Genes affected", genes_plain(genes_ranked(r))],
                  ["Dosage sensitivity", dosage_plain(r)],
                  ["Population frequency", f"{pan}; {cat}"],
                  ["Score", f"{float(r['Total_Score']):.2f} (Tier 1 threshold {thr:.2f})"],
                  ["Contributing features", "; ".join(f"{f.replace('_', ' ')} {p:+.2f} ({short_value(f, v, r)})"
                                                      for f, p, v in all_points(r))]]
            g = gate_notes(r)
            if g:
                kv.append(["Tier gates", "; ".join(g)])
            B.append(("kv2", kv))
            B.append(("interp", interpretation(r, D)))
        if figure:
            cap = (f"Tier 1 points. Left: total points per locus (dot; colour = panel status) with the tier "
                   f"thresholds (dashed) and the complexity class. Right: points per feature category, "
                   f"diamond = total.")
            if figure["n_shown"] < figure["n_total"]:
                cap += f" Showing the {figure['n_shown']} highest-scoring of {figure['n_total']} Tier 1 loci."
            B.append(("figure", figure, cap))
        if recons:
            B += recon_blocks(recons, D.get("recon_total", len(recons)))
    else:
        B.append(("p", "Tier 1: none. (No Tier 1 locus, so no Tier 1 points graphic.)"))
    if t2:
        B.append(("h2", f"Tier 2 ({len(t2)} loci)"))
        tab = [["Locus (hs1)", "Size / class", "Genes", "Population", "Score", "Note"]]
        for r in t2:
            size = "size not reported" if malformed(r) else size_txt(r["Size"])
            tab.append([coords(r), f"{size}; {r['Cmplx_Class']}", genes_plain(genes_ranked(r), max_genes=3),
                        short_population(r), f"{float(r['Total_Score']):.2f}", interpretation(r, D)])
        B.append(("table", tab, [3.9, 3.3, 2.9, 2.0, 0.9, 4.5]))
    else:
        B.append(("p", "Tier 2: none."))
    B.append(("p", f"Tier 3: {counts['Tier3']} loci; likely benign: {counts['Benign_likely']} loci "
                   f"(of {len(rows)} CxSV loci scored). Not listed individually."
              + (f" Copy-number-only candidates (not scored): {len(D['cn_only'])}." if D["cn_only"] else "")))

    B.append(("h", "3. What the tiers mean"))
    thr_txt = tier_threshold_text(P)
    meaning = T.get("tier_meaning") or {}
    B.append(("table", [["Tier", "Total points", "Meaning"]] +
              [[sc.tier_label(k).replace("likely benign", "Likely benign"), thr_txt.get(k, "NA"), meaning.get(k, "")]
               for k in (TOP, "Tier2", "Tier3", "Benign_likely") if k in thr_txt], [2.0, 2.4, 13.1]))
    if T.get("tier_rules"):
        B.append(("small_bullets", list(T["tier_rules"])))
    B.append(("p", "A locus is first required to meet CxSV criteria (at least 3 breakpoints, nested SVs, overlapping "
                   "SV types, or a copy-number change with an inversion). Points are then summed over explicit "
                   "features modelled on ACMG/ClinGen CNV scoring; the main ones are:"))
    B.append(("table", [["Feature", "Points", "What it measures"]] + [list(f) for f in T.get("features") or []],
              [3.0, 6.5, 8.0]))
    B.append(("small", "ClinGen dosage scores: HI/TS 3 = sufficient evidence, 2 = emerging, 1 = little evidence, "
                       "0 = no evidence, 30 = autosomal recessive gene, 40 = dosage sensitivity unlikely. LOEUF = "
                       "gnomAD loss-of-function tolerance (low = intolerant). AF = allele frequency."))

    B.append(("h", "4. Limitations"))
    lim = list(T.get("limitations") or [])
    off = []
    if rows:
        r0 = rows[0]
        if str(r0.get("hpo_match_value", "")).startswith("no patient HPO"):
            off.append("phenotype (HPO) terms")
        if str(r0.get("omim_value", "")).startswith("off"):
            off.append("OMIM annotation")
        if str(r0.get("inheritance_value", "")).startswith("parents not"):
            off.append("parental samples")
    if off:
        lim.append(f"No {', '.join(off[:-1]) + (' or ' if len(off) > 1 else '') + off[-1]} were used, so those "
                   "scoring features are zero.")
    B.append(("bullets", lim))
    return B


# ─── structural reconstruction ───────────────────────────────────────────────
RECON_CAPTION = "Inferred from split-read junctions and read depth in this sample; not orthogonally confirmed."


def page_limit(n_panels, per_page=2, base_pages=2):
    """Pages a report may take: the report without reconstructions (2 pages for the accepted
    layout; more only when many Tier 1 loci are listed), plus one page per started group of
    per_page reconstruction panels."""
    return max(2, int(base_pages)) + (-(-int(n_panels) // max(1, int(per_page))) if n_panels else 0)


def recon_title(k, row, res):
    return f"Tier 1 #{k}  {coords(row)} (hs1)"


def recon_subtitle(res):
    genes = []
    for n, gl in (res.get("segment_genes") or {}).items():
        for g, _ in gl:
            if g not in genes:
                genes.append(g)
    for g in res.get("broken_genes") or ():
        if g not in genes:
            genes.append(g)
    return ", ".join(genes) if genes else "no gene in the reconstructed segments"


def recon_rows(res):
    """kv2 rows describing one reconstruction (segments, junctions with support, inference)."""
    rows = []
    paths = res.get("paths") or {}
    rows.append(["Inferred structure", (res.get("summary") or "not determined") + "; derivative "
                 + "; ".join(f"{k}: {' '.join(v)}" for k, v in paths.items()) if paths else
                 (res.get("summary") or "not determined")])
    segs = []
    for s in res.get("segments") or ():
        if s["kind"] != "internal":
            continue
        info = (res.get("segment_info") or {}).get(s["name"], {})
        visits = None
        if res.get("best"):
            idx = [x["name"] for x in res["segments"]].index(s["name"])
            visits = res["best"]["visits"][idx]
        cn = (f"CN {info['cn']:.2f} ({info['source']})" if info.get("cn") is not None
              else "CN not measurable (segment too short)")
        fate = ("deleted" if visits == 0 else (f"{visits} copies" if visits and visits > 1 else "kept")) \
            if visits is not None else "not placed"
        gl = (res.get("segment_genes") or {}).get(s["name"], [])
        gtxt = ", ".join(f"{g} ({w})" for g, w in gl) if gl else "no gene"
        segs.append(f"{s['name']} {s['chrom']}:{s['start']:,}-{s['end']:,} ({size_txt(s['end'] - s['start'])}): "
                    f"{cn}; on the derivative: {fate}; genes: {gtxt}")
    if segs:
        rows.append(["Segments", "; ".join(segs)])
    js = []
    for j in res.get("junctions") or ():
        sup = ", ".join(f"{k} {v} reads" for k, v in sorted(j["callers"].items())) or "support not recorded"
        types = "; ".join(f"{k} {'/'.join(v)}" for k, v in j["types"].items())
        js.append(f"{j['id']} {j['end1']} | {j['end2']}: join {', '.join(j['orientations'])}"
                  f"{'' if j['used'] else ' (not on the drawn derivative)'}; {sup}; types: {types}")
    if js:
        rows.append(["Junctions (support)", "; ".join(js)])
    if res.get("status") != "resolved":
        alts = res.get("alternatives") or []
        txt = "; ".join(res.get("reasons") or [])
        if alts:
            txt += ". Alternatives: " + " || ".join(" | ".join(" ".join(p) for p in a) for a in alts)
        rows.append(["Ambiguous" if res.get("status") == "ambiguous" else "Not reconstructed", txt])
    if res.get("notes"):
        rows.append(["Notes", "; ".join(res["notes"])])
    return rows


def make_reconstructions(a, D, t1):
    """Reconstruct the Tier 1 loci (highest score first, at most reconstruction_max_panels) and
    draw one panel each. Returns [{row, res, figure}] (empty without the inputs or Tier 1)."""
    if not getattr(a, "junctions", None) or not t1:
        return []
    import cxsv_recon as R
    import cxsv_recon_figure as RF
    T = D["T"]
    n_max = int(T.get("reconstruction_max_panels", 6))
    rows = sorted(t1, key=lambda r: -float(r["Total_Score"]))
    shown = rows[:n_max]
    callers = dict(x.split("=", 1) for x in (a.caller_vcf or []))
    P = D["P"]
    par_bed = (P.get("scoring_resources") or {}).get("par_bed")
    par = [(c, s, e) for c, s, e, _ in sc.read_bed(par_bed)] if par_bed and os.path.exists(par_bed) else []
    sex = sample_sex(D)["sex"]
    params = {**P, "copy_number": P.get("copy_number") or {}}
    recs = R.reconstruct_run(shown, a.junctions, callers, a.cxsv_vcf, a.locus_depth, a.cn_segments, a.genes,
                             params, sex, par)
    out = []
    if a.out_recon_dir:
        os.makedirs(a.out_recon_dir, exist_ok=True)
    for k, (row, res) in enumerate(zip(shown, recs), 1):
        fig = None
        if res.get("segments") and res.get("paths") and a.out_recon_dir:
            ev = R.layout(res, title=recon_title(k, row, res), subtitle=recon_subtitle(res))
            safe = re.sub(r"[^A-Za-z0-9]+", "_", row["Locus_ID"]).strip("_")
            fig = RF.recon_figure(ev, os.path.join(a.out_recon_dir, f"{k:02d}_{safe}.png"))
        out.append({"rank": k, "row": row, "res": res, "figure": fig})
    if a.out_recon_tsv:
        with open(a.out_recon_tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=R.TSV_COLUMNS, delimiter="\t")
            w.writeheader()
            for x in out:
                w.writerows(R.flat_rows(D["sample"], x["rank"], x["res"]))
    D["recon_total"] = len(rows)
    return out


def recon_blocks(recons, n_total):
    B = [("h2", "Structural reconstruction of the Tier 1 loci")]
    lead = ("For each Tier 1 locus the reference segments between its junction breakpoints are shown above the "
            "rearranged sequence inferred from this run's junction calls and read depth (to scale; flanks "
            "shortened; J1, J2 ... label the junctions listed under the panel).")
    if n_total > len(recons):
        lead += f" Shown: the {len(recons)} highest-scoring of {n_total} Tier 1 loci; all are in the reconstruction table."
    B.append(("small_kwn", lead))
    for x in recons:
        res = x["res"]
        cap = RECON_CAPTION
        if res.get("status") == "ambiguous":
            cap += " Ambiguous: " + "; ".join(res.get("reasons") or []) + "."
        if x["figure"]:
            B.append(("recon_figure", x["figure"], cap))
        else:
            B.append(("small", f"Tier 1 #{x['rank']} {coords(x['row'])}: no structure drawn. {cap}"))
        B.append(("kv2", recon_rows(res)))
    return B


# ─── renderers ───────────────────────────────────────────────────────────────
TONE = {"neutral": "E8F1F8", "caution": "FFF4D6", "alert": "FDE2D6"}


def strip_tags(t):
    return re.sub(r"</?b>", "", t)


def to_docx(blocks, path, footer):
    from docx import Document
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Cm, Pt, RGBColor

    def shade(cell, hexcolor):
        tcPr = cell._tc.get_or_add_tcPr()
        shd = OxmlElement("w:shd")
        shd.set(qn("w:val"), "clear")
        shd.set(qn("w:color"), "auto")
        shd.set(qn("w:fill"), hexcolor)
        tcPr.append(shd)

    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.59), Cm(27.94)        # US letter, as the PDF
    sec.left_margin = sec.right_margin = Cm(1.6)
    sec.top_margin = sec.bottom_margin = Cm(1.4)
    st = doc.styles["Normal"]
    st.font.name = "Calibri"
    st.font.size = Pt(9)
    fp = sec.footer.paragraphs[0]
    fp.text = footer
    fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    fp.runs[0].font.size = Pt(7.5)
    fp.runs[0].italic = True
    for b in blocks:
        k = b[0]
        if k == "title":
            r = doc.add_paragraph().add_run(b[1])
            r.bold, r.font.size = True, Pt(13)
        elif k == "h":
            r = doc.add_paragraph().add_run(b[1])
            r.bold, r.font.size = True, Pt(10.5)
            r.font.color.rgb = RGBColor(0x1F, 0x3A, 0x5F)
        elif k == "h2":
            r = doc.add_paragraph().add_run(b[1])
            r.bold, r.font.size = True, Pt(9.5)
        elif k in ("p", "small", "small_kwn"):
            p = doc.add_paragraph()
            for part in re.split(r"(<b>.*?</b>)", b[1]):
                rr = p.add_run(part[3:-4] if part.startswith("<b>") else part)
                rr.bold = part.startswith("<b>")
                if k in ("small", "small_kwn"):
                    rr.font.size = Pt(7.5)
        elif k == "interp":
            p = doc.add_paragraph()
            p.add_run("Interpretation: ").bold = True
            p.add_run(b[1])
        elif k in ("bullets", "small_bullets"):
            for x in b[1]:
                bp = doc.add_paragraph(x, style="List Bullet")
                if k == "small_bullets":
                    for rr in bp.runs:
                        rr.font.size = Pt(7.5)
        elif k == "box":
            tb = doc.add_table(rows=1, cols=1)
            tb.style = "Table Grid"
            c = tb.cell(0, 0)
            shade(c, TONE[b[3]])
            p = c.paragraphs[0]
            p.add_run(b[1] + ": ").bold = True
            p.add_run(b[2])
            doc.add_paragraph()
        elif k in ("figure", "recon_figure"):
            doc.add_picture(b[1]["path"], width=Cm(16.0 if k == "figure" else 17.5))
            r = doc.add_paragraph().add_run(b[2])
            r.font.size, r.italic = Pt(7.5), True
        elif k in ("kv", "kv2", "table"):
            rows = b[1]
            tb = doc.add_table(rows=len(rows), cols=len(rows[0]))
            tb.style = "Table Grid"
            for i, row in enumerate(rows):
                for j, v in enumerate(row):
                    c = tb.cell(i, j)
                    c.text = ""
                    rr = c.paragraphs[0].add_run(strip_tags(str(v)))
                    rr.font.size = Pt(8)
                    bold = (k in ("kv", "kv2") and j % 2 == 0) or (k == "table" and i == 0)
                    rr.bold = bold
                    if bold:
                        shade(c, "EEF2F6")
            if k == "table" and len(b) > 2:
                for j, w in enumerate(b[2]):
                    for i in range(len(rows)):
                        tb.cell(i, j).width = Cm(w * 17.5 / sum(b[2]))
            doc.add_paragraph()
    doc.save(path)


def to_pdf(blocks, path, footer):
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import cm
    from reportlab.platypus import Image, KeepTogether, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle

    ss = getSampleStyleSheet()
    base = ParagraphStyle("b", parent=ss["Normal"], fontName="Helvetica", fontSize=8.5, leading=10.5, alignment=TA_LEFT)
    small = ParagraphStyle("s", parent=base, fontSize=7, leading=8.5)
    cap = ParagraphStyle("cap", parent=small, fontName="Helvetica-Oblique")
    cellst = ParagraphStyle("c", parent=base, fontSize=7.5, leading=9)
    cellb = ParagraphStyle("cb", parent=cellst, fontName="Helvetica-Bold")
    tcell = ParagraphStyle("tc", parent=base, fontSize=7, leading=8.4)        # multi-column tables
    tcellb = ParagraphStyle("tcb", parent=tcell, fontName="Helvetica-Bold")
    sbl = ParagraphStyle("sbl", parent=small, leftIndent=10, firstLineIndent=-7)
    small_kwn = ParagraphStyle("skwn", parent=small, keepWithNext=1, spaceAfter=3)
    title = ParagraphStyle("t", parent=base, fontName="Helvetica-Bold", fontSize=12.5, leading=15, spaceAfter=4)
    h = ParagraphStyle("h", parent=base, fontName="Helvetica-Bold", fontSize=10, leading=12, spaceBefore=6,
                       spaceAfter=3, textColor=colors.HexColor("#1F3A5F"), keepWithNext=1)
    h2 = ParagraphStyle("h2", parent=base, fontName="Helvetica-Bold", fontSize=9, leading=11, spaceBefore=4,
                        spaceAfter=2, keepWithNext=1)
    bl = ParagraphStyle("bl", parent=base, leftIndent=10, firstLineIndent=-7)
    W = letter[0] - 3.2 * cm
    story = []

    def esc(t):
        return (str(t).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
                .replace("&lt;b&gt;", "<b>").replace("&lt;/b&gt;", "</b>"))

    for b in blocks:
        k = b[0]
        if k == "title":
            story.append(Paragraph(esc(b[1]), title))
        elif k in ("h", "h2"):
            story.append(Paragraph(esc(b[1]), h if k == "h" else h2))
        elif k in ("p", "small"):
            story += [Paragraph(esc(b[1]), base if k == "p" else small), Spacer(1, 3)]
        elif k == "small_kwn":                     # stays on the page of the flowable after it
            story.append(Paragraph(esc(b[1]), small_kwn))
        elif k == "interp":
            story += [Paragraph("<b>Interpretation:</b> " + esc(b[1]), base), Spacer(1, 3)]
        elif k in ("bullets", "small_bullets"):
            story += [Paragraph("\u2022 " + esc(x), bl if k == "bullets" else sbl) for x in b[1]] + [Spacer(1, 3)]
        elif k == "box":
            t = Table([[Paragraph(f"<b>{esc(b[1])}:</b> {esc(b[2])}",
                                  ParagraphStyle("bx", parent=base, fontSize=9.5, leading=12))]], colWidths=[W])
            t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#" + TONE[b[3]])),
                                   ("BOX", (0, 0), (-1, -1), 0.6, colors.grey),
                                   ("LEFTPADDING", (0, 0), (-1, -1), 6), ("TOPPADDING", (0, 0), (-1, -1), 5),
                                   ("BOTTOMPADDING", (0, 0), (-1, -1), 5)]))
            story += [t, Spacer(1, 5)]
        elif k == "figure":
            f = b[1]
            fw = W * 0.8
            story += [Image(f["path"], width=fw, height=fw * f["height_px"] / f["width_px"]),
                      Paragraph(esc(b[2]), cap), Spacer(1, 5)]
        elif k == "recon_figure":
            f = b[1]
            fw = W
            story += [KeepTogether([Image(f["path"], width=fw, height=fw * f["height_px"] / f["width_px"]),
                                    Paragraph(esc(b[2]), cap)]), Spacer(1, 3)]
        elif k in ("kv", "kv2", "table"):
            rows = b[1]
            if k == "table":
                data = [[Paragraph(esc(v), tcellb if i == 0 else tcell) for v in row] for i, row in enumerate(rows)]
            else:
                data = [[Paragraph(esc(v), cellb if j % 2 == 0 else cellst) for j, v in enumerate(row)]
                        for row in rows]
            if k == "kv":
                cw = [W * x for x in (0.15, 0.35, 0.15, 0.35)]
            elif k == "kv2":
                cw = [W * 0.2, W * 0.8]
            else:
                cw = [W * x / sum(b[2]) for x in b[2]]
            t = Table(data, colWidths=cw, repeatRows=1 if k == "table" else 0)
            sty = [("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#B0B7C0")), ("VALIGN", (0, 0), (-1, -1), "TOP"),
                   ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                   ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2)]
            if k == "table":
                sty.append(("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#EEF2F6")))
            else:
                sty += [("BACKGROUND", (j, 0), (j, -1), colors.HexColor("#EEF2F6")) for j in range(0, len(rows[0]), 2)]
            t.setStyle(TableStyle(sty))
            story += [t, Spacer(1, 5)]

    pages = [0]

    def on_page(canvas, doc_):
        pages[0] = max(pages[0], doc_.page)
        canvas.saveState()
        canvas.setFont("Helvetica-Oblique", 7)
        canvas.drawCentredString(letter[0] / 2, 0.8 * cm, footer)
        canvas.drawRightString(letter[0] - 1.6 * cm, 0.8 * cm, f"page {doc_.page}")
        canvas.restoreState()

    SimpleDocTemplate(path, pagesize=letter, leftMargin=1.6 * cm, rightMargin=1.6 * cm, topMargin=1.3 * cm,
                      bottomMargin=1.4 * cm, title=strip_tags(blocks[0][1]), author="CxSV pipeline",
                      invariant=1).build(story, onFirstPage=on_page, onLaterPages=on_page)
    return pages[0]


# ─── main ────────────────────────────────────────────────────────────────────
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("sample", "scores", "params", "text", "out-docx", "out-pdf", "out-figure"):
        ap.add_argument(f"--{k}", required=True)
    for k in ("cn-stats", "cn-only", "library-qc", "depth-summary", "compatibility", "date", "quarantine",
              "junctions", "cxsv-vcf", "locus-depth", "cn-segments", "genes", "out-recon-dir", "out-recon-tsv",
              "rendered-from"):
        ap.add_argument(f"--{k}", default=None)
    ap.add_argument("--caller-vcf", action="append", default=[], help="CALLER=filtered VCF (Sniffles2, cuteSV, SVIM)")
    a = ap.parse_args(argv)
    D = load(a)
    D["rendered_from"] = a.rendered_from
    t1 = [r for r in D["scores"] if r.get("Tier") == TOP]
    fig = rfig.tier1_points_figure(t1, D["P"].get("tiers") or [], a.out_figure,
                                   max_rows=int(D["T"].get("figure_max_rows", 15)))
    if fig is None:                 # no Tier 1: the report prints a note; the figure file holds the same note
        rfig.no_tier1_note(a.out_figure)
    recons = make_reconstructions(a, D, t1)
    if a.out_recon_dir:
        os.makedirs(a.out_recon_dir, exist_ok=True)   # declared output even without Tier 1
    if a.out_recon_tsv and not recons:
        import cxsv_recon as R
        with open(a.out_recon_tsv, "w") as fh:
            fh.write("\t".join(R.TSV_COLUMNS) + "\n")
    blocks = sample_blocks(D, fig, recons)
    footer = D["T"].get("footer", "Research report")
    to_docx(blocks, a.out_docx, footer)
    pages = to_pdf(blocks, a.out_pdf, footer)
    n_panels = sum(1 for x in recons if x["figure"])
    base = pages
    if recons:                       # the same report without the reconstruction section
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            base = to_pdf(sample_blocks(D, fig, None), os.path.join(td, "base.pdf"), footer)
    per = D["T"].get("reconstruction_panels_per_page", 2)
    limit = page_limit(n_panels, per, base)
    if pages > limit:
        print(f"WARNING: clinician report has {pages} pages, more than the limit {limit} "
              f"({max(2, base)} without reconstructions + 1 per {per} reconstruction panels)", file=sys.stderr)
    print(json.dumps({"sample": a.sample, "tier1": len(t1), "figure_rows": fig["n_shown"] if fig else 0,
                      "reconstructions": len(recons), "reconstruction_panels": n_panels,
                      "ambiguous": sum(1 for x in recons if x["res"].get("status") != "resolved"),
                      "pages": pages, "pages_without_reconstruction": base, "page_limit": limit}))


if __name__ == "__main__":
    main()

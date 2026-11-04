"""Curated literature set of published pathogenic germline complex SVs (roadmap item 11).

Pure functions (standard library only) shared by build_literature_set.py, the
patient-side annotation (workflow/scripts/match_literature.py) and the tests.

Junction model
--------------
A literature junction is one novel adjacency between two reference positions,
each with a SIDE in the CxSV convention (workflow/lib/cxsv_scoring.py,
record_orientation):

    '+'  the sequence ENDING at the position is retained (joined at its 3' end)
    '-'  the sequence STARTING at the position is retained

With the ends ordered (chrom1, pos1) <= (chrom2, pos2) in hs1 .fai order:

    different chromosomes  -> BND
    (+, -)                 -> DEL  (deletion-like join)
    (-, +)                 -> DUP  (tandem-duplication-like join)
    (+, +) or (-, -)       -> INV  (inversion-type join, head-to-head / tail-to-tail)

which is how the three callers report such a join (Sniffles2 DEL/DUP/INV/BND
records; docs/JUNCTIONS.md), so a literature junction can be matched to patient
and panel junctions with the item-3 matcher.

Derivative chains (Collins et al. 2017, Additional file 3) list the fragments
of each derivative chromosome 5'->3' with a strand; consecutive fragments give
one junction each (derivative_junctions).
"""
import csv
import gzip
import io
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict

SIDES = ("+", "-")
PRECISIONS = ("bp", "approximate", "region")
TELOMERES = ("pter", "qter")


# ─── minimal .xlsx reader (stdlib) ───────────────────────────────────────────
_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
_R_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"


def _col_index(ref):
    letters = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_xlsx(path):
    """{sheet name: [[cell text, ...], ...]} for every sheet of an .xlsx file.

    Values are returned as text exactly as stored (shared strings resolved,
    numbers as written in the file); empty cells are ''."""
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall("m:si", _NS):
            shared.append("".join(t.text or "" for t in si.iter("{%s}t" % _NS["m"])))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels}
    out = {}
    for sh in wb.find("m:sheets", _NS):
        t = target[sh.get(_R_NS)].lstrip("/")
        t = t if t.startswith("xl/") else "xl/" + t
        root = ET.fromstring(z.read(t))
        rows = []
        for row in root.iter("{%s}row" % _NS["m"]):
            cells = {}
            for c in row.findall("m:c", _NS):
                v = c.find("m:v", _NS)
                if c.get("t") == "s" and v is not None:
                    val = shared[int(v.text)]
                elif c.get("t") == "inlineStr":
                    val = "".join(x.text or "" for x in c.iter("{%s}t" % _NS["m"]))
                else:
                    val = v.text if v is not None else ""
                cells[_col_index(c.get("r"))] = val
            rows.append([cells.get(i, "") for i in range(max(cells) + 1)] if cells else [])
        out[sh.get("name")] = rows
    return out


# ─── tables ──────────────────────────────────────────────────────────────────
def read_tsv(path):
    """Rows of a tab-separated file with '#' comment lines and one header row."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#") and ln.strip()]
    return list(csv.DictReader(io.StringIO("".join(lines)), delimiter="\t"))


def write_tsv(path, rows, columns, header_comment=None):
    with open(path, "w", newline="") as fh:
        if header_comment:
            for line in header_comment.rstrip("\n").split("\n"):
                fh.write(f"# {line}\n" if line else "#\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        for r in rows:
            w.writerow(["." if r.get(c) in (None, "") else r.get(c) for c in columns])


def ucsc(chrom):
    """'7', 'chr7', 'X', 23 -> 'chr7' / 'chrX' (UCSC names, used for all builds)."""
    c = str(chrom).strip()
    c = c[3:] if c.lower().startswith("chr") else c
    return {"23": "chrX", "24": "chrY", "MT": "chrM", "M": "chrM"}.get(c, "chr" + c)


# ─── junction model ──────────────────────────────────────────────────────────
def flip(side):
    return "-" if side == "+" else "+"


def order_key(order, chrom, pos):
    """hs1 .fai contig rank (unknown contigs after known ones, by name), then position."""
    return (order.get(chrom, 1 << 20), chrom if chrom not in order else "", int(pos))


def oriented(c1, p1, s1, c2, p2, s2, order):
    """Ends ordered (chrom1, pos1) <= (chrom2, pos2); sides travel with their ends."""
    if order_key(order, c2, p2) < order_key(order, c1, p1):
        return c2, int(p2), s2, c1, int(p1), s1
    return c1, int(p1), s1, c2, int(p2), s2


def svtype_of(c1, p1, s1, c2, p2, s2, order):
    """Junction-site svtype of a join, from its two sides (ends ordered first)."""
    c1, p1, s1, c2, p2, s2 = oriented(c1, p1, s1, c2, p2, s2, order)
    if c1 != c2:
        return "BND"
    if s1 == s2:
        return "INV"
    return "DEL" if (s1, s2) == ("+", "-") else "DUP"


def make_junction(event_id, jid, c1, p1, s1, c2, p2, s2, precision, order, **extra):
    if s1 not in SIDES or s2 not in SIDES:
        raise ValueError(f"{event_id} J{jid}: sides must be +/- (got {s1!r}, {s2!r})")
    if precision not in PRECISIONS:
        raise ValueError(f"{event_id} J{jid}: bad precision {precision!r}")
    c1, p1, s1, c2, p2, s2 = oriented(c1, p1, s1, c2, p2, s2, order)
    t = svtype_of(c1, p1, s1, c2, p2, s2, order)
    j = {"event_id": event_id, "junction_id": f"{event_id}.J{jid}", "chrom1": c1, "pos1": p1, "side1": s1,
         "chrom2": c2, "pos2": p2, "side2": s2, "svtype": t,
         "svlen": abs(p2 - p1) if t != "BND" else None, "precision": precision}
    j.update(extra)
    return j


def _is_qter_end(f):
    return (f[4] == "+" and f[3] == "qter") or (f[4] == "-" and f[2] == "pter")


def _is_pter_start(f):
    return (f[4] == "+" and f[2] == "pter") or (f[4] == "-" and f[3] == "qter")


def fix_derivative_boundaries(ders):
    """Re-attach telomeric end fragments listed under the next derivative.

    In UTR22 (Collins et al. 2017, Additional file 3) the last, qter-bearing
    fragment of der(7) and of der(8) (AE: 7q-qter; AI: 8q-qter) carries the label
    of the following derivative, so that derivative would start with a qter end
    and the previous one would end inside a chromosome. When a derivative starts
    with a fragment that ends at a telomere and the previous derivative does not
    end at one, the fragment is moved to the end of the previous derivative.
    Returns (fixed {derivative: fragments}, [moved fragment labels])."""
    names = list(ders)
    fixed = {n: list(v) for n, v in ders.items()}
    moved = []
    for prev, cur in zip(names, names[1:]):
        f = fixed[cur]
        if f and _is_qter_end(f[0]) and not _is_pter_start(f[0]) and fixed[prev] \
                and not _is_qter_end(fixed[prev][-1]):
            fixed[prev].append(f.pop(0))
            moved.append(f"{fixed[prev][-1][0]}: {cur} -> {prev}")
    return fixed, moved


def derivative_junctions(fragments):
    """Novel adjacencies along one derivative chromosome.

    fragments: [(label, chrom, start, end, strand)] 5'->3' on the derivative,
    start/end 1-based reference coordinates ('pter'/'qter' allowed for the
    terminal fragments), strand '+' (reference orientation) or '-'.
    Returns [(left_label, right_label, chrom1, pos1, side1, chrom2, pos2, side2)];
    joins that restore a reference adjacency (consecutive '+' fragments with
    end + 1 == next start, or the '-' equivalent) are skipped."""
    out = []
    for (la, ca, sa, ea, xa), (lb, cb, sb, eb, xb) in zip(fragments, fragments[1:]):
        # leaving fragment A: '+' leaves at its end (segment ending there retained -> '+'),
        # '-' leaves at its start (segment starting there retained -> '-')
        pa, side_a = (ea, "+") if xa == "+" else (sa, "-")
        # entering fragment B: '+' enters at its start ('-'), '-' enters at its end ('+')
        pb, side_b = (sb, "-") if xb == "+" else (eb, "+")
        if pa in TELOMERES or pb in TELOMERES:
            raise ValueError(f"telomere inside a derivative between {la} and {lb}")
        pa, pb = int(pa), int(pb)
        if ca == cb and xa == xb == "+" and pb == pa + 1:
            continue
        if ca == cb and xa == xb == "-" and pa == pb + 1:
            continue
        out.append((la, lb, ca, pa, side_a, cb, pb, side_b))
    return out


def breakends(junctions, merge_bp=300):
    """Distinct breakpoint positions of a list of junction dicts (merged within merge_bp)."""
    by = defaultdict(list)
    for j in junctions:
        by[j["chrom1"]].append(int(j["pos1"]))
        by[j["chrom2"]].append(int(j["pos2"]))
    n = 0
    for c, ps in by.items():
        last = None
        for p in sorted(ps):
            if last is None or p - last > merge_bp:
                n += 1
            last = p
    return n


# ─── source parsers ──────────────────────────────────────────────────────────
def _int(x):
    x = str(x).strip().replace(",", "")
    if x in ("", "nan", "NaN", "NA"):
        return None
    return int(float(x))


NAZARYAN_NAMES = {"Patient 72": "P72", "Patient 81": "P81", "Patient 00": "P00", "Patient 06": "P06",
                  "Patient 74": "P74"}


def parse_nazaryan2018(rows, types, order, chrom_of_case=None):
    """S2 Table rows (read_xlsx sheet) + curated junction types -> (junctions, skipped).

    rows: header row then one row per junction; the case name, category and
    chromosome are given on the first junction of a case only. types:
    {(case, junction): row of curated/nazaryan2018_junction_types.tsv}."""
    head = [h.replace("\n", " ").strip() for h in rows[0]]
    col = {h: i for i, h in enumerate(head)}
    i_name, i_chr, i_j = col["Name"], col["Chromosome"], col["Junction"]
    i_s1, i_s2 = col["Side 1"], col["Side 2"]
    case, chrom = None, None
    out, skipped = [], []
    for r in rows[1:]:
        r = r + [""] * (len(head) - len(r))
        if not str(r[i_j]).strip():
            continue                       # blank rows / the duplicated block at the end of the sheet
        if str(r[i_name]).strip():
            case = NAZARYAN_NAMES.get(r[i_name].strip(), r[i_name].strip())
            chrom = ucsc(r[i_chr])
        jn = _int(r[i_j])
        t = types[(case, str(jn))]
        p1, p2 = _int(r[i_s1]), _int(r[i_s2])
        if t["svtype"] == "." or p1 is None or p2 is None:
            skipped.append({"event_id": f"NP18_{case}", "junction": jn, "reason": t["basis"] or "one_sided"})
            continue
        prec = "approximate" if "approximate" in t["basis"] else "bp"
        j = make_junction(f"NP18_{case}", jn, chrom, p1, t["side1"], chrom, p2, t["side2"], prec, order,
                          type_basis=t["basis"], evidence=f"S2 Table {case} junction {jn}")
        if j["svtype"] != t["svtype"]:
            raise ValueError(f"NP18_{case} J{jn}: sides give {j['svtype']}, curated {t['svtype']}")
        out.append(j)
    return out, skipped


def parse_collins_derivatives(sheet_rows):
    """One Additional-file-3 sheet -> ({derivative: [fragment tuples]}, [CNV dicts])."""
    head = [h.strip() for h in sheet_rows[0]]
    ders = defaultdict(list)
    cnvs = []
    i_cnv = head.index("CNV") if "CNV" in head else None
    for r in sheet_rows[1:]:
        r = r + [""] * (len(head) - len(r))
        if r[0].strip():
            der, frag, strand, chrom, s, e = (x.strip() for x in r[:6])
            ders[der].append((frag, ucsc(chrom), s if s in TELOMERES else int(float(s)),
                              e if e in TELOMERES else int(float(e)), strand))
        if i_cnv is not None and r[i_cnv].strip():
            c = r[i_cnv:]
            cn = {"label": c[0].strip(), "chrom": ucsc(str(int(float(c[1])))) if c[1].strip() not in ("X", "Y")
                  else ucsc(c[1]), "start": int(float(c[2])), "end": int(float(c[3])), "class": c[5].strip()}
            if len(c) > 6 and "±" in c[6]:
                cn["copies"] = float(c[6].split("±")[0])
            cnvs.append(cn)
    return dict(ders), cnvs


def parse_grochowski2024(rows, order):
    """DUP-TRP/INV-DUP from segment coordinates (DUP1, TRP, DUP2; hg19 array positions).

    Model (Carvalho et al. 2011; the positive-control E01 structure): the derivative
    is ... DUP1 TRP DUP2 | TRP(inverted) | DUP1 TRP DUP2 ..., i.e. two novel joins:
    junction 1 DUP2 end -> TRP end (head-to-head, in the distal inverted repeat) and
    junction 2 TRP start -> DUP1 start (tail-to-tail). Both are array positions."""
    js, cn, flags = [], [], []
    for r in rows:
        e = f"GRO24_{r['case']}"
        d1s, d1e, ts, te, d2s, d2e = (int(r[k]) for k in ("dup1_start", "dup1_end", "trp_start", "trp_end",
                                                           "dup2_start", "dup2_end"))
        if not (d1s <= d1e <= ts <= te <= d2s <= d2e):
            flags.append({"event_id": e, "flag": "segments not in DUP1<TRP<DUP2 order; model applied as given"})
        js.append(make_junction(e, 1, "chrX", te, "+", "chrX", d2e, "+", "approximate", order,
                                evidence="Table S1 TRP end / DUP2 end"))
        js.append(make_junction(e, 2, "chrX", d1s, "-", "chrX", ts, "-", "approximate", order,
                                evidence="Table S1 DUP1 start / TRP start"))
        cn += [dict(event_id=e, chrom="chrX", start=d1s, end=d1e, delta=1),
               dict(event_id=e, chrom="chrX", start=ts, end=te, delta=2),
               dict(event_id=e, chrom="chrX", start=d2s, end=d2e, delta=1)]
    return js, cn, flags


CLINVAR_TYPES = {"complex substitution": "complex_substitution", "inversion": "inversion",
                 "interchromosomal translocation": "translocation",
                 "intrachromosomal translocation": "translocation"}
CLINVAR_PATHOGENIC = ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")


def parse_clinvar_nstd102(path, order):
    """dbVar nstd102 (ClinVar) GRCh38 variant calls -> (events, junctions, counts).

    Kept: germline-or-unknown-origin records of type complex substitution,
    inversion or translocation with a Pathogenic / Likely pathogenic
    clinical_significance and numeric coordinates. Region-level:
      complex substitution -> event region only (no junction),
      inversion            -> one INV junction start..stop (precision 'region'),
      translocation        -> one BND junction between from_ and to_ coordinates
                              when both are given (precision 'region'), else region only.
    Translocation orientation is not given, so sides are written '+'/'-' placeholders
    and only the svtype (BND) is used by matching."""
    counts = defaultdict(int)
    events, junctions = [], []
    with gzip.open(path, "rt") as fh:
        header = None
        for line in fh:
            if line.startswith("#variant_call_accession"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#") or header is None:
                continue
            d = dict(zip(header, line.rstrip("\n").split("\t")))
            t = CLINVAR_TYPES.get(d["variant_call_type"])
            if not t:
                continue
            counts[f"{t}: all"] += 1
            if d.get("clinical_significance") not in CLINVAR_PATHOGENIC:
                continue
            counts[f"{t}: pathogenic/likely"] += 1
            if "somatic" in d.get("origin", "").lower():
                counts[f"{t}: somatic dropped"] += 1
                continue
            start = _int(d.get("start") or d.get("outer_start") or d.get("inner_start") or "")
            stop = _int(d.get("stop") or d.get("outer_stop") or d.get("inner_stop") or "")
            if not d.get("chr") or start is None or stop is None:
                counts[f"{t}: no coordinates"] += 1
                continue
            e = f"CV_{d['variant_call_accession']}"
            chrom = ucsc(d["chr"])
            events.append({"event_id": e, "source_id": "clinvar_nstd102", "case_id": d["variant_call_id"],
                           "source_table": "nstd102.GRCh38.variant_call.tsv.gz", "build": "GRCh38",
                           "phenotype": d.get("description", "").split(" AND ", 1)[-1][:120],
                           "genes_reported": ".", "inheritance": d.get("origin") or "not stated",
                           "author_event_type": d["variant_call_type"],
                           "author_classification": d["clinical_significance"],
                           "carrier_sex": "female", "zygosity": d.get("zygosity") or "not stated",
                           "region_chrom": chrom, "region_start": start, "region_end": stop,
                           "notes": "region-level ClinVar record"})
            counts[f"{t}: kept"] += 1
            if t == "inversion" and stop > start:
                junctions.append(make_junction(e, 1, chrom, start, "+", chrom, stop, "+", "region", order,
                                               evidence="nstd102 start/stop"))
            elif t == "translocation" and d.get("from_chr") and d.get("to_chr") and _int(d.get("from_coord") or "") \
                    and _int(d.get("to_coord") or "") and ucsc(d["from_chr"]) != ucsc(d["to_chr"]):
                junctions.append(make_junction(e, 1, ucsc(d["from_chr"]), _int(d["from_coord"]), "+",
                                               ucsc(d["to_chr"]), _int(d["to_coord"]), "-", "region", order,
                                               evidence="nstd102 from/to"))
    return events, junctions, dict(counts)

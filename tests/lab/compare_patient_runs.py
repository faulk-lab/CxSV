#!/usr/bin/env python3
"""Compare two patient-workflow result trees record by record (bundle validation, roadmap item 2).

    compare_patient_runs.py <results_root_A> <results_root_B> [--out table.tsv]
        [--alias-a PATH=TOKEN]... [--alias-b PATH=TOKEN]...

Every file below either root is compared, except run logs (logs/, *.log,
.snakemake/):
  VCF   records only ('##' meta lines carry dates and command lines)
  PDF   after dropping /CreationDate, /ModDate and /ID (render timestamps)
  .snf  Sniffles2 SNF: header without its run fields (config: paths, command,
        start date, run id) + the DECOMPRESSED blocks (each gzip member carries
        its own mtime)
  other byte for byte; if that differs, again after replacing each side's
        --alias paths by their tokens (e.g. the results root /work/results/patient
        inside the container vs the host path), reported as
        "identical after path normalisation"
Two run records are expected to differ and are listed as such: panel_release_check.txt
(release path, runtime image) and panel_compatibility.txt (release path); their
status lines must still agree. The clinician report (<sample>.clinician_report.docx/.pdf,
panel-v1.3) names the run's identity in one row ("Pipeline bundle": the bundle id, or "not
run from an offline bundle" for a run outside a bundle; "Report date ...; code <commit12>"),
so a bundle run and a repository run differ there by design: the report TEXT (docx: every
paragraph, table cell and footer; PDF: every text string of the decompressed page streams,
in order) is compared with those identity values masked, and an equal text is listed as an
expected run-record difference. Exit 1 when any other file differs or exists on one side only.
"""
import argparse
import base64
import json
import os
import re
import sys
import zipfile
import zlib

EXPECTED = {"panel_release_check.txt": "status:", "panel_compatibility.txt": "status:"}
PDF_VOLATILE = re.compile(rb"^(/CreationDate|/ModDate|/ID|<< /CreationDate).*$", re.M)
REPORT = re.compile(r"\.clinician_report\.(docx|pdf)$")
# run identity in the clinician report: bundle id or "not run from an offline bundle"; "; code <commit12>"
REPORT_IDENTITY = [(re.compile(r"cxsv-bundle-\d{8}-[0-9a-f]{12}|not run from an offline bundle"), "<RUN-IDENTITY>"),
                   (re.compile(r";?\s*code [0-9a-f]{12}\b"), "")]


def mask_identity(text):
    for rx, rep in REPORT_IDENTITY:
        text = rx.sub(rep, text)
    return text


def docx_text(path):
    """Text of every w:t run of the document, headers and footers, paragraph by paragraph."""
    out = []
    with zipfile.ZipFile(path) as z:
        for n in sorted(z.namelist()):
            if re.fullmatch(r"word/(document|header\d*|footer\d*)\.xml", n):
                xml = z.read(n).decode("utf-8")
                for par in re.findall(r"<w:p[ >].*?</w:p>", xml, re.S):
                    out.append("".join(re.findall(r"<w:t(?: [^>]*)?>([^<]*)</w:t>", par)))
    return "\n".join(out)


def pdf_text(path):
    """Literal text strings of every FlateDecode stream, in order (layout operators dropped)."""
    data, out = open(path, "rb").read(), []
    for m in re.finditer(rb"stream\r?\n(.*?)endstream", data, re.S):
        body = m.group(1).strip()
        if body.rstrip().endswith(b"~>"):                 # reportlab: /ASCII85Decode then /FlateDecode
            try:
                body = body.strip()
                body = base64.a85decode(body if body.startswith(b"<~") else body[:-2], adobe=body.startswith(b"<~"))
            except ValueError:
                continue
        try:
            body = zlib.decompress(body)
        except zlib.error:
            pass
        if b"BT" not in body:                               # page content streams only (text objects)
            continue
        out += [x.decode("latin-1") for x in re.findall(rb"\(((?:\\.|[^\\)])*)\)\s*Tj", body)]
    return "\n".join(out)


def report_text(path):
    """Masked report text; None when no text could be extracted (never counts as equal)."""
    t = docx_text(path) if path.endswith(".docx") else pdf_text(path)
    return mask_identity(t) if t.strip() else None


def files(root):
    out = set()
    for d, dirs, fs in os.walk(root):
        dirs[:] = [x for x in dirs if x not in ("logs", ".snakemake")]
        for f in fs:
            if not f.endswith(".log"):
                out.add(os.path.relpath(os.path.join(d, f), root))
    return out


def snf_content(data):
    i = data.index(b"\n")
    head = json.loads(data[:i])
    head.pop("config", None)
    rest, blocks = data[i + 1:], []
    while rest:
        d = zlib.decompressobj(31)
        blocks.append(d.decompress(rest))
        rest = d.unused_data
    return json.dumps(head, sort_keys=True).encode() + b"\0" + b"\0".join(blocks)


def content(path):
    data = open(path, "rb").read()
    if path.endswith(".snf"):
        return snf_content(data), "snf header minus run fields + decompressed blocks"
    if path.endswith(".vcf"):
        return b"\n".join(l for l in data.split(b"\n") if not l.startswith(b"##")), "records"
    if path.endswith(".pdf"):
        return PDF_VOLATILE.sub(b"", re.sub(rb"/CreationDate \(D:[^)]*\)|/ModDate \(D:[^)]*\)", b"", data)), \
            "pdf minus timestamps"
    return data, "bytes"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--out")
    ap.add_argument("--alias-a", action="append", default=[], help="PATH=TOKEN for tree A")
    ap.add_argument("--alias-b", action="append", default=[], help="PATH=TOKEN for tree B")
    a = ap.parse_args()

    def aliases(pairs):
        out = [tuple(x.split("=", 1)) for x in pairs]
        return sorted(((p.encode(), t.encode()) for p, t in out), key=lambda x: -len(x[0]))

    def normalise(data, pairs):
        for path, token in pairs:
            data = data.replace(path, token)
        return data
    al_a, al_b = aliases(a.alias_a), aliases(a.alias_b)
    fa, fb = files(a.a), files(a.b)
    rows, bad = [], 0
    for rel in sorted(fa | fb):
        if rel not in fa or rel not in fb:
            rows.append((rel, "only in " + ("A" if rel in fa else "B"), "-"))
            bad += 1
            continue
        ca, how = content(os.path.join(a.a, rel))
        cb, _ = content(os.path.join(a.b, rel))
        base = os.path.basename(rel)
        if ca == cb:
            rows.append((rel, "identical", how))
        elif (al_a or al_b) and normalise(ca, al_a) == normalise(cb, al_b) and base not in EXPECTED:
            rows.append((rel, "identical after path normalisation", how))
        elif REPORT.search(base) and (report_text(os.path.join(a.a, rel)) or 0) == (report_text(os.path.join(a.b, rel)) or 1):
            rows.append((rel, "differs as expected (run identity row); report text identical", "report text, identity masked"))
        elif base in EXPECTED:
            sa = [l for l in ca.decode().splitlines() if l.startswith(EXPECTED[base])]
            sb = [l for l in cb.decode().splitlines() if l.startswith(EXPECTED[base])]
            ok = sa == sb and sa
            rows.append((rel, "differs as expected (paths/runtime); same status" if ok else "STATUS DIFFERS", how))
            bad += 0 if ok else 1
        else:
            rows.append((rel, "DIFFERS", how))
            bad += 1
    text = "file\tresult\tcompared\n" + "".join("\t".join(r) + "\n" for r in rows)
    if a.out:
        open(a.out, "w").write(text)
    n_id = sum(r[1] == "identical" for r in rows)
    n_norm = sum(r[1] == "identical after path normalisation" for r in rows)
    print(text, end="")
    print(f"# {len(rows)} files: {n_id} identical, {n_norm} identical after path normalisation, "
          f"{len(rows) - n_id - n_norm - bad} expected run-record differences, {bad} unexpected differences")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

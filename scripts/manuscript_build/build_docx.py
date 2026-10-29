"""Assemble the Human Genomics manuscript (DOCX) from the markdown parts.

- citation keys [@a; @b] -> BMC numbered citations [1, 2] / [1–3], numbered by first appearance
- reference list in BMC style from references.json
- double spacing, continuous line numbers, page numbers (BMC formatting requirements)
Usage: python build_docx.py   (run from submission/build)
Outputs: ../Manuscript_HumanGenomics.docx, ../Cover_letter.docx, manuscript_numbered.md
"""
import json, re
from pathlib import Path
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, Cm

B = Path(__file__).resolve().parent
OUT = B / "output"
OUT.mkdir(exist_ok=True)
PARTS = ["ms_front.md", "ms_background.md", "ms_methods.md", "ms_results.md", "ms_discussion.md", "ms_back.md"]
TAIL = ["ms_figure_legends.md", "ms_tables.md", "ms_additional_files.md"]

refs = {r["key"]: r for r in json.load(open(B / "references.json", encoding="utf-8"))}


# ---------------------------------------------------------------- citations
def number_citations(text, order=None):
    order = [] if order is None else order
    def repl(m):
        keys = re.findall(r"@([\w\-]+)", m.group(1))
        nums = []
        for k in keys:
            if k not in refs:
                raise KeyError(f"unknown citation key {k}")
            if k not in order:
                order.append(k)
            nums.append(order.index(k) + 1)
        nums = sorted(set(nums))
        out, i = [], 0
        while i < len(nums):            # compress runs of >= 3 into ranges
            j = i
            while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
                j += 1
            out.append(f"{nums[i]}–{nums[j]}" if j - i >= 2 else ", ".join(str(n) for n in nums[i:j + 1]))
            i = j + 1
        return "[" + ", ".join(out) + "]"
    return re.sub(r"\[(@[^\]]+)\]", repl, text), order


BMC_TXT = {}
for line in (B / "references_bmc.txt").read_text(encoding="utf-8").splitlines():
    m = re.match(r"^\[([\w\-]+)\]\s+(.*)$", line)
    if m:
        BMC_TXT[m.group(1)] = m.group(2)


def short_pages(m):           # BMC style: 2562–2573 -> 2562–73 (keep >= 2 digits of the end page)
    a, b = m.group(1), m.group(2)
    if len(a) == len(b) and a.isdigit() and b.isdigit():
        k = 0
        while k < len(a) - 2 and a[k] == b[k]:
            k += 1
        b = b[k:]
    return f":{a}–{b}."


def bmc_ref(r):
    txt = BMC_TXT.get(r["key"])
    if txt is None:
        raise KeyError(f"no formatted reference for {r['key']}")
    return re.sub(r":(\d+)–(\d+)\.", short_pages, txt, count=1)


# ---------------------------------------------------------------- docx helpers
def set_base(doc):
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"; st.font.size = Pt(12)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    pf = st.paragraph_format; pf.line_spacing = 2.0; pf.space_after = Pt(0)
    for lvl, size in ((1, 14), (2, 12), (3, 12)):
        h = doc.styles[f"Heading {lvl}"]
        h.font.name = "Times New Roman"; h.font.size = Pt(size); h.font.bold = True
        h.font.italic = lvl == 3; h.font.color.rgb = None
        h.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
        h.paragraph_format.line_spacing = 2.0; h.paragraph_format.space_before = Pt(12); h.paragraph_format.space_after = Pt(0)
    sec = doc.sections[0]
    sec.left_margin = sec.right_margin = Cm(2.5); sec.top_margin = sec.bottom_margin = Cm(2.5)


def add_line_numbers(doc):
    for sec in doc.sections:
        ln = OxmlElement("w:lnNumType")
        ln.set(qn("w:countBy"), "1"); ln.set(qn("w:restart"), "continuous"); ln.set(qn("w:distance"), "360")
        sec._sectPr.append(ln)


def add_page_numbers(doc):
    p = doc.sections[0].footer.paragraphs[0]; p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    for tag, text in (("begin", None), (None, "PAGE"), ("end", None)):
        if tag:
            el = OxmlElement("w:fldChar"); el.set(qn("w:fldCharType"), tag)
        else:
            el = OxmlElement("w:instrText"); el.set(qn("xml:space"), "preserve"); el.text = text
        run._r.append(el)


INLINE = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)")


def add_runs(par, text):
    for tok in INLINE.split(text):
        if not tok:
            continue
        if tok.startswith("**") and tok.endswith("**"):
            par.add_run(tok[2:-2]).bold = True
        elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2:
            par.add_run(tok[1:-1]).italic = True
        elif tok.startswith("`") and tok.endswith("`"):
            par.add_run(tok[1:-1])
        else:
            par.add_run(tok.replace("\\*", "*"))


def add_table(doc, rows):
    cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
    cells = [r for r in cells if not all(re.fullmatch(r":?-{2,}:?", c) for c in r)]
    ncol = max(len(r) for r in cells)
    t = doc.add_table(rows=len(cells), cols=ncol); t.style = "Table Grid"
    for i, r in enumerate(cells):
        for j in range(ncol):
            c = t.cell(i, j); c.text = ""
            p = c.paragraphs[0]; p.paragraph_format.line_spacing = 1.0
            add_runs(p, (r[j] if j < len(r) else "").replace("\\|", "|"))
            for run in p.runs:
                run.font.size = Pt(9); run.bold = run.bold or i == 0
    doc.add_paragraph()


def md_to_docx(doc, md):
    lines = md.splitlines(); i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("|"):
            block = []
            while i < len(lines) and lines[i].startswith("|"):
                block.append(lines[i]); i += 1
            add_table(doc, block); continue
        if ln.strip() == "\\newpage":
            doc.add_page_break(); i += 1; continue
        m = re.match(r"^(#{1,3})\s+(.*)", ln)
        if m:
            lvl = len(m.group(1))
            if lvl == 1:
                p = doc.add_paragraph(); r = p.add_run(m.group(2)); r.bold = True; r.font.size = Pt(16)
            else:
                doc.add_heading(m.group(2), level=lvl - 1)
        elif re.match(r"^\s*(\d+\.|-)\s+", ln):
            p = doc.add_paragraph(style="List Number" if re.match(r"^\s*\d+\.", ln) else "List Bullet")
            add_runs(p, re.sub(r"^\s*(\d+\.|-)\s+", "", ln))
        elif ln.strip():
            add_runs(doc.add_paragraph(), ln)
        i += 1


# ---------------------------------------------------------------- manuscript
body = "\n\n".join((B / f).read_text(encoding="utf-8") for f in PARTS)
body, order = number_citations(body)
reflist = "## References\n\n" + "\n".join(f"{n}. {bmc_ref(refs[k])}" for n, k in enumerate(order, 1))
tail = "\n\n".join((B / f).read_text(encoding="utf-8") for f in TAIL if (B / f).exists())
n_main = len(order)
tail, order_all = number_citations(tail, list(order)) if "[@" in tail else (tail, order)
extra = order_all[n_main:]
if extra:
    raise SystemExit("figure legends / tables cite new references; move citations into the main text first")
full = body + "\n\n" + reflist + "\n\n" + tail
(B / "manuscript_numbered.md").write_text(full, encoding="utf-8")

doc = Document(); set_base(doc)
md_to_docx(doc, full)
add_line_numbers(doc); add_page_numbers(doc)
doc.save(OUT / "Manuscript_HumanGenomics.docx")
print(f"manuscript: {len(order)} references cited; saved {OUT / 'Manuscript_HumanGenomics.docx'}")

# ---------------------------------------------------------------- cover letter (single spacing, no line numbers)
cl = Document(); set_base(cl)
cl.styles["Normal"].paragraph_format.line_spacing = 1.15
cl.styles["Normal"].paragraph_format.space_after = Pt(6)
md_to_docx(cl, (B / "cover_letter.md").read_text(encoding="utf-8"))
cl.save(OUT / "Cover_letter.docx")
print("cover letter saved")

# ---------------------------------------------------------------- response to reviewers (single spacing)
rr = Document(); set_base(rr)
rr.styles["Normal"].paragraph_format.line_spacing = 1.15
rr.styles["Normal"].paragraph_format.space_after = Pt(6)
resp = "\n\n".join((B / f"response_part{i}.md").read_text(encoding="utf-8") for i in range(1, 5))
md_to_docx(rr, resp)
rr.save(OUT / "Response_to_Reviewers.docx")
print("response letter saved")

rr2 = Document(); set_base(rr2)
rr2.styles["Normal"].paragraph_format.line_spacing = 1.15
rr2.styles["Normal"].paragraph_format.space_after = Pt(6)
md_to_docx(rr2, (B / "response_round2.md").read_text(encoding="utf-8"))
rr2.save(OUT / "Response_to_Reviewers_Round2.docx")
print("round-2 response letter saved")

rr4 = Document(); set_base(rr4)
rr4.styles["Normal"].paragraph_format.line_spacing = 1.15
rr4.styles["Normal"].paragraph_format.space_after = Pt(6)
md_to_docx(rr4, (B / "response_round4.md").read_text(encoding="utf-8"))
rr4.save(OUT / "Response_to_Reviewer_Round4.docx")
print("round-4 response letter saved")

rr5 = Document(); set_base(rr5)
rr5.styles["Normal"].paragraph_format.line_spacing = 1.15
rr5.styles["Normal"].paragraph_format.space_after = Pt(6)
md_to_docx(rr5, (B / "response_round5.md").read_text(encoding="utf-8"))
rr5.save(OUT / "Response_to_Reviewers_Round5.docx")
print("round-5 response letter saved")

rr6 = Document(); set_base(rr6)
rr6.styles["Normal"].paragraph_format.line_spacing = 1.15
rr6.styles["Normal"].paragraph_format.space_after = Pt(6)
md_to_docx(rr6, (B / "response_round6.md").read_text(encoding="utf-8"))
rr6.save(OUT / "Response_to_External_Reviewer_Round6.docx")
print("round-6 response letter saved")

rr7 = Document(); set_base(rr7)
rr7.styles["Normal"].paragraph_format.line_spacing = 1.15
rr7.styles["Normal"].paragraph_format.space_after = Pt(6)
md_to_docx(rr7, (B / "response_round7.md").read_text(encoding="utf-8"))
rr7.save(OUT / "Response_to_External_Reviewer_Round7.docx")
print("round-7 response letter saved")

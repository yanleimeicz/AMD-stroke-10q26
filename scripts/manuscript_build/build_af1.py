"""Additional file 1 (DOCX): Supplementary Methods + Figures S1-S7 with titles and legends.
Also writes ms_figure_legends.md (main-figure titles and legends for the manuscript) from figure_legends_final.md.
Usage: python build_af1.py   (run from submission/build)
"""
import re
from pathlib import Path
from docx import Document
from docx.shared import Pt, Cm
from docx.oxml.ns import qn

B = Path(__file__).resolve().parent
SUB = B
LEG = (B / "figure_legends_final.md").read_text(encoding="utf-8")

# parse "### Figure X. Title" blocks
blocks = {}
for m in re.finditer(r"^### Figure (S?\d+)\. (.+?)\n(.*?)(?=^### Figure |\Z)", LEG, flags=re.S | re.M):
    blocks[m.group(1)] = (m.group(2).strip(), m.group(3).strip())

# ---- main-figure legends for the manuscript (BMC: "Fig. 1 Title" + legend, after the references)
main = ["\\newpage", "", "## Figure legends", ""]
for i in range(1, 6):
    t, l = blocks[str(i)]
    main += [f"**Fig. {i}** {t.rstrip('.')}", "", l, ""]
(B / "ms_figure_legends.md").write_text("\n".join(main), encoding="utf-8")


def runs(par, text, size=None):
    for tok in re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*)", text):
        if not tok:
            continue
        if tok.startswith("**"):
            r = par.add_run(tok[2:-2]); r.bold = True
        elif tok.startswith("*") and len(tok) > 2:
            r = par.add_run(tok[1:-1]); r.italic = True
        else:
            r = par.add_run(tok)
        if size:
            r.font.size = Pt(size)


doc = Document()
st = doc.styles["Normal"]; st.font.name = "Times New Roman"; st.font.size = Pt(11)
st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
st.paragraph_format.line_spacing = 1.15; st.paragraph_format.space_after = Pt(6)
for lvl in (1, 2):
    h = doc.styles[f"Heading {lvl}"]; h.font.name = "Times New Roman"; h.font.color.rgb = None
    h.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
sec = doc.sections[0]; sec.left_margin = sec.right_margin = Cm(2.0)

# ---- citations in the Supplementary Note use the numbering of the main text
PARTS = ["ms_front.md", "ms_background.md", "ms_methods.md", "ms_results.md", "ms_discussion.md", "ms_back.md"]
_main = "\n\n".join((B / f).read_text(encoding="utf-8") for f in PARTS)
_order = []
for m in re.finditer(r"\[(@[^\]]+)\]", _main):
    for k in re.findall(r"@([\w\-]+)", m.group(1)):
        if k not in _order:
            _order.append(k)


_supp_order = []


def cite(text):
    def repl(m):
        nums, snums = [], []
        for k in re.findall(r"@([\w\-]+)", m.group(1)):
            if k in _order:
                nums.append(_order.index(k) + 1)
            else:
                if k not in _supp_order:
                    _supp_order.append(k)
                snums.append(_supp_order.index(k) + 1)
        parts = [str(n) for n in sorted(set(nums))] + [f"S{n}" for n in sorted(set(snums))]
        return "[" + ", ".join(parts) + "]"
    return re.sub(r"\[(@[^\]]+)\]", repl, text)


for ln in cite((B / "supp_methods.md").read_text(encoding="utf-8")).splitlines():
    if ln.startswith("# "):
        p = doc.add_paragraph(); r = p.add_run(ln[2:]); r.bold = True; r.font.size = Pt(16)
    elif ln.startswith("## "):
        doc.add_heading(ln[3:], level=1)
    elif ln.startswith("### "):
        doc.add_heading(ln[4:], level=2)
    elif ln.strip():
        runs(doc.add_paragraph(), ln)

doc.add_page_break()
doc.add_heading("Supplementary Figures", level=1)
for i in range(1, 12):
    key = f"S{i}"
    t, l = blocks[key]
    img = SUB / "figures" / "supp" / f"FigS{i}.png"
    if i > 1:
        doc.add_page_break()
    p = doc.add_paragraph(); r = p.add_run(f"Figure {key}. "); r.bold = True; r2 = p.add_run(t); r2.bold = True
    doc.add_picture(str(img), width=Cm(16.5))
    runs(doc.add_paragraph(), cite(l), size=10)

if _supp_order:
    import json as _json
    _refs = {r["key"]: r for r in _json.load(open(B / "references.json", encoding="utf-8"))}
    _txt = {}
    for _ln in (B / "references_bmc.txt").read_text(encoding="utf-8").splitlines():
        _m = re.match(r"^\[([\w\-]+)\]\s+(.*)$", _ln)
        if _m:
            _txt[_m.group(1)] = _m.group(2)
    doc.add_page_break()
    doc.add_heading("Supplementary references", level=1)
    for _i, _k in enumerate(_supp_order, 1):
        runs(doc.add_paragraph(), f"S{_i}. " + _txt[_k], size=10)
(SUB / "additional_files").mkdir(exist_ok=True)
doc.save(SUB / "additional_files" / "Additional_file_1_Supplementary_Methods_and_Figures.docx")
print("Additional file 1 saved;", len(blocks), "legend blocks")

"""Self-check: list grey (achromatic, not black/white) text spans in the main-figure PDFs."""
import fitz
from rev_common import OUT
for k in range(1, 6):
    doc = fitz.open(OUT / f"Fig{k}.pdf"); bad = {}; cols = {}
    for page in doc:
        for b in page.get_text("dict")["blocks"]:
            for l in b.get("lines", []):
                for s in l["spans"]:
                    c = s["color"]; r, g, bb = (c >> 16) & 255, (c >> 8) & 255, c & 255
                    cols[f"#{c:06x}"] = cols.get(f"#{c:06x}", 0) + 1
                    if max(r, g, bb) - min(r, g, bb) < 10 and 0 < max(r, g, bb) < 243:
                        bad[s["text"]] = f"#{c:06x}"
    print(f"Fig{k}: span colors {cols}; grey spans: {list(bad.items())[:8]}")

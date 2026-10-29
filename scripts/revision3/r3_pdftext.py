import fitz, sys
src, dst = sys.argv[1], sys.argv[2]
d = fitz.open(src)
with open(dst, 'w', encoding='utf-8') as f:
    for i, p in enumerate(d):
        f.write(f"\n=== page {i+1}\n"); f.write(p.get_text())
print(len(d))

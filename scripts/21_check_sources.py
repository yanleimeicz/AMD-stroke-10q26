import re, sys, requests
w = lambda s: sys.stdout.buffer.write((s + "\n").encode("utf-8"))
u = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003219/"
w("IAMDGC: " + " ".join(re.findall(r'href="([^"?/][^"]*)"', requests.get(u, timeout=120).text)))
for hu in [u + "harmonised/"]:
    try: w("  harmonised: " + " ".join(re.findall(r'href="([^"?/][^"]*)"', requests.get(hu, timeout=120).text)))
    except Exception as e: w(str(e))
for url in ["https://open.win.ox.ac.uk/ukbiobank/big40/release2/IDPs.html", "https://open.win.ox.ac.uk/ukbiobank/big40/BIG40-IDPs_v4/IDPs.html"]:
    try:
        t = requests.get(url, timeout=120).text
        rows = re.findall(r"<tr>(.*?)</tr>", t, flags=re.S)
        hits = [re.sub(r"<[^>]+>", " | ", r) for r in rows if "WMH" in r or "BIANCA" in r]
        w(f"{url}: {len(rows)} rows; WMH rows:"); [w("   " + " ".join(h.split())[:200]) for h in hits[:10]]
    except Exception as e: w(f"{url} ERR {e}")

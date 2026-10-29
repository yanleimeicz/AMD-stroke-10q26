"""Look up GWAS Catalog accessions (with retry) and list summary-statistics FTP files."""
import sys, time, re, requests
API = "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId"
FTP = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
w = lambda s: sys.stdout.buffer.write((s + "\n").encode("utf-8"))


def get_json(url, **kw):
    for i in range(6):
        try:
            return requests.get(url, timeout=120, **kw).json()
        except Exception:
            time.sleep(5 * (i + 1))
    return {}


def ftp_dir(acc):
    n = int(acc[4:]); lo = (n - 1) // 1000 * 1000 + 1
    return f"{FTP}/GCST{lo:08d}-GCST{lo + 999:08d}/{acc}/" if n >= 90000000 else f"{FTP}/GCST{lo:06d}-GCST{lo + 999:06d}/{acc}/"


def list_ftp(url):
    try:
        html = requests.get(url, timeout=120).text
        return [h for h in re.findall(r'href="([^"?/][^"]*)"', html)]
    except Exception as e:
        return [f"ERR {e}"]


for name, pmid in {"WMH (Sargurupremraj 2020)": 32358547, "IAMDGC AMD (Fritsche 2016)": 26691988}.items():
    st = get_json(API, params={"pubmedId": pmid, "size": 200}).get("_embedded", {}).get("studies", [])
    w(f"\n=== {name}: {len(st)} studies")
    for s in st:
        w(f"{s['accessionId']} | ss={s.get('fullPvalueSet')} | {s.get('diseaseTrait', {}).get('trait', '')[:60]} | {s.get('initialSampleSize', '')[:100]}")
    time.sleep(2)

for acc in sys.argv[1:]:
    u = ftp_dir(acc); w(f"\n--- {acc} {u}")
    for f in list_ftp(u):
        w("   " + f)
        if f.endswith("/") and "harmon" in f:
            for g in list_ftp(u + f):
                w("      " + g)

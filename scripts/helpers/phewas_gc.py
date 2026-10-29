import json, urllib.request, time, csv, sys
SNPS = ["rs4611812", "rs6795735", "rs17676309", "rs4607103"]
rows = []
for s in SNPS:
    url = f"https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id={s}&size=200&page=0"
    while url:
        for k in range(5):
            try:
                d = json.load(urllib.request.urlopen(url, timeout=90)); break
            except Exception as e:
                time.sleep(5)
        for a in d.get("_embedded", {}).get("associations", []):
            ea = ";".join(x.get("effect_allele", "") for x in a.get("snp_allele", []) if x.get("rs_id") == s)
            rows.append(dict(snp=s, effect_allele=ea, p=a.get("p_value"), beta=a.get("beta"), or_value=a.get("or_value"),
                             ci=a.get("range"), risk_freq=a.get("risk_frequency"),
                             efo=";".join(t["efo_trait"] for t in a.get("efo_traits", [])),
                             reported=";".join(a.get("reported_trait", [])), accession=a.get("accession_id"),
                             pmid=a.get("pubmed_id"), author=a.get("first_author"), location=";".join(a.get("locations", []))))
        nxt = d.get("_links", {}).get("next", {}).get("href")
        url = nxt
    print(s, sum(r["snp"] == s for r in rows), file=sys.stderr)
with open("gwascatalog_phewas.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)

"""Round-4 item C: Pan-UK Biobank (EUR) 3p14 region fetch (remote tabix via HTTP range; GRCh37) + UKB British LD (PolyFun npz),
then SuSiE/coloc inputs. rs4611812 GRCh37 = 3:64,699,445 (GRCh38 64,713,769; offset 14,324 verified with rs6795735 and rs6445425
via grch37.rest.ensembl.org). Outputs: derived/panukb3p14/region_<tag>.tsv, derived/panukb3p14/ld/*, results/revision4/c_input/*"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import csv, gzip, os, struct, sys, time, zlib
import numpy as np, pandas as pd, requests
from scipy import sparse
from joblib import Parallel, delayed
SES = requests.Session(); SES.trust_env = False          # stale host proxy 127.0.0.1:7890

P = PROJ_ROOT; D = f"{P}/derived/panukb3p14"; O = f"{P}/results/revision4"; I = f"{O}/c_input"
for d in (D, f"{D}/ld", O, I): os.makedirs(d, exist_ok=True)
BASE = "https://pan-ukb-us-east-1.s3.amazonaws.com"
CHR, START, END = "3", 64_200_000, 65_200_000
PHENOS = [("AMD_phe362.2", "phecode", "362.2", ""), ("AMD_sr6148_5", "categorical", "6148", "5"), ("T2D_phe250", "phecode", "250", "")]

def get(url, rng=None, tries=6):
    for t in range(tries):
        try:
            hd = {"Range": f"bytes={rng[0]}-{rng[1]}"} if rng else {}
            r = SES.get(url, headers=hd, timeout=120); r.raise_for_status(); return r.content
        except Exception as e:
            print("retry", url, e, flush=True); time.sleep(3 * (t + 1))
    raise RuntimeError(url)

def tbi_offsets(tbi_bytes, chrom, beg):
    b = gzip.decompress(tbi_bytes); assert b[:4] == b"TBI\1"
    n_ref, fmt, cs, cb, ce, meta, skip, l_nm = struct.unpack("<8i", b[4:36])
    names = b[36:36 + l_nm].split(b"\0")[:n_ref]; p = 36 + l_nm
    for i in range(n_ref):
        n_bin, = struct.unpack("<i", b[p:p + 4]); p += 4
        for _ in range(n_bin):
            binid, n_chunk = struct.unpack("<Ii", b[p:p + 8]); p += 8 + 16 * n_chunk
        n_intv, = struct.unpack("<i", b[p:p + 4]); p += 4
        ioff = struct.unpack(f"<{n_intv}Q", b[p:p + 8 * n_intv]); p += 8 * n_intv
        if names[i].decode() in (chrom, "chr" + chrom):
            return next(o for o in ioff[beg >> 14:] if o > 0)
    raise KeyError(chrom)

def fetch(tag, url, tbi_url):
    out = f"{D}/region_{tag}.tsv"
    if os.path.exists(out): return tag, -1
    voff = tbi_offsets(get(tbi_url), CHR, START); coff, uoff = voff >> 16, voff & 0xFFFF
    d0 = get(url, (0, 200_000)); bs = struct.unpack("<H", d0[16:18])[0] + 1
    header = zlib.decompress(d0[:bs], 31).split(b"\n")[0].decode()
    start, text, done, lines, first = coff, b"", False, [], True
    while not done:
        data = get(url, (start, start + 4_000_000 - 1)); blocks, used = [], 0
        while used < len(data) and len(data) - used >= 18:
            bsize = struct.unpack("<H", data[used + 16:used + 18])[0] + 1
            if used + bsize > len(data): break
            blocks.append(zlib.decompress(data[used:used + bsize], 31)); used += bsize
        if first and blocks: blocks[0] = blocks[0][uoff:]; first = False
        text += b"".join(blocks); start += used
        rows = text.split(b"\n"); text = rows[-1]
        for r in rows[:-1]:
            x = r.split(b"\t", 2)
            if x[0].decode() not in (CHR, "chr" + CHR): done = True; break
            pos = int(x[1])
            if pos > END: done = True; break
            if pos >= START: lines.append(r.decode())
        if used == 0: break
    with open(out, "w") as f: f.write(header + "\n" + "\n".join(lines) + "\n")
    return tag, len(lines)

man = list(csv.reader(gzip.open(f"{P}/derived/panukb/phenotype_manifest.tsv.bgz", "rt"), delimiter="\t")); h = man[0]; ix = {k: i for i, k in enumerate(h)}
to_url = lambda s: BASE + "/" + s.split("pan-ukb-us-east-1/", 1)[1]
sel = []
for tag, tt, pc, cd in PHENOS:
    x = next(r for r in man[1:] if r[ix["trait_type"]] == tt and r[ix["phenocode"]] == pc and r[ix["coding"]] == cd and r[ix["pheno_sex"]] == "both_sexes")
    sel.append(dict(tag=tag, trait_type=tt, phenocode=pc, description=x[ix["description"]], n_cases_EUR=int(x[ix["n_cases_EUR"]]),
                    n_controls_EUR=int(x[ix["n_controls_EUR"]]), qc=x[ix["phenotype_qc_EUR"]], url=to_url(x[ix["aws_path"]]), url_tabix=to_url(x[ix["aws_path_tabix"]])))
res = Parallel(n_jobs=len(sel), backend="threading")(delayed(fetch)(s["tag"], s["url"], s["url_tabix"]) for s in sel)
print(res, flush=True)
pd.DataFrame(sel).drop(columns=["url_tabix"]).to_csv(f"{O}/3p14_panukb_phenotypes_extra.tsv", sep="\t", index=False)


import math
for s_ in sel:
    d = pd.read_csv(f"{D}/region_{s_['tag']}.tsv", sep="\t", dtype={"chr": str})
    for pos, nm in ((64699445, "rs4611812"), (64705365, "rs6795735")):
        x = d[d.pos == pos]
        for _, r in x.iterrows():
            print(s_["tag"], s_["n_cases_EUR"], nm, r.ref, r.alt, r.beta_EUR, r.se_EUR, r.beta_EUR / r.se_EUR, flush=True)

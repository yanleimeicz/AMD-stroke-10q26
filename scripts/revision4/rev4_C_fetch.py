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
PHENOS = [("AMD_phe362.29", "phecode", "362.29", ""), ("T2D_phe250.2", "phecode", "250.2", ""), ("T2D_E11", "icd10", "E11", ""),
          ("IS_I63", "icd10", "I63", ""), ("CAD_phe411.4", "phecode", "411.4", "")]

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
pd.DataFrame(sel).drop(columns=["url_tabix"]).to_csv(f"{O}/3p14_panukb_phenotypes.tsv", sep="\t", index=False)

LDB = "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD"; LDN = "chr3_64000001_67000001"
for ext in ("gz", "npz"):
    fn = f"{D}/ld/{LDN}.{ext}"
    if not os.path.exists(fn):
        with SES.get(f"{LDB}/{LDN}.{ext}", stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(fn + ".part", "wb") as f:
                for c in r.iter_content(8 << 20): f.write(c)
        os.replace(fn + ".part", fn)
    print(fn, os.path.getsize(fn), flush=True)

# ---------------------------------------------------------------- SuSiE / coloc inputs on a common variant set
ldv = pd.read_csv(f"{D}/ld/{LDN}.gz", sep=r"\s+"); ldv.columns = [c.lower() for c in ldv.columns]
M = sparse.load_npz(f"{D}/ld/{LDN}.npz").tocsr()
posmap = {}
for i, (pos, a1, a2) in enumerate(zip(ldv.position, ldv.allele1, ldv.allele2)): posmap.setdefault(int(pos), []).append((i, a1, a2))
def ld_index(pos, ref, alt):
    for i, a1, a2 in posmap.get(int(pos), []):
        if a1 == ref and a2 == alt: return i, 1
        if a1 == alt and a2 == ref: return i, -1
    return None, 0
REG = {}
for s in sel:
    d = pd.read_csv(f"{D}/region_{s['tag']}.tsv", sep="\t", dtype={"chr": str})
    d = d[d.beta_EUR.notna() & d.se_EUR.notna()].copy()
    if "low_confidence_EUR" in d: d = d[d.low_confidence_EUR.astype(str) != "True"]
    d = d[d.af_controls_EUR.between(0.005, 0.995)]
    d["key"] = d.pos.astype(str) + ":" + d.ref + ":" + d.alt
    REG[s["tag"]] = d.set_index("key")
common = set.intersection(*[set(d.index) for d in REG.values()])
base = REG[sel[0]["tag"]].loc[sorted(common, key=lambda k: int(k.split(":")[0]))]
idx = [ld_index(p, r, a) for p, r, a in zip(base.pos, base.ref, base.alt)]
base = base.assign(ld_i=[i for i, _ in idx], ld_s=[s for _, s in idx])
base = base[base.ld_i.notna()]
base = base[~base.ld_i.duplicated()]
keys = base.index.tolist(); li = base.ld_i.astype(int).tolist(); sg = base.ld_s.values.astype(float)
A = M[li][:, li].toarray(); A = A + A.T; np.fill_diagonal(A, 1.0); A = A * np.outer(sg, sg)
A.astype("<f8").tofile(f"{I}/R.f64")
for s in sel:
    d = REG[s["tag"]].loc[keys]
    n = s["n_cases_EUR"] + s["n_controls_EUR"]
    pd.DataFrame(dict(key=keys, pos=d.pos.values, ref=d.ref.values, alt=d.alt.values, af=d.af_controls_EUR.values, beta=d.beta_EUR.values,
                      se=d.se_EUR.values, z=(d.beta_EUR / d.se_EUR).values, n=n, s=s["n_cases_EUR"] / n)).to_csv(f"{I}/{s['tag']}_z.tsv", sep="\t", index=False)
print("variants", len(keys), flush=True)
for k, pos in (("rs4611812", 64699445), ("rs6795735", 64705365), ("rs6445425", 64719091)):
    m = [kk for kk in keys if kk.startswith(f"{pos}:")]
    print(k, m, flush=True)

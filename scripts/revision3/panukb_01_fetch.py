"""Round-3: Pan-UK Biobank 10q26 region fetch (remote tabix via HTTP range) + UKB British LD (Weissbrod/PolyFun npz).
Outputs: derived/panukb/region_<tag>.tsv, derived/panukb/ld/chr10_123000001_126000001.{npz,gz}, results/revision3/panukb/panukb_phenotypes.tsv"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import gzip, io, os, struct, zlib, csv, sys, time
import requests
SES = requests.Session(); SES.trust_env = False   # host has a stale system proxy (127.0.0.1:7890); connect directly
from joblib import Parallel, delayed

P = PROJ_ROOT; D = f"{P}/derived/panukb"; O = f"{P}/results/revision3/panukb"
os.makedirs(D, exist_ok=True); os.makedirs(O, exist_ok=True); os.makedirs(f"{D}/ld", exist_ok=True)
BASE = "https://pan-ukb-us-east-1.s3.amazonaws.com"
CHR, START, END = "10", 123_900_000, 124_600_000           # GRCh37
PHENOS = [  # tag, trait_type, phenocode, coding, role
    ("AMD_phe362.29", "phecode", "362.29", "", "AMD (hospital ICD, phecode 362.29) - primary"),
    ("AMD_phe362.2", "phecode", "362.2", "", "AMD (phecode 362.2, broader)"),
    ("AMD_sr6148_5", "categorical", "6148", "5", "AMD (self-reported eye problem, touchscreen subset)"),
    ("IS_I63", "icd10", "I63", "", "Ischemic stroke (ICD-10 I63, hospital) - primary"),
    ("IS_phe433.2", "phecode", "433.2", "", "Occlusion of cerebral arteries (phecode 433.2)"),
    ("CI_phe433.3", "phecode", "433.3", "", "Cerebral ischemia (phecode 433.3)"),
    ("AS_sr6150_3", "categorical", "6150", "3", "Any stroke (self-reported doctor-diagnosed) - AS proxy"),
    ("AS_sr20002_1081", "categorical", "20002", "1081", "Any stroke (self-reported non-cancer illness)"),
    ("CVD_phe433", "phecode", "433", "", "Cerebrovascular disease (phecode 433, broad)"),
    ("TIA_G45", "icd10", "G45", "", "TIA (ICD-10 G45)"),
    ("ICH_phe430.2", "phecode", "430.2", "", "Intracerebral hemorrhage (phecode 430.2)"),
    ("CAD_phe411.4", "phecode", "411.4", "", "Coronary atherosclerosis (phecode 411.4) - comparison"),
    ("IHD_phe411", "phecode", "411", "", "Ischemic heart disease (phecode 411) - comparison"),
]

def get(url, rng=None, tries=5):
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
            w = beg >> 14
            return next(o for o in ioff[w:] if o > 0)
    raise KeyError(chrom)

def bgzf_iter(data):
    while data:
        d = zlib.decompressobj(31); out = d.decompress(data); yield out; data = d.unused_data

def fetch(tag, url, tbi_url):
    voff = tbi_offsets(get(tbi_url), CHR, START)
    coff, uoff = voff >> 16, voff & 0xFFFF
    header = next(bgzf_iter(get(url, (0, 200_000)))).split(b"\n")[0].decode()
    buf, chunk, pos_done, lines = b"", 4_000_000, False, []
    start = coff
    text = b""
    while not pos_done:
        data = get(url, (start, start + chunk - 1))
        # keep only complete bgzf blocks
        blocks, used = [], 0
        while used < len(data):
            if len(data) - used < 18: break
            bsize = struct.unpack("<H", data[used + 16:used + 18])[0] + 1
            if used + bsize > len(data): break
            blocks.append(zlib.decompress(data[used:used + bsize], 31)); used += bsize
        if start == coff and blocks: blocks[0] = blocks[0][uoff:]
        text += b"".join(blocks); start += used
        rows = text.split(b"\n"); text = rows[-1]
        for r in rows[:-1]:
            x = r.split(b"\t", 2)
            if x[0].decode() not in (CHR, "chr" + CHR): pos_done = True; break
            pos = int(x[1])
            if pos > END: pos_done = True; break
            if pos >= START: lines.append(r.decode())
        if used == 0: break
    with open(f"{D}/region_{tag}.tsv", "w") as f:
        f.write(header + "\n" + "\n".join(lines) + "\n")
    return tag, len(lines), header

# manifest
man = list(csv.reader(gzip.open(f"{D}/phenotype_manifest.tsv.bgz", "rt"), delimiter="\t")); h = man[0]; ix = {k: i for i, k in enumerate(h)}
sel = []
for tag, tt, pc, cd, role in PHENOS:
    x = next(r for r in man[1:] if r[ix["trait_type"]] == tt and r[ix["phenocode"]] == pc and r[ix["coding"]] == cd and r[ix["pheno_sex"]] == "both_sexes")
    to_url = lambda s: BASE + "/" + s.split("pan-ukb-us-east-1/", 1)[1]
    sel.append(dict(tag=tag, role=role, trait_type=tt, phenocode=pc, coding=cd, description=x[ix["description"]], coding_description=x[ix["coding_description"]],
                    n_cases_EUR=x[ix["n_cases_EUR"]], n_controls_EUR=x[ix["n_controls_EUR"]], phenotype_qc_EUR=x[ix["phenotype_qc_EUR"]],
                    lambda_gc_EUR=x[ix["lambda_gc_EUR"]], url=to_url(x[ix["aws_path"]]), url_tabix=to_url(x[ix["aws_path_tabix"]])))
res = Parallel(n_jobs=len(sel), backend="threading")(delayed(fetch)(s["tag"], s["url"], s["url_tabix"]) for s in sel)
for s, (tag, n, hdr) in zip(sel, res):
    s["n_region_variants"] = n; s["header"] = hdr; print(tag, n, flush=True)
with open(f"{O}/panukb_phenotypes.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(sel[0]), delimiter="\t"); w.writeheader(); w.writerows(sel)

# UKB British-ancestry LD (Weissbrod et al. 2020, PolyFun), public bucket
LDB = "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD"
for ext in ("gz", "npz"):
    fn = f"{D}/ld/chr10_123000001_126000001.{ext}"
    if not os.path.exists(fn):
        with SES.get(f"{LDB}/chr10_123000001_126000001.{ext}", stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(fn + ".part", "wb") as f:
                for c in r.iter_content(8 << 20): f.write(c)
        os.replace(fn + ".part", fn)
    print(fn, os.path.getsize(fn), flush=True)

"""Round-3: ocular eQTL look-ups at 10q26 and 3p14 (Voigt 2026 sc-eQTL SI Table 3; Orozco 2020 eye-eqtl.com).
Builds results/revision3/eye_eqtl/{sc_eqtl_10q26_3p14.tsv, bulk_eqtl_10q26.tsv, rpe_htra1_bound_conditional.tsv}."""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import csv, json, re, html, glob, math, os
from statistics import NormalDist
P = PROJ_ROOT
D = f'{P}/derived/eye_eqtl'; O = f'{P}/results/revision3/eye_eqtl'; os.makedirs(O, exist_ok=True)
N = NormalDist()
# ---- 1000G EUR signed r (A1 = minor allele = AMD-risk allele for rs3750848-G, rs2672592-T, rs72631113-T)
R = {}
for f in ('ld_eur_rsigned.ld', 'ld2_eur_r.ld'):
    for i, l in enumerate(open(f'{D}/{f}')):
        if i == 0: continue
        x = l.split(); R[(x[2], x[5])] = R[(x[5], x[2])] = float(x[6])
A1 = {}
for l in open(f'{D}/ld_eur_frq.frq').readlines()[1:]:
    x = l.split(); A1[x[1]] = (x[2], x[3], float(x[4]))
def r(a, b): return 1.0 if a == b else R.get((a, b))
bim = {}
for l in open(f'{P}/data/raw/ref/g1000_eur/g1000_eur.bim'):
    c, rs, _, p, a1, a2 = l.split()
    if (c == '10' and 123_600_000 < int(p) < 124_400_000) or (c == '3' and 63_900_000 < int(p) < 65_500_000): bim[(c, int(p))] = rs
OFF = {'10': 1759516}   # hg38 -> hg19 offset verified at rs10490924/rs3750846/rs2672592/rs72631113
# ---- Voigt 2026 SI Table 3 (all eSNPs with FDR < 0.05)
rows = []
for x in csv.DictReader(open(f'{D}/web/v2_SI_Table3_media-1.csv')):
    c, p, ref, alt = x['snps'].split(':')[:4]; c = c.replace('chr', ''); p = int(p)
    keep = (c == '10' and 121_500_000 <= p <= 123_000_000 and x['gene'] in ('HTRA1', 'ARMS2', 'PLEKHA1')) or \
           (c == '3' and x['gene'] in ('ADAMTS9', 'ADAMTS9-AS1', 'ADAMTS9-AS2'))
    if not keep: continue
    rs = bim.get((c, p + OFF[c])) if c in OFF else ''
    n = sum(int(x[k]) for k in ('count_allele_1_1', 'count_allele_1_2', 'count_allele_2_2'))
    amd_risk = {'rs10490924': 'T', 'rs3750848': 'G', 'rs3750846': 'C', 'rs36212733': 'C', 'rs11200638': 'A'}.get(rs or '', '')
    rows.append(dict(source='Voigt2026_scEQTL_SI_Table3', cell_type=x['celltype'], tissue=x['tissue'], gene=x['gene'], variant_hg38=x['snps'],
                     rsid=rs or '', effect_allele=alt, amd_risk_allele=amd_risk,
                     beta_per_effect_allele=float(x['beta']), t=float(x['statistic']), p=float(x['pvalue']), fdr=float(x['FDR']),
                     maf=re.search(r'MAF=([\d.]+)', x['INFO']).group(1), snp_class=x['snp_class'], n_donors=n,
                     r_rs3750848=(r(rs, 'rs3750848') if rs else None) or '', r_rs2672592=(r(rs, 'rs2672592') if rs else None) or ''))
# explicit absence rows for the key variants
for rs, pos in [('rs2672592', 'chr10:122471234'), ('rs72631113', 'chr10:122453933')]:
    rows.append(dict(source='Voigt2026_scEQTL_SI_Table3', cell_type='all 30 cell types', tissue='Retina/RPE_Choroid', gene='HTRA1/ARMS2/PLEKHA1 (any gene)',
                     variant_hg38=pos, rsid=rs, effect_allele='T', amd_risk_allele='T', beta_per_effect_allele='', t='', p='', fdr='>=0.05 (not reported)',
                     maf='', snp_class='', n_donors='', r_rs3750848=r(rs, 'rs3750848'), r_rs2672592=r(rs, 'rs2672592')))
for rs, pos in [('rs4611812', 'chr3:64713769'), ('rs6795735', 'chr3:64719689'), ('3p14 13-variant aortic/AMD/T2D credible set', 'chr3:64713769-64744445')]:
    rows.append(dict(source='Voigt2026_scEQTL_SI_Table3', cell_type='all 30 cell types', tissue='Retina/RPE_Choroid', gene='any gene',
                     variant_hg38=pos, rsid=rs, effect_allele='', amd_risk_allele='', beta_per_effect_allele='', t='', p='', fdr='>=0.05 (not reported)',
                     maf='', snp_class='', n_donors='', r_rs3750848='', r_rs2672592=''))
keys = list(rows[0].keys())
with open(f'{O}/sc_eqtl_10q26_3p14.tsv', 'w') as f:
    f.write('\t'.join(keys) + '\n')
    for x in rows: f.write('\t'.join(str(x[k]) for k in keys) + '\n')
# ---- bounded conditional analysis for RPE HTRA1 (lead set = AMD haplotype; rs2672592 absent => |t| < smallest listed |t| for HTRA1-RPE)
rpe = [x for x in rows if x['gene'] == 'HTRA1' and x['cell_type'] == 'RPE']
t_lead = min(x['t'] for x in rpe)                    # -5.19 per AMD-risk allele (alt allele = risk allele)
t_floor = min(abs(x['t']) for x in rpe)              # 5.06: any variant with |t| above this would have been listed (per-gene FDR monotone)
out = []
for rr_ in (r('rs3750848', 'rs2672592'), 0.6187, 0.664):   # 1000G EUR (plink), revision-2 EUR value, FinnGen SISu
    z1 = t_lead
    for z2, lab in [(z1 * rr_, 'LD-only expectation'), (-t_floor, 'max admissible (same direction)'), (0.0, 'null marginal'), (t_floor, 'max admissible (opposite)')]:
        z1c = (z1 - rr_ * z2) / math.sqrt(1 - rr_ ** 2); z2c = (z2 - rr_ * z1) / math.sqrt(1 - rr_ ** 2)
        out.append(dict(r_lead_rs2672592=round(rr_, 4), scenario=lab, t_lead_marginal=round(z1, 3), t_rs2672592_marginal=round(z2, 3),
                        z_lead_given_rs2672592=round(z1c, 3), p_lead_given_rs2672592='%.3g' % (2 * N.cdf(-abs(z1c))),
                        z_rs2672592_given_lead=round(z2c, 3), p_rs2672592_given_lead='%.3g' % (2 * N.cdf(-abs(z2c)))))
with open(f'{O}/rpe_htra1_bound_conditional.tsv', 'w') as f:
    k = list(out[0].keys()); f.write('\t'.join(k) + '\n')
    for x in out: f.write('\t'.join(str(x[i]) for i in k) + '\n')
# ---- Orozco 2020 bulk eQTL (eye-eqtl.com, FDR <= 0.05 only)
bulk = []
for fn in sorted(glob.glob(f'{D}/web/eyeeqtl/*.html')):
    s = open(fn, encoding='utf-8', errors='ignore').read()
    m = re.search(r'id="div-data"[^>]*data="([^"]*)"', s) or re.search(r'data="([^"]*)"[^>]*id="div-data"', s)
    q = os.path.basename(fn)[:-5]
    d = json.loads(html.unescape(m.group(1))) if m else []
    if not d:
        bulk.append(dict(source='Orozco2020_eye-eqtl.com', query=q, tissue='all 4 (retina/RPE x macula/non-macula)', gene='', snp='', p='', fdr='none at FDR<=0.05',
                         beta='', ref='', alt='', minor='', maf='', snp_bp_hg38='', r_rs3750848='', r_rs2672592='', r_rs4611812=''))
        continue
    for x in sorted(d, key=lambda x: x['pvalue']):
        s_ = x['snps']
        bulk.append(dict(source='Orozco2020_eye-eqtl.com', query=q, tissue=x['tissue'], gene=x['symbol'], snp=s_, p='%.3g' % x['pvalue'], fdr='%.3g' % x['fdr'],
                         beta=round(x['beta'], 3), ref=x['ref'], alt=x['alt'], minor=x['minor'], maf=round(x['maf'], 3), snp_bp_hg38=x['snp_bp'],
                         r_rs3750848=r(s_, 'rs3750848') or '', r_rs2672592=r(s_, 'rs2672592') or '', r_rs4611812=r(s_, 'rs4611812') or ''))
with open(f'{O}/bulk_eqtl_10q26.tsv', 'w') as f:
    k = list(bulk[0].keys()); f.write('\t'.join(k) + '\n')
    for x in bulk: f.write('\t'.join(str(x[i]) for i in k) + '\n')
print(len(rows), len(out), len(bulk)); print(t_lead, t_floor)
for x in out: print(x)

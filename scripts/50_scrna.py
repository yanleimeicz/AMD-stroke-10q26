"""Single-cell localisation of candidate shared genes (CELLxGENE h5ad: retina/RPE/choroid and brain vasculature).

Usage: python 50_scrna.py <gene_list.txt>
Outputs -> results/scrna/<dataset>/ : UMAP by cell type, dot plot, violin plots, per-cell-type mean/fraction table,
and (if both AMD and normal donors exist) donor-level pseudo-bulk comparison per cell type.
"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import sys, os, warnings, faulthandler
faulthandler.enable()
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
warnings.filterwarnings("ignore")
sc.settings.n_jobs = 30
RAW = Path(f"{PROJ_ROOT}/data/raw/sc"); OUTR = Path(f"{PROJ_ROOT}/results/scrna")
genes_req = [g.strip() for g in open(sys.argv[1]) if g.strip()]
DATASETS = {
    "rpe_choroid": "rpe_choroid_scrna_all.h5ad",
    "retina_rpe_cowan": "retina_rpe_periphery_cowan.h5ad",
    "brain_vascular_siletti": "brain_vascular_siletti.h5ad",
    "brain_vascular_multiomics": "brain_vascular_multiomics_snrna.h5ad",
}


def lognorm(ad):
    """CELLxGENE: X may already be log-normalised; raw counts live in .raw.X. Return log1p(CP10k)."""
    X = ad.raw.X if ad.raw is not None else ad.X
    Xm = X[:2000].toarray() if sp.issparse(X) else np.asarray(X[:2000])
    is_counts = np.allclose(Xm, np.round(Xm)) and Xm.max() > 20
    if is_counts:
        a = sc.AnnData(X=X.copy(), obs=ad.obs.copy(), var=(ad.raw.var if ad.raw is not None else ad.var).copy())
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    else:
        a = sc.AnnData(X=ad.X.copy(), obs=ad.obs.copy(), var=ad.var.copy())
    for k in ad.obsm.keys():                                   # normalise embedding key names (e.g. X_UMAP -> X_umap)
        a.obsm["X_umap" if k.lower() == "x_umap" else k] = np.asarray(ad.obsm[k])
    sym = a.var["feature_name"].astype(str) if "feature_name" in a.var else a.var_names.astype(str)
    a.var_names = pd.Index(sym); a.var_names_make_unique()
    return a


ONLY = os.environ.get("SC_ONLY")
step = lambda m: print("   ..", m, flush=True)
for name, fn in DATASETS.items():
    if ONLY and name != ONLY: continue
    f = RAW / fn
    if not f.exists():
        print("missing", name); continue
    out = OUTR / name; out.mkdir(parents=True, exist_ok=True)
    print(f"== {name}", flush=True)
    ad = sc.read_h5ad(f)
    a = lognorm(ad); del ad
    ct = "cell_type" if "cell_type" in a.obs else a.obs.columns[0]
    a.obs[ct] = a.obs[ct].astype(str)
    keep_ct = a.obs[ct].value_counts(); keep_ct = keep_ct[keep_ct >= 30].index
    a = a[a.obs[ct].isin(keep_ct)].copy()
    genes = [g for g in genes_req if g in a.var_names]
    print(f"   cells {a.n_obs}, cell types {a.obs[ct].nunique()}, genes found {len(genes)}/{len(genes_req)}", flush=True)
    if "X_umap" not in a.obsm:
        sub = a if a.n_obs <= 60000 else a[np.random.default_rng(1).choice(a.n_obs, 60000, replace=False)].copy()
        sc.pp.highly_variable_genes(sub, n_top_genes=2000); sc.pp.pca(sub, n_comps=30, use_highly_variable=True)
        sc.pp.neighbors(sub); sc.tl.umap(sub); a = sub
    sc.settings.figdir = str(out); sc.settings.set_figure_params(dpi_save=300, vector_friendly=True)
    # raster (PNG) for cell-level scatter plots; subsample very large datasets for plotting only
    step("subsample")
    ap = a if a.n_obs <= 60000 else a[np.random.default_rng(7).choice(a.n_obs, 60000, replace=False)].copy()
    step("umap celltype")
    sc.pl.umap(ap, color=ct, frameon=False, legend_loc="on data", legend_fontsize=5, show=False, save=f"_{name}_celltype.png")
    sc.pl.umap(ap, color=ct, frameon=False, legend_fontsize=6, show=False, save=f"_{name}_celltype_legend.png")
    if genes:
        for i in range(0, len(genes), 12):
            sc.pl.umap(ap, color=genes[i:i + 12], frameon=False, cmap="Reds", ncols=4, show=False, save=f"_{name}_genes_{i // 12 + 1}.png")
        step("dotplot")
        sc.pl.dotplot(a, genes, groupby=ct, standard_scale="var", show=False, save=f"{name}_dotplot.pdf")
        step("violin")
        sc.pl.stacked_violin(ap, genes[:20], groupby=ct, swap_axes=False, show=False, save=f"{name}_violin.pdf")
        step("tables")
        # per cell-type mean expression and detection fraction
        Xg = a[:, genes].X; Xg = Xg.toarray() if sp.issparse(Xg) else np.asarray(Xg)
        df = pd.DataFrame(Xg, columns=genes); df["ct"] = a.obs[ct].values
        mean = df.groupby("ct").mean(); frac = (df.set_index("ct") > 0).groupby(level=0).mean()
        mean.to_csv(out / "celltype_mean_expr.tsv", sep="\t"); frac.to_csv(out / "celltype_fraction_expr.tsv", sep="\t")
        top = pd.DataFrame({"top_celltype": mean.idxmax(), "top_mean": mean.max(), "frac_in_top": [frac.loc[mean[g].idxmax(), g] for g in genes]})
        top.to_csv(out / "gene_top_celltype.tsv", sep="\t"); print(top.to_string(), flush=True)
        # donor-level pseudo-bulk AMD vs normal per cell type (if available)
        dis = a.obs.get("disease"); donor = a.obs.get("donor_id")
        if dis is not None and donor is not None:
            dis = dis.astype(str); amd = dis.str.contains("macular degeneration", case=False)
            if amd.sum() > 0 and (dis == "normal").sum() > 0:
                df["donor"] = donor.astype(str).values; df["grp"] = np.where(amd, "AMD", np.where(dis == "normal", "normal", "other"))
                pb = df[df.grp != "other"].groupby(["ct", "donor", "grp"]).mean().reset_index()
                rows = []
                for (c, g) in [(c, g) for c in pb.ct.unique() for g in genes]:
                    x = pb[(pb.ct == c) & (pb.grp == "AMD")][g]; y = pb[(pb.ct == c) & (pb.grp == "normal")][g]
                    if len(x) >= 2 and len(y) >= 2:
                        rows.append(dict(cell_type=c, gene=g, n_amd=len(x), n_normal=len(y), mean_amd=x.mean(), mean_normal=y.mean(),
                                         log2fc=np.log2((x.mean() + 1e-3) / (y.mean() + 1e-3)), p_wilcox=stats.mannwhitneyu(x, y).pvalue))
                if rows:
                    r = pd.DataFrame(rows); r["fdr"] = stats.false_discovery_control(r.p_wilcox)
                    r.sort_values("p_wilcox").to_csv(out / "pseudobulk_AMD_vs_normal.tsv", sep="\t", index=False)
                    print("   pseudo-bulk AMD vs normal:", (r.p_wilcox < 0.05).sum(), "nominal hits", flush=True)
    del a
print("scRNA done")

"""Single-cell localisation of ADAMTS9 / ADAMTS9-AS2 (3p14 deep-dive) in retina / RPE-choroid and brain vasculature
(CELLxGENE h5ad, same files and log1p(CP10k) normalisation as 50_scrna.py; library size from all genes, then gene subset).

Outputs -> results/adamts9/scrna/: per-cell-type mean / fraction table, dot plot and violin plots per dataset, and
donor-level pseudo-bulk AMD vs normal comparison per cell type when the dataset has both AMD and normal donors.
"""
import os as _os, pathlib as _pl   # ---- project root: set COMORBIDITY_ROOT to override; default = parent of the 'scripts' folder
PROJ_ROOT = (_os.environ.get("COMORBIDITY_ROOT") or next((str(p.parent) for p in _pl.Path(__file__).resolve().parents if p.name == "scripts"), str(_pl.Path.cwd()))).replace("\\", "/")
import warnings
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
warnings.filterwarnings("ignore")
RAW = Path(f"{PROJ_ROOT}/data/raw/sc"); OUT = Path(f"{PROJ_ROOT}/results/adamts9/scrna"); OUT.mkdir(parents=True, exist_ok=True)
GENES = ["ADAMTS9", "ADAMTS9-AS2"]
DATASETS = {
    "rpe_choroid": "rpe_choroid_scrna_all.h5ad",
    "retina_rpe_cowan": "retina_rpe_periphery_cowan.h5ad",
    "brain_vascular_siletti": "brain_vascular_siletti.h5ad",
    "brain_vascular_multiomics": "brain_vascular_multiomics_snrna.h5ad",
}


def gene_expr(a):
    """log1p(CP10k) of GENES; library size from all genes of raw counts (or X if already normalised)."""
    use_raw = a.raw is not None
    X = a.raw.X if use_raw else a.X
    var = a.raw.var if use_raw else a.var
    sym = var["feature_name"].astype(str).values if "feature_name" in var else np.asarray(var.index.astype(str))
    Xs = X[:2000]; Xs = Xs.toarray() if sp.issparse(Xs) else np.asarray(Xs)
    is_counts = np.allclose(Xs, np.round(Xs)) and Xs.max() > 20
    idx = {g: np.where(sym == g)[0] for g in GENES}
    cols = [i[0] for g, i in idx.items() if len(i)]
    found = [g for g, i in idx.items() if len(i)]
    sub = X[:, cols]; sub = sub.toarray() if sp.issparse(sub) else np.asarray(sub)
    if is_counts:
        tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1
        sub = np.log1p(sub / tot[:, None] * 1e4)
    return pd.DataFrame(sub, columns=found, index=a.obs_names), is_counts


def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out


summary = []
for name, fn in DATASETS.items():
    f = RAW / fn
    if not f.exists():
        print("missing", name); continue
    print(f"== {name}", flush=True)
    a = ad.read_h5ad(f)
    df, is_counts = gene_expr(a)
    ct = "cell_type" if "cell_type" in a.obs else a.obs.columns[0]
    df["cell_type"] = a.obs[ct].astype(str).values
    for c in ("disease", "donor_id", "tissue"):
        if c in a.obs: df[c] = a.obs[c].astype(str).values
    keep = df.cell_type.value_counts(); df = df[df.cell_type.isin(keep[keep >= 30].index)]
    genes = [g for g in GENES if g in df]
    print(f"   cells {len(df)}, cell types {df.cell_type.nunique()}, genes {genes}, counts={is_counts}", flush=True)
    if not genes: continue
    tab = df.groupby("cell_type")[genes].agg(["mean", lambda x: (x > 0).mean()])
    tab.columns = [f"{g}_{'mean' if s == 'mean' else 'frac'}" for g, s in tab.columns]
    tab["n_cells"] = df.cell_type.value_counts()
    tab = tab.sort_values(f"{genes[0]}_mean", ascending=False)
    tab.to_csv(OUT / f"{name}_celltype_expr.tsv", sep="\t")
    for g in genes:
        top = tab.iloc[0]
        summary.append(dict(dataset=name, gene=g, n_cells=len(df), top_cell_type=tab[f"{g}_mean"].idxmax(),
                            top_mean=tab[f"{g}_mean"].max(), top_frac=tab.loc[tab[f"{g}_mean"].idxmax(), f"{g}_frac"],
                            overall_frac=(df[g] > 0).mean(),
                            top3=", ".join(f"{c} ({tab.loc[c, f'{g}_mean']:.2f}/{tab.loc[c, f'{g}_frac']:.0%})"
                                           for c in tab[f"{g}_mean"].sort_values(ascending=False).index[:3])))
    # dot plot (all cell types with >= 30 cells)
    order = tab.index[::-1]
    fig, ax = plt.subplots(figsize=(2.2 + 1.3 * len(genes), max(3, 0.24 * len(order))))
    for j, g in enumerate(genes):
        m = tab.loc[order, f"{g}_mean"].values; fr = tab.loc[order, f"{g}_frac"].values
        sc_ = ax.scatter(np.full(len(order), j), np.arange(len(order)), s=fr * 260 + 2, c=m, cmap="Reds", vmin=0,
                         vmax=max(tab[[f"{x}_mean" for x in genes]].values.max(), 0.1), edgecolor="grey", linewidth=0.3)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=7)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes); ax.set_xlim(-0.6, len(genes) - 0.4)
    plt.colorbar(sc_, ax=ax, label="mean log1p(CP10k)", shrink=0.5)
    for fr_ in (0.1, 0.4, 0.8): ax.scatter([], [], s=fr_ * 260 + 2, c="grey", label=f"{fr_:.0%}")
    ax.legend(title="fraction", loc="lower right", bbox_to_anchor=(1.9, 0), fontsize=7, frameon=False)
    ax.set_title(name, fontsize=9); plt.tight_layout()
    plt.savefig(OUT / f"{name}_dotplot.pdf"); plt.savefig(OUT / f"{name}_dotplot.png", dpi=200); plt.close()
    # violin of ADAMTS9 in the top 12 cell types
    g = genes[0]; top12 = tab[f"{g}_mean"].sort_values(ascending=False).index[:12]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.violinplot([df.loc[df.cell_type == c, g].values for c in top12], showmedians=True, widths=0.8)
    ax.set_xticks(range(1, len(top12) + 1)); ax.set_xticklabels(top12, rotation=50, ha="right", fontsize=7)
    ax.set_ylabel(f"{g} log1p(CP10k)"); ax.set_title(name, fontsize=9); plt.tight_layout()
    plt.savefig(OUT / f"{name}_violin_{g}.pdf"); plt.savefig(OUT / f"{name}_violin_{g}.png", dpi=200); plt.close()
    # donor pseudo-bulk AMD vs normal
    if "disease" in df and "donor_id" in df:
        dis = df.disease.str.lower()
        grp = np.where(dis.str.contains("macular degeneration"), "AMD", np.where(dis == "normal", "normal", "other"))
        print("   disease groups:", pd.Series(grp).value_counts().to_dict(), flush=True)
        if (grp == "AMD").any() and (grp == "normal").any():
            d2 = df.assign(grp=grp)[lambda x: x.grp != "other"]
            pb = d2.groupby(["cell_type", "donor_id", "grp"])[genes].mean().reset_index()
            ncell = d2.groupby(["cell_type", "donor_id"]).size().rename("n").reset_index()
            pb = pb.merge(ncell, on=["cell_type", "donor_id"]); pb = pb[pb.n >= 10]
            rows = []
            for (c, gg) in [(c, gg) for c in pb.cell_type.unique() for gg in genes]:
                x = pb[(pb.cell_type == c) & (pb.grp == "AMD")][gg]; y = pb[(pb.cell_type == c) & (pb.grp == "normal")][gg]
                if len(x) >= 3 and len(y) >= 3:
                    rows.append(dict(cell_type=c, gene=gg, n_amd=len(x), n_normal=len(y), mean_amd=x.mean(), mean_normal=y.mean(),
                                     log2fc=np.log2((x.mean() + 1e-3) / (y.mean() + 1e-3)),
                                     p_wilcox=stats.mannwhitneyu(x, y, alternative="two-sided").pvalue))
            if rows:
                r = pd.DataFrame(rows); r["fdr"] = r.groupby("gene").p_wilcox.transform(bh)
                r.sort_values("p_wilcox").to_csv(OUT / f"{name}_pseudobulk_AMD_vs_normal.tsv", sep="\t", index=False)
                print("   pseudo-bulk tests:", len(r), "; P<0.05:", int((r.p_wilcox < 0.05).sum()), flush=True)
    del a, df
pd.DataFrame(summary).to_csv(OUT / "summary.tsv", sep="\t", index=False)
print("scrna adamts9 done", flush=True)

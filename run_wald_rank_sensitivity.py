import numpy as np
import pandas as pd
import decoupler as dc
import anndata as ad
from scipy.stats import spearmanr
from pathlib import Path

GENE_COL = "gene_symbol"  
STAT_COL = "stat"         
NES_THRESHOLD = 1.5
BASE = Path(__file__).resolve().parent
DATA = BASE / "geo_data"
OUT = BASE / "results"

FILES = {
    "spaceflight": DATA / "spaceflight_deg_full.csv",
    "sarcopenia": DATA / "sarcopenia_deg.csv",
    "disuse": DATA / "disuse_deg.csv",
}



def load_wald_matrix():
    """Load and merge per-gene Wald statistics for all three contrasts."""
    frames = []
    for cond, path in FILES.items():
        df = pd.read_csv(path)
        sub = df[[GENE_COL, STAT_COL]].rename(columns={STAT_COL: f"wald_{cond}"})
        sub = sub.dropna(subset=[GENE_COL, f"wald_{cond}"])
        sub["abs_stat"] = sub[f"wald_{cond}"].abs()
        sub = sub.sort_values("abs_stat", ascending=False).drop_duplicates(subset=GENE_COL)
        sub = sub.drop(columns="abs_stat")
        frames.append(sub.set_index(GENE_COL))
    wald = pd.concat(frames, axis=1)
    return wald


def rank_transform(wald_df):
    """Within-dataset rank-z transform, scale-free alternative to raw Wald."""
    out = wald_df.copy()
    for col in wald_df.columns:
        ranks = wald_df[col].rank(method="average")
        n = ranks.notna().sum()
        out[col] = (ranks - (n + 1) / 2) / np.sqrt((n**2 - 1) / 12)
    return out


def run_ulm(signature_df, net):
    """Score a gene x condition signature matrix with decoupleR ULM."""
    adata = ad.AnnData(X=signature_df.T.values,
                        obs=pd.DataFrame(index=signature_df.columns),
                        var=pd.DataFrame(index=signature_df.index))
    dc.mt.ulm(adata, net=net, tmin=5)
    nes = adata.obsm["score_ulm"].T
    return nes


def classify_pan_atrophy(nes_df, threshold=NES_THRESHOLD):
    active = (nes_df.abs() >= threshold)
    n_active = active.sum(axis=1)
    return n_active[n_active == nes_df.shape[1]].index.tolist()


def compare_rankings(nes_primary, nes_rank_based, top_n=47):
    combined_primary = nes_primary.abs().sum(axis=1).rank(ascending=False, method="min")
    combined_rankbased = nes_rank_based.abs().sum(axis=1).rank(ascending=False, method="min")
    common = combined_primary.index.intersection(combined_rankbased.index)
    rho, p = spearmanr(combined_primary[common], combined_rankbased[common])

    pan_primary = set(classify_pan_atrophy(nes_primary))
    pan_rankbased = set(classify_pan_atrophy(nes_rank_based))
    jaccard = len(pan_primary & pan_rankbased) / len(pan_primary | pan_rankbased)

    top_primary = set(combined_primary.nsmallest(top_n).index)
    top_rankbased = set(combined_rankbased.nsmallest(top_n).index)
    top_jaccard = len(top_primary & top_rankbased) / len(top_primary | top_rankbased)

    return {
        "spearman_rho_full_ranking": rho,
        "spearman_p": p,
        "pan_atrophy_jaccard": jaccard,
        "n_pan_atrophy_primary": len(pan_primary),
        "n_pan_atrophy_rankbased": len(pan_rankbased),
        f"top{top_n}_rank_jaccard": top_jaccard,
    }


if __name__ == "__main__":
    net = dc.op.collectri(organism="mouse")

    wald = load_wald_matrix()
    wald_imputed = wald.fillna(0.0)

    nes_primary = run_ulm(wald_imputed, net)

    wald_rankz = rank_transform(wald_imputed)
    nes_rankbased = run_ulm(wald_rankz, net)

    results = compare_rankings(nes_primary, nes_rankbased)
    OUT.mkdir(exist_ok=True)
    pd.Series(results).to_csv(OUT / "wald_rank_sensitivity_results.csv")
    print(results)

import os
import numpy as np
import pandas as pd
import decoupler as dc
import anndata as ad
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from scipy.stats import pearsonr

GENE_COL = "gene_symbol"
STAT_COL = "stat"
DATA_DIR = "geo_data"
N_PERM = 10000
N_BOOT = 10000
SEED = 42
NET_CACHE_PATH = "collectri_mouse_cache.csv"
ZENODO_URL = "https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1"
ZENODO_RAW_PATH = "CollecTRI_regulons_raw.csv"
ZENODO_TIMEOUT = (30, 600)

FILES = {
    "spaceflight": f"{DATA_DIR}/spaceflight_deg_full.csv",
    "sarcopenia":  f"{DATA_DIR}/sarcopenia_deg.csv",
    "disuse":      f"{DATA_DIR}/disuse_deg.csv",
}

PAIRS = [("spaceflight", "disuse"), ("sarcopenia", "disuse"), ("spaceflight", "sarcopenia")]

def try_robust_zenodo_download():
    print("Attempting robust direct download from Zenodo...")
    session = requests.Session()
    retries = Retry(total=5, backoff_factor=10, status_forcelist=[502, 503, 504],
                     allowed_methods=["GET"])
    session.mount("https://", HTTPAdapter(max_retries=retries))

    if not os.path.exists(ZENODO_RAW_PATH):
        with session.get(ZENODO_URL, stream=True, timeout=ZENODO_TIMEOUT) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            downloaded = 0
            with open(ZENODO_RAW_PATH, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total:
                            print(f"\r  {downloaded/1e6:.1f}/{total/1e6:.1f} MB", end="")
        print("\nDownload complete, parsing...")
    else:
        print(f"Reusing already-downloaded {ZENODO_RAW_PATH}")

    raw = pd.read_csv(ZENODO_RAW_PATH)
    cols_lower = {c.lower(): c for c in raw.columns}

    def find_col(cands):
        for c in cands:
            if c in cols_lower:
                return cols_lower[c]
        return None

    src_col = find_col(["source", "tf", "source_genesymbol"])
    tgt_col = find_col(["target", "target_genesymbol", "gene"])
    wt_col = find_col(["weight", "mor", "sign"])

    if not all([src_col, tgt_col, wt_col]):
        raise ValueError(f"Could not identify source/target/weight columns. "
                          f"Columns found: {list(raw.columns)}")

    net = raw[[src_col, tgt_col, wt_col]].rename(
        columns={src_col: "source", tgt_col: "target", wt_col: "weight"})

    net["source"] = net["source"].astype(str).str.capitalize()
    net["target"] = net["target"].astype(str).str.capitalize()
    net = net.drop_duplicates(subset=["source", "target"])
    return net, "Zenodo (direct, capitalize-approximated mouse orthology)"


def build_from_omnipath():
    print("Falling back to OmniPath REST API...")
    import omnipath as op

    def determine_weight(row):
        if row["consensus_stimulation"] and not row["consensus_inhibition"]:
            return 1
        elif row["consensus_inhibition"] and not row["consensus_stimulation"]:
            return -1
        elif row["is_stimulation"] and not row["is_inhibition"]:
            return 1
        elif row["is_inhibition"] and not row["is_stimulation"]:
            return -1
        else:
            return 1

    raw = op.interactions.CollecTRI.get(genesymbols=True, organism="mouse")
    net = raw[["source_genesymbol", "target_genesymbol", "is_stimulation",
               "is_inhibition", "consensus_stimulation", "consensus_inhibition"]].copy()
    net = net.rename(columns={"source_genesymbol": "source", "target_genesymbol": "target"})
    net["weight"] = net.apply(determine_weight, axis=1)
    net = net[["source", "target", "weight"]].drop_duplicates(subset=["source", "target"])
    return net, "OmniPath (fallback, smaller regulon, native mouse symbols)"


def get_collectri_net(force_refresh=False):
    if os.path.exists(NET_CACHE_PATH) and not force_refresh:
        print(f"Loading cached CollecTRI regulon from {NET_CACHE_PATH}")
        net = pd.read_csv(NET_CACHE_PATH)
        source = "local cache"
    else:
        try:
            net, source = try_robust_zenodo_download()
        except Exception as e:
            print(f"Zenodo path failed ({e}); trying OmniPath fallback...")
            net, source = build_from_omnipath()
        net.to_csv(NET_CACHE_PATH, index=False)

    print(f"Regulon source: {source}")
    print(f"Regulon size: {net['source'].nunique()} TFs, "
          f"{net['target'].nunique()} targets, {len(net)} edges")
    print("Manuscript reference size: 1,165 TFs / 16,883 targets / 43,226 edges")
    return net

def load_wald_matrix():
    frames = []
    for cond, path in FILES.items():
        df = pd.read_csv(path)
        sub = df[[GENE_COL, STAT_COL]].rename(columns={STAT_COL: f"wald_{cond}"})
        sub = sub.dropna(subset=[GENE_COL, f"wald_{cond}"])
        sub["abs_stat"] = sub[f"wald_{cond}"].abs()
        sub = sub.sort_values("abs_stat", ascending=False).drop_duplicates(subset=GENE_COL)
        sub = sub.drop(columns="abs_stat")
        frames.append(sub.set_index(GENE_COL))
    return pd.concat(frames, axis=1)


def run_ulm(signature_df, net):
    adata = ad.AnnData(X=signature_df.T.values,
                        obs=pd.DataFrame(index=signature_df.columns),
                        var=pd.DataFrame(index=signature_df.index))
    dc.mt.ulm(adata, net=net, tmin=5)
    return adata.obsm["score_ulm"].T


def permutation_bootstrap_test(nes_df, cond_a, cond_b, n_perm=N_PERM, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    x = nes_df[f"wald_{cond_a}"].to_numpy()
    y = nes_df[f"wald_{cond_b}"].to_numpy()
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    n = len(x)

    r_obs, analytic_p = pearsonr(x, y)

    perm_rs = np.empty(n_perm)
    for i in range(n_perm):
        perm_rs[i] = np.corrcoef(x, rng.permutation(y))[0, 1]
    perm_p = (np.sum(np.abs(perm_rs) >= np.abs(r_obs)) + 1) / (n_perm + 1)

    idx = np.arange(n)
    boot_rs = np.empty(n_boot)
    for i in range(n_boot):
        samp = rng.choice(idx, size=n, replace=True)
        boot_rs[i] = np.corrcoef(x[samp], y[samp])[0, 1]
    ci_lo, ci_hi = np.percentile(boot_rs, [2.5, 97.5])

    return {
        "pair": f"{cond_a}-{cond_b}",
        "n_TFs": n,
        "r_observed": r_obs,
        "analytic_p": analytic_p,
        "permutation_p": perm_p,
        "bootstrap_95CI_low": ci_lo,
        "bootstrap_95CI_high": ci_hi,
    }


if __name__ == "__main__":
    net = get_collectri_net(force_refresh=True)
    wald = load_wald_matrix().fillna(0.0)
    nes = run_ulm(wald, net)

    results = [permutation_bootstrap_test(nes, a, b) for a, b in PAIRS]
    out = pd.DataFrame(results)
    out.to_csv("correlation_permutation_results.csv", index=False)
    print()
    print(out.to_string(index=False))
    print()
    print("Compare r_observed above against manuscript values: "
          "0.457 (spaceflight-disuse), 0.313 (sarcopenia-disuse), "
          "0.146 (spaceflight-sarcopenia) before using these numbers.")

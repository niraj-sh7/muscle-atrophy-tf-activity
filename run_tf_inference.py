from pathlib import Path

import pandas as pd
import anndata as ad
import decoupler as dc

BASE = Path(__file__).resolve().parent
DATA = BASE / "geo_data"
OUT = BASE / "results"
OUT.mkdir(exist_ok=True)

paths = {
    "spaceflight": DATA / "spaceflight_deg_full.csv",
    "sarcopenia": DATA / "sarcopenia_deg.csv",
    "disuse": DATA / "disuse_deg.csv",
}
missing = [f"{name}: {path}" for name, path in paths.items() if not path.exists()]
if missing:
    raise FileNotFoundError(
        "Missing DE inputs required for TF inference:\n" + "\n".join(missing)
    )

sf = pd.read_csv(paths["spaceflight"])
sa = pd.read_csv(paths["sarcopenia"])
du = pd.read_csv(paths["disuse"])

print("Spaceflight n=", len(sf))
print("Sarcopenia  n=", len(sa))
print("Disuse      n=", len(du))


def collapse(df, label):
    """Collapse to gene_symbol -> Wald stat, taking max |stat| per symbol."""
    d = df[["gene_symbol", "stat"]].dropna()
    d = d.loc[d.groupby("gene_symbol")["stat"].apply(lambda s: s.abs().idxmax())]
    s = d.set_index("gene_symbol")["stat"]
    s.name = label
    return s


s_sf = collapse(sf, "Spaceflight")
s_sa = collapse(sa, "Sarcopenia")
s_du = collapse(du, "Disuse_Atrophy")

mat = pd.concat([s_sf, s_sa, s_du], axis=1)
print("Union signature matrix:", mat.shape)
mat = mat.fillna(0.0)
mat.to_csv(OUT / "signature_matrix_genes.csv")

net_path = BASE / "data" / "collectri_mouse.csv"
if not net_path.exists():
    raise FileNotFoundError(f"Missing regulon file: {net_path}")
net = pd.read_csv(net_path)
net["mor"] = net["weight"]
print("Regulon TFs:", net["source"].nunique(), "edges:", len(net))

adata = ad.AnnData(mat.T.astype(float))
adata.var_names = mat.index.astype(str)

dc.mt.ulm(adata, net=net, tmin=5, verbose=True)

print(list(adata.obsm.keys()))
es = adata.obsm["score_ulm"]
pv = adata.obsm["padj_ulm"]
es.index = adata.obs_names
pv.index = adata.obs_names
print("ES shape:", es.shape)

es_long = es.reset_index().melt(id_vars="index", var_name="TF", value_name="NES")
es_long = es_long.rename(columns={"index": "condition"})
pv_long = pv.reset_index().melt(id_vars="index", var_name="TF", value_name="padj")
pv_long = pv_long.rename(columns={"index": "condition"})
acts = es_long.merge(pv_long, on=["condition", "TF"])
acts.to_csv(OUT / "tf_activities_long.csv", index=False)

nes_wide = es.T
pv_wide = pv.T
nes_wide.to_csv(OUT / "tf_activity_NES_wide.csv")
pv_wide.to_csv(OUT / "tf_activity_padj_wide.csv")
print("Conditions:", nes_wide.columns.tolist())
print("Number of TFs:", nes_wide.shape[0])
print("Top by absolute NES (Spaceflight):")
print(nes_wide.reindex(nes_wide["Spaceflight"].abs().sort_values(ascending=False).index)
      .head(15).round(2).to_string())

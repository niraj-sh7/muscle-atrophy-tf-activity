from pathlib import Path

import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import mygene

BASE = Path(__file__).resolve().parent
DATA = BASE / "data"
OUT = BASE / "geo_data"
src = DATA / "GSE273092_featureCounts.txt"
if not src.exists():
    raise FileNotFoundError(
        f"Missing {src}. This script needs the raw GSE273092 featureCounts table."
    )

cnt = pd.read_csv(src, sep="\t", comment="#")
print("FC raw shape:", cnt.shape)
mat = cnt.set_index("Geneid").drop(columns=["Chr", "Start", "End", "Strand", "Length"])
mat.columns = [c.replace(".bam", "") for c in mat.columns]
mat = mat.astype(int)
print("Matrix:", mat.shape)
print("Samples:", mat.columns.tolist())

meta_dict = {}
for i in range(1, 19):
    sid = f"JR-JC-{i}"
    if i in (1, 2, 3):
        gt, intv = "flox", "ctrl"
    elif i in (4, 5, 6):
        gt, intv = "mKO", "ctrl"
    elif i in (7, 8, 9):
        gt, intv = "flox", "unloading"
    elif i in (10, 11, 12):
        gt, intv = "mKO", "unloading"
    elif i in (13, 14, 15):
        gt, intv = "flox", "reloading"
    else:
        gt, intv = "mKO", "reloading"
    meta_dict[sid] = (gt, intv)

meta = pd.DataFrame(
    {"sample": list(meta_dict.keys()),
     "genotype": [v[0] for v in meta_dict.values()],
     "intervention": [v[1] for v in meta_dict.values()]}
).set_index("sample")

keep = meta[(meta["genotype"] == "flox") &
            (meta["intervention"].isin(["ctrl", "unloading"]))].index.tolist()
mat_sub = mat[keep]
meta_sub = meta.loc[keep].copy()
meta_sub["group"] = meta_sub["intervention"].map(
    {"ctrl": "control", "unloading": "unloading"})
print(meta_sub)

mat_sub = mat_sub[mat_sub.sum(axis=1) >= 10]
print("Filtered matrix:", mat_sub.shape)

dds = DeseqDataSet(counts=mat_sub.T, metadata=meta_sub, design="~group", quiet=True)
dds.deseq2()
stat_res = DeseqStats(
    dds, contrast=["group", "unloading", "control"], quiet=True)
stat_res.summary()
res = stat_res.results_df.copy()
res = res.reset_index()
if "Geneid" in res.columns:
    res = res.rename(columns={"Geneid": "ensembl_id"})
else:
    res = res.rename(columns={"index": "ensembl_id"})

mg = mygene.MyGeneInfo()
ensembl_ids = res["ensembl_id"].dropna().astype(str).tolist()
symbol_map = {}
batch_size = 1000
for i in range(0, len(ensembl_ids), batch_size):
    batch = ensembl_ids[i:i + batch_size]
    try:
        hits = mg.querymany(
            batch,
            scopes="ensembl.gene",
            fields="symbol",
            species="mouse",
            as_dataframe=False,
            verbose=False,
        )
    except Exception as exc:
        print(f"Warning: MyGene lookup failed for batch {i // batch_size + 1}: {exc}")
        continue
    for hit in hits:
        if not hit.get("notfound") and hit.get("symbol"):
            symbol_map[hit["query"]] = hit["symbol"]

res["gene_symbol"] = res["ensembl_id"].map(symbol_map)
missing_symbols = res["gene_symbol"].isna().sum()
if missing_symbols:
    print(f"Warning: {missing_symbols} disuse genes could not be mapped to symbols.")
    res["gene_symbol"] = res["gene_symbol"].fillna(res["ensembl_id"])

OUT.mkdir(exist_ok=True)
res.to_csv(OUT / "disuse_deg.csv", index=False)
print("Saved disuse DEGs:", res.shape)
print(res.dropna(subset=["padj"]).sort_values("padj").head(15).to_string())
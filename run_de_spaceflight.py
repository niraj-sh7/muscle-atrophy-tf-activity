from pathlib import Path

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE = Path(__file__).resolve().parent
DATA = BASE / "data"
OUT = BASE

cnt = pd.read_csv(DATA / "expr_df_factors.csv", index_col=0)
fac = pd.read_csv(DATA / "factors.csv", index_col=0)
print("Counts (samples x genes):", cnt.shape)
cnt = cnt.round().astype(int)
fac = fac.loc[cnt.index]
fac["group"] = np.where(fac["factor"] == 0, "FLT", "GC")
print(fac["group"].value_counts())

dds = DeseqDataSet(counts=cnt, metadata=fac, design="~group", quiet=True)
dds.deseq2()
stat_res = DeseqStats(dds, contrast=["group", "FLT", "GC"], quiet=True)
stat_res.summary()
res = stat_res.results_df.copy()

sigs = pd.read_csv(DATA / "sigs.csv")
sigs = sigs.rename(columns={sigs.columns[0]: "ensembl_id"})
gmap1 = sigs[["ensembl_id", "gene_symbol"]]
gse145480 = DATA / "GSE145480_counts.txt"
if gse145480.exists():
    gmap2 = pd.read_csv(gse145480, sep="\t", usecols=["ensembl_id", "gene_symbol"])
    gmap = pd.concat([gmap1, gmap2]).drop_duplicates("ensembl_id").set_index("ensembl_id")
else:
    print(f"Warning: {gse145480} not found; using sigs.csv gene symbols only.")
    gmap = gmap1.drop_duplicates("ensembl_id").set_index("ensembl_id")
res["gene_symbol"] = gmap.reindex(res.index)["gene_symbol"]
res = res.dropna(subset=["gene_symbol"]).reset_index().rename(
    columns={"index": "ensembl_id"})
print("Spaceflight all-gene results:", res.shape)
geo = OUT / "geo_data"
geo.mkdir(exist_ok=True)
res.to_csv(geo / "spaceflight_deg_full.csv", index=False)
print(res.dropna(subset=["padj"]).sort_values("padj").head(15)[
    ["gene_symbol", "log2FoldChange", "stat", "padj"]
].to_string())

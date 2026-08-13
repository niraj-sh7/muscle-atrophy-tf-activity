from pathlib import Path

import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

BASE = Path(__file__).resolve().parent
DATA = BASE / "data"
OUT = BASE / "geo_data"
src = DATA / "GSE145480_counts.txt"
if not src.exists():
    raise FileNotFoundError(
        f"Missing {src}. This script needs the raw GSE145480 count table."
    )

cnt = pd.read_csv(src, sep="\t")
print("Raw counts shape:", cnt.shape)
gene_map = cnt[["ensembl_id", "gene_symbol"]].drop_duplicates("ensembl_id").set_index("ensembl_id")
mat = cnt.drop(columns=["gene_symbol"]).set_index("ensembl_id").astype(int)

samples = mat.columns.tolist()
ages = [int(s.split("m_")[0]) for s in samples]
meta = pd.DataFrame({"sample": samples, "age": ages}).set_index("sample")
keep_samples = [s for s, a in zip(samples, ages) if a in (8, 28)]
meta_sub = meta.loc[keep_samples].copy()
meta_sub["group"] = np.where(meta_sub["age"] == 28, "old", "young")
mat_sub = mat[keep_samples]

mat_sub = mat_sub[mat_sub.sum(axis=1) >= 10]
print("Filtered matrix:", mat_sub.shape, "samples:", meta_sub.shape)

dds = DeseqDataSet(counts=mat_sub.T, metadata=meta_sub, design="~group", quiet=True)
dds.deseq2()
stat_res = DeseqStats(dds, contrast=["group", "old", "young"], quiet=True)
stat_res.summary()
res = stat_res.results_df.copy()
res["gene_symbol"] = gene_map.reindex(res.index)["gene_symbol"]
res = res.dropna(subset=["gene_symbol"])
res = res.reset_index().rename(columns={"index": "ensembl_id"})
OUT.mkdir(exist_ok=True)
res.to_csv(OUT / "sarcopenia_deg.csv", index=False)
print("Saved sarcopenia DEGs:", res.shape)
print(res[["gene_symbol", "log2FoldChange", "stat", "padj"]].dropna().sort_values("padj").head(15).to_string())

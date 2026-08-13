import os
import pandas as pd
import numpy as np
from scipy import stats

OUT = "/home/user/workspace/results"
                                           
import os
if not os.path.exists(OUT):
    OUT = os.path.join(os.path.dirname(__file__), "results")
nes = pd.read_csv(os.path.join(OUT, "tf_activity_NES_wide.csv"), index_col=0)
pv = pd.read_csv(os.path.join(OUT, "tf_activity_padj_wide.csv"), index_col=0)
print("Loaded NES:", nes.shape)

                                          
NES_CUT = 1.5
P_CUT = 0.05
conditions = nes.columns.tolist()
sig = nes.abs() >= NES_CUT                       
sig_strict = (nes.abs() >= NES_CUT) & (pv <= P_CUT)
n_sig = sig.sum(axis=1)

def classify(row):
    n = row.sum()
    if n == 3:
        return "Pan-Atrophy"
    if n == 2:
        return "Partial"
    if n == 1:
        return "Condition-Specific"
    return "None"

cls = sig.apply(classify, axis=1)
class_df = nes.copy()
class_df["padj_min"] = pv.min(axis=1)
class_df["n_sig_conditions"] = n_sig
class_df["category"] = cls
                                     
def direction_str(row):
    parts = []
    for c in conditions:
        if sig.loc[row.name, c]:
            sign = "+" if nes.loc[row.name, c] > 0 else "-"
            parts.append(f"{sign}{c}")
    return ",".join(parts)
class_df["sig_directions"] = class_df.apply(direction_str, axis=1)

                                      
def concordance(tf_name):
    pan = class_df[class_df["category"] == "Pan-Atrophy"]
    if tf_name not in pan.index:
        return "N/A"
    signs = [np.sign(nes.loc[tf_name, c]) for c in conditions]
    if len(set(signs)) == 1:
        return "Concordant"
    return "Sign-Discordant"

class_df["concordance"] = class_df.index.map(concordance)

sign_discord = class_df[
    (class_df["category"] == "Pan-Atrophy") & (class_df["concordance"] == "Sign-Discordant")
][conditions + ["concordance"]]
sign_discord.to_csv(os.path.join(OUT, "pan_atrophy_sign_discordant.csv"))
print("\nSign-discordant Pan-Atrophy TFs:")
print(sign_discord.round(2).to_string())

concordant_pan = class_df[
    (class_df["category"] == "Pan-Atrophy") & (class_df["concordance"] == "Concordant")
]
print(f"\nConcordant Pan-Atrophy TFs: {len(concordant_pan)}")
print(f"Sign-Discordant Pan-Atrophy TFs: {len(sign_discord)}")

class_df.to_csv(os.path.join(OUT, "tf_classification.csv"))

print("\nCategory counts:")
counts = class_df["category"].value_counts()
print(counts)
counts.to_csv(os.path.join(OUT, "category_counts.csv"))

                                                      
corr = nes.corr(method="pearson")
corr.to_csv(os.path.join(OUT, "correlation_matrix.csv"))
print("\nPearson correlation between conditions:")
print(corr.round(3))

                                                
ATROPHY_TFS = {
    "Foxo1", "Foxo3", "Foxo4", "Foxo6",
    "Nfkb1", "Nfkb2", "Rela", "Relb",
    "Atf3", "Atf4", "Ddit3",                                              
    "Tp53", "Trp53",
    "Stat3",
    "Smad2", "Smad3", "Smad4",
    "Foxk1", "Foxk2",
    "Myod1", "Myog",
    "Igf1",                                                                               
}

pan = set(class_df[class_df["category"] == "Pan-Atrophy"].index)
all_tfs = set(class_df.index)
known_in_pan = ATROPHY_TFS & pan
known_in_other = ATROPHY_TFS & (all_tfs - pan)
unknown_in_pan = pan - ATROPHY_TFS
unknown_in_other = (all_tfs - pan) - ATROPHY_TFS

table = [[len(known_in_pan), len(unknown_in_pan)],
         [len(known_in_other), len(unknown_in_other)]]
odds, p = stats.fisher_exact(table, alternative="greater")
report_lines = [
    "Fisher exact test: enrichment of canonical atrophy-pathway TFs",
    f"  TFs tested: {len(all_tfs)}",
    f"  Pan-atrophy TFs: {len(pan)}",
    f"  Canonical atrophy TFs (FoxO/NFkB/ATF/IGF axis) in regulon: {len(ATROPHY_TFS & all_tfs)}",
    f"  Of those, in Pan-Atrophy class: {sorted(known_in_pan)}",
    f"  Contingency table (known/unknown x pan/other): {table}",
    f"  Odds ratio: {odds:.3f}",
    f"  One-sided p-value: {p:.3e}",
]
report = "\n".join(report_lines)
print("\n" + report)
with open(os.path.join(OUT, "fisher_atrophy_test.txt"), "w") as fh:
    fh.write(report + "\n")

                                               
pan_df = class_df.loc[class_df["category"] == "Pan-Atrophy"].sort_values(
    by=conditions, key=lambda s: s.abs(), ascending=False)
pan_df.to_csv(os.path.join(OUT, "pan_atrophy_TFs.csv"))
print("\nPan-atrophy TFs:")
print(pan_df.head(25).round(2).to_string())
import os
import pandas as pd
import numpy as np

OUT = "/home/user/workspace/results"
if not os.path.exists(OUT):
    OUT = os.path.join(os.path.dirname(__file__), "results")

nes = pd.read_csv(os.path.join(OUT, "tf_activity_NES_wide.csv"), index_col=0)
conditions = nes.columns.tolist()

PRIMARY_CUT = 1.5
THRESHOLDS = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5]

def classify_at_threshold(nes_df, cut):
    sig = nes_df.abs() >= cut
    n_sig = sig.sum(axis=1)
    cls = n_sig.map(lambda n: {3: "Pan-Atrophy", 2: "Partial",
                                1: "Condition-Specific"}.get(n, "None"))
    return cls

rows = []
pan_sets = {}
for cut in THRESHOLDS:
    cls = classify_at_threshold(nes, cut)
    counts = cls.value_counts().reindex(
        ["Pan-Atrophy", "Partial", "Condition-Specific", "None"], fill_value=0
    )
    pan_sets[cut] = set(cls[cls == "Pan-Atrophy"].index)
    rows.append({
        "NES_threshold": cut,
        "Pan-Atrophy": counts["Pan-Atrophy"],
        "Partial": counts["Partial"],
        "Condition-Specific": counts["Condition-Specific"],
        "None": counts["None"],
        "Pan-Atrophy_pct": round(100 * counts["Pan-Atrophy"] / len(cls), 1),
    })

sens_df = pd.DataFrame(rows)
sens_df.to_csv(os.path.join(OUT, "nes_threshold_sensitivity.csv"), index=False)
print("NES threshold sensitivity:")
print(sens_df.to_string(index=False))

primary_pan = pan_sets[PRIMARY_CUT]
overlap_rows = []
for cut, pan in pan_sets.items():
    inter = len(pan & primary_pan)
    union = len(pan | primary_pan)
    jaccard = inter / union if union else np.nan
    overlap_rows.append({
        "NES_threshold": cut,
        "n_pan_atrophy": len(pan),
        "n_overlap_with_primary(1.5)": inter,
        "jaccard_vs_primary": round(jaccard, 3),
    })

overlap_df = pd.DataFrame(overlap_rows)
overlap_df.to_csv(
    os.path.join(OUT, "nes_threshold_pan_atrophy_overlap.csv"), index=False
)
print("\nPan-Atrophy set overlap vs primary threshold (|NES| >= 1.5):")
print(overlap_df.to_string(index=False))

CORE_TFS = ["Nr3c1", "Foxo1", "Smarca4", "Pml", "Ing4"]
print("\nCore pan-atrophy TF robustness across thresholds:")
for tf in CORE_TFS:
    if tf not in nes.index:
        continue
    present_at = [cut for cut, pan in pan_sets.items() if tf in pan]
    print(f"  {tf}: Pan-Atrophy at thresholds {present_at}")

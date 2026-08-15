"""Pairwise TF-activity correlation scatter plots with 95% confidence intervals.

Addresses Reviewer #2 Point 11: "Correlation analyses would be strengthened by
reporting confidence intervals and providing scatter plots to better illustrate
the relationships among conditions."

Produces a 3-panel supplementary figure (Spaceflight vs Sarcopenia, Spaceflight
vs Disuse Atrophy, Sarcopenia vs Disuse Atrophy), each showing:
  - scatter of NES values across all 732 TFs
  - OLS regression line with shaded 95% CI band (seaborn regplot)
  - Pearson r, 95% CI (Fisher z-transform), and p-value annotated on the panel

Output:
  figures/SFig_correlation_scatter.png   (300 DPI)
  results/correlation_CIs.csv            (r, CI bounds, p-value per pair)
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

OUT = "/home/user/workspace/results"
FIG = "/home/user/workspace/figures"
if not os.path.exists(OUT):
    OUT = os.path.join(os.path.dirname(__file__), "results")
    FIG = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(FIG, exist_ok=True)

nes = pd.read_csv(os.path.join(OUT, "tf_activity_NES_wide.csv"), index_col=0)
conditions = nes.columns.tolist()  # e.g. Spaceflight, Sarcopenia, Disuse_Atrophy

def pearson_ci(x, y, alpha=0.05):
    """Pearson r with Fisher z-transform 95% CI."""
    r, p = stats.pearsonr(x, y)
    n = len(x)
    z = np.arctanh(r)
    se = 1 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf(1 - alpha / 2)
    lo, hi = np.tanh(z - z_crit * se), np.tanh(z + z_crit * se)
    return r, lo, hi, p, n

pairs = [
    (conditions[0], conditions[1]),
    (conditions[0], conditions[2]),
    (conditions[1], conditions[2]),
]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
ci_rows = []

for ax, (c1, c2) in zip(axes, pairs):
    x, y = nes[c1].values, nes[c2].values
    r, lo, hi, p, n = pearson_ci(x, y)
    ci_rows.append({
        "Comparison": f"{c1} vs {c2}",
        "Pearson_r": round(r, 4),
        "CI_lower": round(lo, 4),
        "CI_upper": round(hi, 4),
        "p_value": p,
        "n_TFs": n,
    })

    sns.regplot(
        x=x, y=y, ax=ax,
        scatter_kws={"s": 10, "alpha": 0.35, "color": "steelblue"},
        line_kws={"color": "firebrick", "lw": 1.5},
        ci=95,
    )
    ax.axhline(0, color="grey", lw=0.5, ls="--")
    ax.axvline(0, color="grey", lw=0.5, ls="--")
    ax.set_xlabel(f"{c1} NES")
    ax.set_ylabel(f"{c2} NES")
    p_str = f"{p:.2e}" if p < 0.001 else f"{p:.3f}"
    ax.set_title(
        f"{c1} vs {c2}\nr = {r:.3f} (95% CI: {lo:.3f}-{hi:.3f}), p = {p_str}",
        fontsize=10
    )

plt.tight_layout()
fig_path = os.path.join(FIG, "SFig_correlation_scatter.png")
plt.savefig(fig_path, dpi=300, bbox_inches="tight")
print(f"Saved figure to {fig_path}")

ci_df = pd.DataFrame(ci_rows)
ci_df.to_csv(os.path.join(OUT, "correlation_CIs.csv"), index=False)
print("\nPearson correlations with 95% CI:")
print(ci_df.to_string(index=False))

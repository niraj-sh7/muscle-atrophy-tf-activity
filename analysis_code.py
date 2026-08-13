"""
Analysis code: comparative TF activity inference across three muscle atrophy
models (spaceflight / sarcopenia / disuse atrophy).

This script reproduces every result, table and figure used in the
manuscript "Transcription Factor Activity Divergence Across Muscle Atrophy
Conditions: A Comparative Analysis of Spaceflight, Aging, and Disuse Using
decoupleR" (Shah, in preparation).

Workflow:

  Step 1. Differential expression with PyDESeq2 for each condition
          - Spaceflight  : OSD-576 / RR-23 (FLT vs GC, mouse tibialis anterior)
          - Sarcopenia   : GSE145480 (28-mo vs 8-mo, mouse gastrocnemius)
          - Disuse       : GSE273092 (10-d unloading vs control, gastrocnemius)
  Step 2. Build a gene-by-condition Wald-statistic signature matrix
  Step 3. Download the CollecTRI mouse regulon via decoupler-py
  Step 4. Run the ULM (univariate linear model) method on the signature
  Step 5. Classify TFs (Pan-Atrophy / Partial / Condition-Specific / None)
  Step 6. Compute Pearson/Spearman correlations and Fisher exact enrichment
  Step 7. Generate Figures 1-6 and Supplementary Tables 1-2

Software environment:
  python 3.12; pandas 2.x; numpy 2.x; pydeseq2 0.5.4; decoupler 2.1.6;
  anndata 0.10.x; matplotlib 3.9; seaborn 0.13; networkx 3.3;
  matplotlib_venn 1.1; adjustText 1.3; scipy 1.13
Inputs (workspace/):
  - sigs.csv                  spaceflight DEG (significant subset)
  - expr_df_factors.csv       18 sample x ~12k gene OSD-576 count matrix
  - factors.csv               experimental factor assignments for OSD-576
  - geo_data/GSE145480_counts.txt           sarcopenia raw counts
  - geo_data/GSE273092_featureCounts.txt    disuse raw counts

Outputs (workspace/):
  - results/tf_activity_NES_wide.csv          724 TF x 3 condition NES matrix
  - results/tf_activity_padj_wide.csv         BH-adjusted ULM p-values
  - results/tf_classification.csv             classified TFs with directions
  - results/correlation_matrix.csv            Pearson r between conditions
  - results/category_counts.csv               category sizes
  - results/fisher_atrophy_test.txt           Fisher exact test report
  - figures/Figure1.png ... Figure6.png       300-DPI publication figures
  - supplementary_table_1.csv                 full TF activity table
  - supplementary_table_2.csv                 pan-atrophy TFs only

"""

import subprocess
import sys
from pathlib import Path

WORKSPACE = Path("/home/user/workspace")

PIPELINE = [
    "run_de_spaceflight.py",  
    "run_de_sarcopenia.py",   
    "run_de_disuse.py",     
    "run_tf_inference.py", 
    "run_classify.py", 
    "make_figures.py", 
]


def run(script: str) -> None:
    """Run a sub-script and stream its stdout/stderr."""
    print(f"\n=== Running {script} ===", flush=True)
    rc = subprocess.run(
        [sys.executable, str(WORKSPACE / script)],
        check=True,
    ).returncode
    print(f"=== Finished {script} (rc={rc}) ===")


if __name__ == "__main__":
    for s in PIPELINE:
        run(s)
    print("\nPipeline complete.")

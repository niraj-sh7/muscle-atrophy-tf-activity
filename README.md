# Muscle Atrophy TF Activity

Code for **"Transcription Factor Activity Divergence Across Muscle Atrophy Conditions: A Comparative Analysis of Spaceflight, Aging, and Disuse Using decoupleR"** (Shah, submitted to *PLOS ONE*).

This repository reproduces every result, table, and figure in the manuscript: a comparative transcription-factor (TF) activity analysis of mouse skeletal-muscle atrophy across spaceflight, sarcopenia, and disuse, using PyDESeq2 for differential expression and decoupleR (ULM method) against the CollecTRI mouse regulon for TF activity inference.

## Data

This study used only publicly available data. No new data were generated. Raw counts are **not included** in this repository — download them from their original sources below before running the pipeline.

| Condition | Dataset | Source |
|---|---|---|
| Spaceflight | OSD-576 (Rodent Research-23, RR-23) | [NASA Open Science Data Repository](https://osdr.nasa.gov/bio/repo/data/studies/OSD-576) |
| Sarcopenia | GSE145480 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145480) |
| Disuse atrophy | GSE273092 | [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273092) |

Place downloaded files into local `data/` (OSD-576) and `geo_data/` (GSE145480, GSE273092) folders matching the paths expected by the `run_de_*.py` scripts.

## Setup

```bash
pip install -r requirements.txt
```

Requires Python 3.12. Key dependencies: `pydeseq2`, `decoupler`, `anndata`, `pandas`, `numpy`, `scipy`.

## Reproduce

Run the full pipeline with:

```bash
python analysis_code.py
```

Or run each step manually, in order:

1. `run_de_spaceflight.py` — spaceflight DEGs (PyDESeq2, FLT vs GC)
2. `run_de_sarcopenia.py` — sarcopenia DEGs (28-mo vs 8-mo)
3. `run_de_disuse.py` — disuse DEGs (10-day hindlimb unloading vs control)
4. `run_tf_inference.py` — TF activity inference (decoupleR ULM against CollecTRI mouse regulon)
5. `run_classify.py` — TF classification (Pan-Atrophy / Partial / Condition-Specific), correlation analysis, Fisher exact test

Results are saved to `results/` (not included in this repo — regenerate by running the pipeline against the public data above).

## Citation

If you use this code, please cite the associated manuscript (citation to be updated upon publication).

## License

MIT

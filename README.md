# Gene Expression Analyzer

A Python pipeline for simulating and analyzing differential gene expression between tumor and normal samples, built as a bioinformatics portfolio project.

## What it does
- Simulates gene expression data for 10 genes across 6 tumor and 6 normal samples
- Calculates fold change (Log2FC) between groups
- Runs independent t-tests (SciPy) to test statistical significance per gene
- Applies Benjamini-Hochberg multiple testing correction (statsmodels) to control false positives
- Identifies top differentially expressed candidate genes
- Visualizes results with a volcano plot

## Tools used
Python, NumPy, Pandas, Matplotlib, SciPy, statsmodels

## How to run
```bash
python3 notebooks/04_final_analysis.py
```

## Results
See `results/final_result.csv` for the full results table and `results/volcano_plot.png` for the visualization.
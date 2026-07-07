import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
#Recreate the same dataset
np.random.seed(42)
genes = ["BRCA1", "TP53", "MYC", "EGFR", "KRAS", "PTEN", "AKT1", "VEGFA", "CDKN2A", "MDM2"]
tumor_data = np.random.normal(loc=8, scale=1.5, size=(10, 6))
normal_data = np.random.normal(loc=5, scale=1.5, size=(10, 6))
tumor_columns = [f"Tumor_{i+1}"for i in range(6)]
normal_columns = [f"Normal_{i+1}"for i in range(6)]
tumor_df = pd.DataFrame(tumor_data, index=genes, columns=tumor_columns)
normal_df = pd.DataFrame(normal_data, index=genes, columns=normal_columns)
expression_df = pd.concat([tumor_df, normal_df], axis=1)
#Fold change
expression_df["Tumor_Mean"] = expression_df[tumor_columns].mean(axis=1)
expression_df["Normal_Mean"] = expression_df[normal_columns].mean(axis=1)
expression_df["Log2FC"] = np.log2(expression_df["Tumor_Mean"]/ expression_df["Normal_Mean"])
#T-test + p-values
p_values = []
for gene in genes:
    tumor_values = expression_df.loc[gene, tumor_columns]
    normal_values = expression_df.loc[gene, normal_columns]
    t_stat, p_val = stats.ttest_ind(tumor_values, normal_values)
    p_values.append(p_val)
expression_df["P_Value"] = p_values
#Multiple testing correction
reject, p_adjusted, _, _ = multipletests(expression_df["P_Value"], method="fdr_bh")
expression_df["P_Adjusted"] = p_adjusted
expression_df["Significant"] = reject
#Build final clean results table
results = expression_df[["Tumor_Mean", "Normal_Mean", "Log2FC", "P_Value", "P_Adjusted", "Significant"]]
#Sort by adjusted and p-value(most significant genes first)
results_sorted = results.sort_values(by="P_Adjusted")
print(results_sorted)
#save results table to CSV file
results_sorted.to_csv("results/final_result.csv", index=True)
#identify top differentially expressed genes
top_genes = results_sorted[(results_sorted["Significant"] == True) & (abs(results_sorted["Log2FC"]) > 0.5)]
print("\nTop Differentially Expressed Genes:")
print(top_genes)
import matplotlib.pyplot as plt
#Calculate -log10(p-value) for the y-axis
results_sorted["NegLog10P"] = -np.log10(results_sorted["P_Value"])
#set colors: red if significant, gray if not
colors = ["crimson" if sig else "gray" for sig in results_sorted["Significant"]]
plt.figure(figsize=(9, 7))
plt.scatter(results_sorted["Log2FC"], results_sorted["NegLog10P"], c=colors)
#Label each point with its gene names
for gene in results_sorted.index:
    plt.annotate(gene, 
                 (results_sorted.loc[gene, "Log2FC"], results_sorted.loc[gene, "NegLog10P"]), 
                 fontsize=8, xytext=(5, 5), textcoords="offset points")
plt.axhline(-np.log10(0.05), color="blue", linestyle="--", linewidth=1, label="p 0.05 threshold")
plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10(P-Value)")
plt.title("Volcano Plot: Tumor vs Normal Gene Expression")
plt.legend()
plt.tight_layout()
plt.savefig("results/volcano_plot.png")
plt.show()

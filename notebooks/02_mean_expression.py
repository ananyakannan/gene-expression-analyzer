import numpy as np
import pandas as pd
#Recreate the same dataset (same seed = same numbers as before)
np.random.seed(42)
genes = ["BRCA1", "TP53", "MYC", "EGFR", "KRAS", "PTEN", "AKT1", "VEGFA", "CDKN2A", "MDM2"]
tumor_data = np.random.normal(loc=8, scale=1.5, size=(10, 6))
normal_data = np.random.normal(loc=5, scale=1.5, size=(10, 6))
tumor_columns = [f"Tumor_{i+1}" for i in range(6)]
normal_columns = [f"Normal_{i + 1}" for i in range(6)]
tumor_df = pd.DataFrame(tumor_data, index=genes, columns=tumor_columns)
normal_df = pd.DataFrame(normal_data, index=genes, columns=normal_columns)
expression_df = pd.concat([tumor_df, normal_df], axis=1)
#Calculate mean expression per gene, across tumor samples only
expression_df["Tumor_Mean"] = expression_df[tumor_columns].mean(axis=1)
#Calculate mean expression per gene, across normal samples only
expression_df["Normal_Mean"] = expression_df[normal_columns].mean(axis = 1)
print(expression_df[["Tumor_Mean", "Normal_Mean"]])
#Calculate log2 fold change per gene
expression_df["Log2FC"] = np.log2(expression_df["Tumor_Mean"] / expression_df["Normal_Mean"])
print(expression_df[["Tumor_Mean", "Normal_Mean", "Log2FC"]])
import matplotlib.pyplot as plt
#Set up bar positions
x = np.arange(len(genes)) #position for each gene on the x-axis
width = 0.35 #width of each bars
fig, ax = plt.subplots(figsize=(10, 6))
#Plot tumor bars
ax.bar(x - width/2, expression_df["Tumor_Mean"], width, label="Tumor", color="crimson")
#Plot normal bars
ax.bar(x + width/2, expression_df["Normal_Mean"], width, label="Normal", color="steelblue")
#Labeling
ax.set_xlabel("Gene")
ax.set_ylabel("Mean Expression")
ax.set_title("Tumor vs Normal Gene Expression")
ax.set_xticks(x)
ax.set_xticklabels(genes, rotation=45)
ax.legend()
plt.tight_layout()
plt.savefig("results/tumor_vs_normal_expression.png")
plt.show()

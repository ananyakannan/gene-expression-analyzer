import numpy as np
import pandas as pd
from scipy import stats
#Recreate the same dataset
np.random.seed(42)
genes = ["BRCA1", "TP53", "MYC", "EGFR", "KRAS", "PTEN", "AKT1", "VEGFA", "CDKN2", "MDM2"]
tumor_data = np.random.normal(loc=8, scale=1.5, size=(10, 6))
normal_data = np.random.normal(loc=5, scale=1.5, size=(10, 6))
tumor_columns = [f"Tumor_{i+1}"for i in range(6)]
normal_columns = [f"Normal_{i + 1}"for i in range(6)]
tumor_df = pd.DataFrame(tumor_data, index=genes, columns=tumor_columns)
normal_df = pd.DataFrame(normal_data, index=genes, columns=normal_columns)
expression_df = pd.concat([tumor_df, normal_df], axis=1)
#Run for t-test for each gene
p_values = []
for gene in genes:
    tumor_values = expression_df.loc[gene, tumor_columns]
    normal_values = expression_df.loc[gene, normal_columns]
    t_stat, p_val = stats.ttest_ind(tumor_values, normal_values)
    p_values.append(p_val)
expression_df["P_Value"] = p_values
print(expression_df[["P_Value"]])
from statsmodels.stats.multitest import multipletests
#Apply multiple testing correction (Benjamin-Hochberg FDR model)
reject, p_adjusted, _, _ = multipletests(expression_df["P_Value"], method="fdr_bh")
expression_df["P_Adjusted"] = p_adjusted
expression_df["Significant"] = reject
print(expression_df[["P_Value", "P_Adjusted", "Significant"]])



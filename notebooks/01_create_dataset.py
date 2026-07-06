import numpy as np
import pandas as pd
#Set a seed so our "random" numbers are reproducible
np.random.seed(42)
#Define gene names
genes = ["BRCA1", "TP53", "MYC", "EGFR", "KRAS", "PTEN", "AKT1", "VEGFA", "CDKN2A", "MDM2"]
#Simulate expression data: 10 genes x 6 tumor samples
tumor_data = np.random.normal(loc=8, scale=1.5, size=(10, 6))
#Simulate expression data: 10 genes x 6 normal samples
normal_data = np.random.normal(loc=5, scale=1.5, size=(10, 6))
#Create column names for tumor and normal samples
tumor_columns = [f"Tumor_{i+1}" for i in range(6)]
normal_columns = [f"Normal_{i+1}" for i in range(6)]
#Combine tumor and normal data into one DataFrame
tumor_df = pd.DataFrame(tumor_data, index=genes, columns=tumor_columns)
normal_df = pd.DataFrame(normal_data, index=genes, columns=normal_columns)
#Merge both into a single expression table
expression_df = pd.concat([tumor_df, normal_df], axis=1)
print(expression_df)
#Check the shape of the DataFrame
print("Shape:", expression_df.shape)
#Look at the first few rows
print(expression_df.head())
#Get summary statistics
print(expression_df.describe())
    

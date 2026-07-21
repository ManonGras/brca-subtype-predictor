import pandas as pd
df = pd.read_csv("data/processed/tcga_train.csv", nrows=5)
print(f"Colonnes dupliquees : {df.columns.duplicated().sum()}")
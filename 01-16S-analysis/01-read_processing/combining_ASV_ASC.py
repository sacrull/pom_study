import pandas as pd
df1 = pd.read_table("ASV_taxids",delimiter='\t', header= None)
df1.columns =['ASV', 'Taxa_level', 'TaxID']
print(df1.head())

df2 = pd.read_table("all2id.txt",delimiter='\t', header= None)
df2.columns =['Asc', 'TaxID']
print(df2.head())

combined = pd.merge(df1, df2, on ='TaxID')

print(combined['ASV'].unique())

combined.to_csv('ASVs_ASC.tsv', sep="\t")
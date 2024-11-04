import pandas as pd
df1 = pd.read_table("taxids",delimiter='\t', header= None)
df1.columns =['taxid']
print(df1.head())

df2 = pd.read_table("rankedlineage_clean2",delimiter='\t', header= None)
df2.columns =['taxid', 'lineages']
print(df2.head())

combined = pd.merge(df1, df2, on ='taxid', how = 'left')

combined.to_csv('lineage', sep="\t", index=False, header=False)

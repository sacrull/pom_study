import pandas as pd

# Load the data from the file (replace 'ezbiocloud_id_taxonomy.txt' with your actual file path)
df = pd.read_csv('ezbiocloud_id_taxonomy.txt', header=None, names=['ID', 'Taxonomy'], sep='\t')

# Split the taxonomy column by semicolon into separate columns
taxonomy_split = df['Taxonomy'].str.split(';', expand=True)

# Create a new dataframe to store the transformed results
result = []

# Iterate over each row in the original dataframe
for index, row in df.iterrows():
    # Get the base ID (don't change it, just keep it as it is)
    base_id = row['ID']
    # Create an empty string to store the lineage path
    lineage = ""
    
    # For each taxonomy level in the row, create a new row with the ID, corresponding taxon, and lineage
    for i, taxon in enumerate(taxonomy_split.iloc[index]):
        if taxon is not None and taxon != "":
            # The hierarchical ID stays the same as the base ID
            hierarchical_id = base_id * (10 ** i)  # Not needed, base ID remains fixed
            
            # Update the lineage path (add the taxon to the path with the '|' delimiter)
            if lineage == "":
                lineage = taxon
            else:
                lineage += "|" + taxon
            
            # Append the result as a list: [base ID, taxon, lineage]
            result.append([base_id, taxon, lineage])
        else:
            continue

# Create a new DataFrame from the result
result_df = pd.DataFrame(result, columns=['ID', 'Taxon', 'Lineage'])

df2 = pd.read_csv('rankedlineage_clean', header=None, names=['taxID', 'Taxon'], sep='\t')

combined = pd.merge(result_df, df2, on ='Taxon')

# Rearranging the columns and modifying the Lineage
combined = combined[['taxID', 'Lineage']]  # Drop 'ID' column and reorder columns
combined = combined.drop_duplicates()
combined['Lineage'] = combined['Lineage'] + '|'  # Add '|' at the end of the lineage

# Print the first few rows to verify the result
print(combined.head())

# Save the result to a file (for example 'final_taxonomy.txt')
combined.to_csv('rankedlineage_clean2', header=False, index=False, sep='\t')


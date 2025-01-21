# Read the sorted file into R
taxonomy_data_sorted <- read.table("taxonomy_bac.rare.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Extract numeric part of ASV identifier (remove "ASV" prefix)
taxonomy_data_sorted$ASV_numeric <- as.numeric(sub("ASV", "", taxonomy_data_sorted$V1))

# Sort the data by the numeric ASV values
taxonomy_data_sorted_sorted <- taxonomy_data_sorted[order(taxonomy_data_sorted$ASV_numeric), ]

# Drop the temporary numeric column
taxonomy_data_sorted_sorted$ASV_numeric <- NULL

# Write the sorted data back to a new file
write.table(taxonomy_data_sorted_sorted, "taxonomy_rare.sort.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Inspect the first few rows of the sorted data
head(taxonomy_data_sorted_sorted)

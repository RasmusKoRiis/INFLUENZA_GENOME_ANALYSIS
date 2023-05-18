import pandas as pd
import sys

summary_file = sys.argv[1]
mutation_file = sys.argv[2]
output_file = sys.argv[3]

# Read the summary and mutation tables into pandas DataFrames
summary_table = pd.read_csv(summary_file)
mutation_table = pd.read_csv(mutation_file)

# Select rows containing "HA2" in the "Ref_Name" column from the mutation table
ha2_rows = mutation_table[mutation_table['Ref_Name'].str.contains('HA2')]

# Keep only the required columns
ha2_rows = ha2_rows[['sample', 'Ref_Name', 'Differences']]

# Append the selected rows to the summary table using pd.concat()
updated_summary_table = pd.concat([summary_table, ha2_rows], ignore_index=True)

# Save the updated summary table to a new CSV file
updated_summary_table.to_csv(output_file, index=False)


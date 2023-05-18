import pandas as pd
import sys
import os.path

csv_file = sys.argv[1]
tab_file = sys.argv[2]
output_file_long = sys.argv[3]
output_file_short = sys.argv[4]

# Load the tab-separated file
df_tab = pd.read_csv(tab_file , sep='\t')

# Remove "_bam" from the Sample column
df_tab['Sample'] = df_tab['Sample'].str.replace('_bam', '')

# Rearrange columns
df_tab = df_tab[['Sample'] + list(df_tab.columns[:-1])]
df_tab = df_tab.rename(columns={'Sample': 'sample'})

# Load the comma-separated file
df_comma = pd.read_csv(csv_file)

# Merge the two dataframes based on the Sample column
df_final = pd.merge(df_tab, df_comma[['sample', 'subtype']], on='sample', how='left')

# Merge the two dataframes based on the Sample column
df_final = pd.merge(df_final, df_comma[['sample', 'type']], on='sample', how='left')

# Rename the subtype column to Subtype
df_final = df_final.rename(columns={'subtype': 'Subtype'})

# Rename the subtype column to Subtype
df_final = df_final.rename(columns={'type': 'Type'})

# Remove columns for short summary file
df_short = df_final.drop(['Ref_Len', 'Breadth', 'Min_Depth', "Max_Depth", "Std_Dev", "Above_0.2_Depth", "Above_1_Depth","Above_1.8_Depth", "Variation_Coefficient"], axis=1)

# Save the final dataframe to a CSV file
df_final.to_csv(output_file_long, index=False)
df_short.to_csv(output_file_short, index=False)

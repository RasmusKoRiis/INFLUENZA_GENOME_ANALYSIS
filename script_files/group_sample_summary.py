import pandas as pd
import sys
import numpy as np

csv_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(csv_file)

# Pivot the table to create one line per sample
df_pivot = df.pivot(index='sample', columns='Ref_Name')

# Flatten the column names to create unique column labels
df_pivot.columns = ['_'.join(col).rstrip('_') for col in df_pivot.columns.values]

# List of columns that should be in the dataframe
required_columns = ['Subtype', 
                    'Matches', 
                    'Avg_Depth_A_HA_H1', 'Avg_Depth_A_NA_N1', 'Avg_Depth_A_H1_MP', 'Avg_Depth_A_H1_NP',
                'Thymine_ratio_823_N1','Cytosine_ratio_823_N1','Avg_Depth_A_H1_NS', 'Avg_Depth_A_H1_PA', 'Avg_Depth_A_H1_PB1', 'Avg_Depth_A_H1_PB2',
                'Avg_Depth_A_HA_H3', 'Avg_Depth_A_NA_N2', 'Avg_Depth_A_H3_MP', 'Avg_Depth_A_H3_NP',
                'Avg_Depth_A_H3_NS', 'Avg_Depth_A_H3_PA', 'Avg_Depth_A_H3_PB1', 'Avg_Depth_A_H3_PB2',
                'Avg_Depth_B_VIC_HA', 'Avg_Depth_B_VIC_NA', 
                'Differences_A_HA_H1', 'Differences_A_NA_N1', 'Differences_A_H1_MP', 'Differences_A_H1_NP',
                'Differences_A_H1_NS', 'Differences_A_H1_PA', 'Differences_A_H1_PB1', 'Differences_A_H1_PB2',
                'Differences_A_HA_H3', 'Differences_A_NA_N2', 'Differences_A_H3_MP', 'Differences_A_H3_NP',
                'Differences_A_H3_NS', 'Differences_A_H3_PA', 'Differences_A_H3_PB1', 'Differences_A_H3_PB2',
                'Differences_B_VIC_HA', 'Differences_B_VIC_NA']

# Add missing columns with NaN values
for col in required_columns:
    if not any(col in s for s in df_pivot.columns):
        df_pivot[col] = np.nan

#SUBTYPE

# Reset the index to make sample a regular column
df_pivot = df_pivot.reset_index()

# Find all columns containing 'Subtype'
subtype_cols = [col for col in df_pivot.columns if 'Subtype' in col]

# Merge columns containing 'Subtype'
merged_subtype_col = df_pivot[subtype_cols].mode(axis=1)

# Rename the new column to 'Subtype'
merged_subtype_col.columns = ['Subtype']

# Drop the original columns containing 'Subtype'
df_pivot = df_pivot.drop(subtype_cols, axis=1)

# Concatenate the original data with the new merged column
df_pivot = pd.concat([df_pivot, merged_subtype_col], axis=1)

#MATCHES

# Find all columns containing 'Matches'
matches_cols = [col for col in df_pivot.columns if 'Matches' in col]

# Merge columns containing 'Matches'
merged_matches_col = df_pivot[matches_cols].mode(axis=1)

# Rename the new column to 'Matches'
merged_matches_col.columns = ['Matches']

# Drop the original columns containing 'Matches'
df_pivot = df_pivot.drop(matches_cols, axis=1)

# Concatenate the original data with the new merged column
df_pivot = pd.concat([df_pivot, merged_matches_col], axis=1)

#CLADE

# Reset the index to make sample a regular column
df_pivot = df_pivot.reset_index()

# Find all columns containing 'Clade'
clade_cols = [col for col in df_pivot.columns if 'clade' in col]

# Merge columns containing 'Clade'
merged_clade_col = df_pivot[clade_cols].mode(axis=1)

print(merged_clade_col)

# Rename the new column to 'Clade'
merged_clade_col.columns = ['Clade']

# Drop the original columns containing 'Clade'
df_pivot = df_pivot.drop(clade_cols, axis=1)

# Concatenate the original data with the new merged column
df_pivot = pd.concat([df_pivot, merged_clade_col], axis=1)

#THYMINE

# Reset the index to make sample a regular column
#df_pivot = df_pivot.reset_index()

# Find all columns containing 'Thy'
thy_cols = [col for col in df_pivot.columns if 'thymine_ratio__823' in col]

# Merge columns containing 'Thy'
merged_thy_col = df_pivot[thy_cols].mode(axis=1)

# Rename the new column to 'Thy'
merged_thy_col.columns = ['Thymine_ratio_823_N1']

# Drop the original columns containing 'Thy
df_pivot = df_pivot.drop(thy_cols, axis=1)

# Concatenate the original data with the new merged column
df_pivot = pd.concat([df_pivot, merged_thy_col], axis=1)

#CYTOSINE

# Reset the index to make sample a regular column
#df_pivot = df_pivot.reset_index()

# Find all columns containing 'Thy'
cyt_cols = [col for col in df_pivot.columns if 'cytosine_ratio__823' in col]

# Merge columns containing 'Thy'
merged_cyt_col = df_pivot[cyt_cols].mode(axis=1)
print(merged_thy_col)

# Rename the new column to 'Thy'
merged_cyt_col.columns = ['Cytosine_ratio_823_N1']

# Drop the original columns containing 'Thy
df_pivot = df_pivot.drop(cyt_cols, axis=1)

# Concatenate the original data with the new merged column
df_pivot = pd.concat([df_pivot, merged_cyt_col], axis=1)

print(df_pivot)

# select only columns to keep - complete list
df_pivot = df_pivot.loc[:, ['sample', 'Subtype', 'Matches', 'Clade','Avg_Depth_A_HA_H1', 'Avg_Depth_A_NA_N1', 'Avg_Depth_A_H1_MP', 'Avg_Depth_A_H1_NP',
                'Thymine_ratio_823_N1','Cytosine_ratio_823_N1','Avg_Depth_A_H1_NS', 'Avg_Depth_A_H1_PA', 'Avg_Depth_A_H1_PB1', 'Avg_Depth_A_H1_PB2',
                'Avg_Depth_A_HA_H3', 'Avg_Depth_A_NA_N2', 'Avg_Depth_A_H3_MP', 'Avg_Depth_A_H3_NP',
                'Avg_Depth_A_H3_NS', 'Avg_Depth_A_H3_PA', 'Avg_Depth_A_H3_PB1', 'Avg_Depth_A_H3_PB2',
                'Avg_Depth_B_VIC_HA', 'Avg_Depth_B_VIC_NA', 
                'Differences_A_HA_H1', 'Differences_A_NA_N1', 'Differences_A_H1_MP', 'Differences_A_H1_NP',
                'Differences_A_H1_NS', 'Differences_A_H1_PA', 'Differences_A_H1_PB1', 'Differences_A_H1_PB2',
                'Differences_A_HA_H3', 'Differences_A_NA_N2', 'Differences_A_H3_MP', 'Differences_A_H3_NP',
                'Differences_A_H3_NS', 'Differences_A_H3_PA', 'Differences_A_H3_PB1', 'Differences_A_H3_PB2',
                'Differences_B_VIC_HA', 'Differences_B_VIC_NA']]

# create a list of columns in the desired order
new_order = ['sample', 'Subtype', 'Matches', 'Clade','Thymine_ratio_823_N1','Cytosine_ratio_823_N1','Avg_Depth_A_HA_H1', 'Avg_Depth_A_NA_N1', 'Avg_Depth_A_H1_MP', 'Avg_Depth_A_H1_NP',
             'Avg_Depth_A_H1_NS', 'Avg_Depth_A_H1_PA', 'Avg_Depth_A_H1_PB1', 'Avg_Depth_A_H1_PB2',
             'Avg_Depth_A_HA_H3', 'Avg_Depth_A_NA_N2', 'Avg_Depth_A_H3_MP', 'Avg_Depth_A_H3_NP',
             'Avg_Depth_A_H3_NS', 'Avg_Depth_A_H3_PA', 'Avg_Depth_A_H3_PB1', 'Avg_Depth_A_H3_PB2',
             'Avg_Depth_B_VIC_HA', 'Avg_Depth_B_VIC_NA', 
             'Differences_A_HA_H1', 'Differences_A_NA_N1', 'Differences_A_H1_MP', 'Differences_A_H1_NP',
             'Differences_A_H1_NS', 'Differences_A_H1_PA', 'Differences_A_H1_PB1', 'Differences_A_H1_PB2',
             'Differences_A_HA_H3', 'Differences_A_NA_N2', 'Differences_A_H3_MP', 'Differences_A_H3_NP',
             'Differences_A_H3_NS', 'Differences_A_H3_PA', 'Differences_A_H3_PB1', 'Differences_A_H3_PB2',
             'Differences_B_VIC_HA', 'Differences_B_VIC_NA']

# use the reindex method to rearrange the columns
#df_pivot = df_pivot.reindex(columns=new_order)

# create a dictionary to map the old column names to the new ones
column_names_map = {col: col.replace('Differences', 'Mutations') if 'Differences' in col else col for col in df_pivot.columns}

# use the rename method to change the column names
df_pivot = df_pivot.rename(columns=column_names_map)

# Rename Column
df_pivot = df_pivot.rename(columns={'Matches': 'Fluserver_matches'})

# clade_fix
for index, row in df_pivot.iterrows():
    # Check if the SubType column is empty
    if row['Avg_Depth_A_HA_H1'] < 30 and row['Avg_Depth_A_HA_H3'] < 30 and row['Avg_Depth_B_VIC_HA'] < 30:
         # Add the Ref_Name string to the SubType column
        df.loc[index, 'Clade'] = 'Low_Depth'

df_pivot.to_csv(output_file, index=False)        
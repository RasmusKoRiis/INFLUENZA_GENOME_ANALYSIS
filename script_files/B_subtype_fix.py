import pandas as pd
import sys

csv_file= sys.argv[1]
output_file= sys.argv[2]

# Load the CSV file into a pandas DataFrame
df = pd.read_csv(csv_file)

# B-Subtype_fix
for index, row in df.iterrows():
    # Check if the SubType column is empty
    if pd.isnull(row['Subtype']):
        # Check if the Avg_Depth column is greater than 20
        if row['Avg_Depth'] > 60 and 'B_VIC' in row['Ref_Name'] and 'mixed' not in row['Type']:
            # Add the Ref_Name string to the SubType column
            subtype_text = row['Ref_Name']
            if 'B_VIC' in subtype_text:
                subtype_text = 'B_Victoria'
            df.loc[index, 'Subtype'] = subtype_text    
            df.loc[index, 'Type'] = 'Type_B'
        if 'mixed' in row['Type']:
            df.loc[index, 'Subtype'] = 'mixed'  


# B-Subtype_fix
for index, row in df.iterrows():
    # Check if the SubType column is empty
    if row['Type'] == 'mixed':
        df.loc[index, 'clade'] = 'undetermined'



# Write the modified DataFrame back to a new CSV file
df.to_csv(output_file, index=False)
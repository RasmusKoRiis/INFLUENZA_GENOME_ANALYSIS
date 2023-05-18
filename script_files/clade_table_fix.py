import pandas as pd
import sys

orginal_summary = sys.argv[1]
new_summary = sys.argv[2]

df = pd.read_csv(orginal_summary, sep=";")

# Split "seqName" column into "sample" and "Ref_Name" columns
df[['sample', 'Ref_Name']] = df['seqName'].str.split('_', n=1, expand=True)

# Combine the last character of the "sample" column with the first character of the "Ref_Name" column
df['sample'] = df['sample'] + '_' + df['Ref_Name'].str[:1]

# Remove the first character of the "Ref_Name" column
df['Ref_Name'] = df['Ref_Name'].str[2:]

# Drop the original "ID" column
df.drop(columns=['seqName'], inplace=True)

# Save the final dataframe to a CSV file
df.to_csv(new_summary, index=False)

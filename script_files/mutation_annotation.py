import pandas as pd
import sys

summary_file = sys.argv[1]
flueserver = sys.argv[2]
output_file = sys.argv[3]
lookup = sys.argv[4]

# Read in the dataframes
df2 = pd.read_csv(summary_file)
df1 = pd.read_csv(flueserver)

# Create a new column in df2
df2[lookup] = ''

# Loop through each row in the second DataFrame
for i, row in df2.iterrows():
    # Get the reference name for this row
    ref_name = row['Ref_Name']
    # Check if the reference name is a column in the first DataFrame
    if ref_name in df1.columns:
        # Get the differences for this row
        differences = str(row['Differences'])
        print(type(differences))
        # Split the differences into a list
    
        differences_list = differences.split(';')
           # Initialize an empty list to store the matches for this row
        matches = []
           # Check each difference string against the corresponding column in the first DataFrame
        for difference in differences_list:
                # Split the reference string into a list of substrings
            reference_list = df1[ref_name].str.split(';')
                # Check if any of the substrings match the current difference
            for mut in reference_list:
                for s in mut:
                    if s in difference:
                            # If there is a match, append it to the matches list
                        matches.append(s)
            # If there were any matches, add them to the Matches column for this row
        if len(matches) > 0:
            df2.at[i, lookup] = ';'.join(matches)
        else:
            df2.at[i, lookup] = 'NO FLUSERVER MATCH'

df2.to_csv(output_file, index=False)

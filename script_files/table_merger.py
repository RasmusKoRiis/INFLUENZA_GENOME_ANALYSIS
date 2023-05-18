import sys
import os
import glob
import pandas as pd

# Get the file extension and name from command-line arguments
ext = sys.argv[1]
name = sys.argv[2]
folder_path = sys.argv[3]

# Find all files with the given extension in the specified directory
files = glob.glob(os.path.join(folder_path, f'*{ext}'))

# Read all files into separate dataframes and concatenate them
dfs = [pd.read_csv(f) for f in files]
merged_df = pd.concat(dfs)

# Write the merged dataframe to a file with the same extension in the specified directory
merged_df.to_csv(os.path.join(folder_path, f'{name}'), index=False)

import csv
import sys
import os.path
import pandas as pd

main_csv = sys.argv[1]
adding_csv = sys.argv[2]
output_file = sys.argv[3]
sep_main = sys.argv[4]
sep_adding = sys.argv[5]


look_up_column_main, lookup_column_main_2 = sys.argv[6].split(',')
look_up_column_adding, look_up_column_adding_2 = sys.argv[7].split(',')
adding_column = sys.argv[8]

# Load main file
main_df = pd.read_csv(main_csv, sep=sep_main)

# Load adding file
adding_df = pd.read_csv(adding_csv, sep=sep_adding)

# Rename the subtype column to Subtype
adding_df  = adding_df .rename(columns={look_up_column_adding: look_up_column_main})
adding_df  = adding_df .rename(columns={look_up_column_adding_2: lookup_column_main_2})

# Merge the two dataframes based on the Sample column and Ref_Name
df_final = pd.merge(main_df, adding_df[[look_up_column_main, lookup_column_main_2, adding_column]], on=[look_up_column_main, lookup_column_main_2], how='left')

# Save the final dataframe to a CSV file
df_final.to_csv(output_file, index=False)
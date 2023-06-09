

import pandas as pd
import sys
import numpy as np
from datetime import datetime

# Get the current date and time
start_time = datetime.now()

print("Start time:", start_time)

def add_missing_columns(df, required_columns):
    for col in required_columns:
        if not any(col in s for s in df.columns):
            df[col] = np.nan
    return df

def merge_and_rename_cols(df, search_str, new_col_name):
    cols = [col for col in df.columns if search_str in col]
    merged_col = df[cols].mode(axis=1)
    print(merged_col.head())
    print(merged_col.columns)

    # If more than one mode, write 'MIXED' in the cell
    merged_col = merged_col.apply(lambda row: 'MIXED' if row.count() > 1 else row[0], axis=1)
    
    df = df.drop(cols, axis=1)
    df = pd.concat([df, merged_col.rename(new_col_name)], axis=1)
    return df

def rename_columns(df, column_names_map):
    df = df.rename(columns=column_names_map)
    return df

def classify_quality(df, quality_columns):
    # Select numeric columns only
    numeric_df = df.select_dtypes(include=[np.number])
    
    # Check if each column contains 'Depth' and is in quality_columns
    columns_to_consider = [col for col in numeric_df.columns if 'Depth' in col and any(sub in col for sub in quality_columns)]
    
    df['Average'] = numeric_df[columns_to_consider].mean(axis=1)

    df.loc[df['Average'] < 50, 'Quality'] = 'Bad'
    df.loc[(df['Average'] >= 50) & (df['Average'] < 100), 'Quality'] = 'OK'
    df.loc[df['Average'] >= 100, 'Quality'] = 'Good'
    
    return df

def process_dataframe(df, runname, required_columns, column_names_map, quality_columns):
    df = add_missing_columns(df, required_columns)
    
    df = merge_and_rename_cols(df, 'Subtype', 'Subtype')
    df = merge_and_rename_cols(df, 'FluserverV1', 'FluserverV1')
    df = merge_and_rename_cols(df, 'clade', 'Clade')
    df = merge_and_rename_cols(df, 'thymine_ratio__823', 'Thymine_ratio_823_N1')
    df = merge_and_rename_cols(df, 'cytosine_ratio__823', 'Cytosine_ratio_823_N1')
    
    df['RunName'] = runname
    df['ScriptTimestamp'] = start_time
    df['SampleApprovalStatus'] = ''
    df['SequenceID'] = df['RunName'].astype(str) + "_" + df['sample'].astype(str)

    df = classify_quality(df, quality_columns)

    df = df.loc[:, list(column_names_map.keys())]
    df = rename_columns(df, column_names_map)
    
    return df

def main(csv_file, output_file, runname):
    df = pd.read_csv(csv_file)
    df = df.pivot(index='sample', columns='Ref_Name')
    df.columns = ['_'.join(col).rstrip('_') for col in df.columns.values]
    df = df.reset_index()

    required_columns = ['Subtype', 
                        'FluserverV1', 
                        'Avg_Depth_A_HA_H1', 
                        'Avg_Depth_A_NA_N1', 
                        'Avg_Depth_A_XX_MP', 
                        'Avg_Depth_A_H1_NP',
                        'thymine_ratio_823_N1',
                        'cytosine_ratio_823_N1',
                        'Avg_Depth_A_H1_NS', 
                        'Avg_Depth_A_H1_PA', 
                        'Avg_Depth_A_H1_PB1', 
                        'Avg_Depth_A_H1_PB2',
                        'Avg_Depth_A_HA_H3', 
                        'Avg_Depth_A_NA_N2', 
                        'Avg_Depth_A_H3_NP',
                        'Avg_Depth_A_H3_NS', 
                        'Avg_Depth_A_H3_PA', 
                        'Avg_Depth_A_H3_PB1', 
                        'Avg_Depth_A_H3_PB2',
                        'Avg_Depth_A_HA_H5',  
                        'Avg_Depth_A_HA_H9',  
                        'Avg_Depth_B_VIC_HA', 
                        'Avg_Depth_B_VIC_NA', 
                        'Differences_A_HA_H1', 
                        'Differences_A_NA_N1', 
                        'Differences_A_XX_MP', 
                        'Differences_A_H1_NP',
                        'Differences_A_H1_NS',
                        'Differences_A_H1_PA', 
                        'Differences_A_H1_PB1', 
                        'Differences_A_H1_PB2',
                        'Differences_A_HA_H3', 
                        'Differences_A_NA_N2', 
                        'Differences_A_H3_NP',
                        'Differences_A_H3_NS', 
                        'Differences_A_H3_PA', 
                        'Differences_A_H3_PB1', 
                        'Differences_A_H3_PB2',
                        'Differences_A_HA_H5',  
                        'Differences_A_HA_H9',  
                        'Differences_B_VIC_HA', 
                        'Differences_B_VIC_NA',
                        'ScriptTimestamp',
                        'SampleApprovalStatus',
                        'Quality']
 
    column_names_map = {
        'SequenceID': 'SequenceID',
        'RunName': 'Run Name',
        'sample': 'Sample',
        'Subtype': 'Subtype',
        'FluserverV1': 'Fluserver Mutations',
        'Clade': 'Clade',
        'thymine_ratio_823_N1': 'Neuraminidase 1 823 Thymine Ratio',
        'cytosine_ratio_823_N1': 'Neuraminidase 1 823 Cytosine Ratio',
        'Avg_Depth_A_HA_H1': 'Average Depth A HA H1',
        'Avg_Depth_A_NA_N1': 'Average Depth A NA N1',
        'Avg_Depth_A_XX_MP': 'Average Depth A MP',
        'Avg_Depth_A_H1_NP': 'Average Depth A H1 NP',
        'Avg_Depth_A_H1_NS': 'Average Depth A H1 NS',
        'Avg_Depth_A_H1_PA': 'Average Depth A H1 PA',
        'Avg_Depth_A_H1_PB1': 'Average Depth A H1 PB1',
        'Avg_Depth_A_H1_PB2': 'Average Depth A H1 PB2',
        'Avg_Depth_A_HA_H3': 'Average Depth A HA H3',
        'Avg_Depth_A_NA_N2': 'Average Depth A NA N2',
        'Avg_Depth_A_H3_NP': 'Average Depth A H3 NP',
        'Avg_Depth_A_H3_NS': 'Average Depth A H3 NS',
        'Avg_Depth_A_H3_PA': 'Average Depth A H3 PA',
        'Avg_Depth_A_H3_PB1': 'Average Depth A H3 PB1',
        'Avg_Depth_A_H3_PB2': 'Average Depth A H3 PB2',
        'Avg_Depth_A_HA_H5': 'Average Depth A HA H5',  
        'Avg_Depth_A_HA_H9': 'Average Depth A HA H9',  
        'Avg_Depth_B_VIC_HA': 'Average Depth B VIC HA',
        'Avg_Depth_B_VIC_NA': 'Average Depth B VIC NA',
        'Differences_A_HA_H1': 'Mutations A HA H1',
        'Differences_A_NA_N1': 'Mutations A NA N1',
        'Differences_A_XX_MP': 'Mutations A MP',
        'Differences_A_H1_NP': 'Mutations A H1 NP',
        'Differences_A_H1_NS': 'Mutations A H1 NS',
        'Differences_A_H1_PA': 'Mutations A H1 PA',
        'Differences_A_H1_PB1': 'Mutations A H1 PB1',
        'Differences_A_H1_PB2': 'Mutations A H1 PB2',
        'Differences_A_HA_H3': 'Mutations A HA H3',
        'Differences_A_NA_N2': 'Mutations A NA N2',
        'Differences_A_H3_NP': 'Mutations A H3 NP',
        'Differences_A_H3_NS': 'Mutations A H3 NS',
        'Differences_A_H3_PA': 'Mutations A H3 PA',
        'Differences_A_H3_PB1': 'Mutations A H3 PB1',
        'Differences_A_H3_PB2': 'Mutations A H3 PB2',
        'Differences_A_HA_H5': 'Mutations A HA H5',  
        'Differences_A_HA_H9': 'Mutations A HA H9',  
        'Differences_B_VIC_HA': 'Mutations B VIC HA',
        'Differences_B_VIC_NA': 'Mutations B VIC NA',
        'ScriptTimestamp': 'Script Timestamp',
        'SampleApprovalStatus': 'Sample Approval Status',
        'Quality' : 'Quality'
    }

    quality_columns = ['H1', 'H3', 'H5', 'H9', 'N1', 'N2', 'B']
    
    df = process_dataframe(df, runname, required_columns, column_names_map, quality_columns)
    
    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    csv_file = "/Users/rasmuskopperudriis/Downloads/app_summary.csv"
    output_file = "/Users/rasmuskopperudriis/Coding/work_folder/test/test02.csv"
    runname = "INFL01"
    main(csv_file, output_file, runname)

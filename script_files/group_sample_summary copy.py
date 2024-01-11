import pandas as pd
import sys
import numpy as np
from datetime import datetime

# Get the current date and time
start_time = datetime.now()

def add_average_depth_columns(df):
    ha_columns = [col for col in df.columns if 'Average Depth' in col and 'HA' in col]
    na_columns = [col for col in df.columns if 'Average Depth' in col and 'NA' in col]

    df['Average Depth HA'] = df[ha_columns].max(axis=1)
    df['Average Depth NA'] = df[na_columns].max(axis=1)

    return df


def add_missing_columns(df, required_columns):
    for col in required_columns:
        if not any(col in s for s in df.columns):
            df[col] = np.nan
    return df

def merge_and_rename_cols(df, search_str, new_col_name):
    cols = [col for col in df.columns if search_str in col]
    merged_col = df[cols].mode(axis=1)

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

    df.loc[df['Average'] < 50, 'Quality'] = 'Low'
    df.loc[(df['Average'] >= 50) & (df['Average'] < 100), 'Quality'] = 'Medium'
    df.loc[df['Average'] >= 100, 'Quality'] = 'High'
    
    return df

def add_mutation_columns(df):
    # Identify HA and NA mutation columns
    ha_mutation_cols = [col for col in df.columns if 'Mutations' in col and 'HA' in col]
    na_mutation_cols = [col for col in df.columns if 'Mutations' in col and 'NA' in col]

    # Function to aggregate mutations
    def aggregate_mutations(row, cols):
        mutations = row[cols].dropna()
        if len(mutations) > 1:
            return 'Mixed'
        elif len(mutations) == 1:
            return mutations.iloc[0]
        else:
            return np.nan

    # Apply the function to each row
    df['Mutations HA'] = df.apply(lambda row: aggregate_mutations(row, ha_mutation_cols), axis=1)
    df['Mutations NA'] = df.apply(lambda row: aggregate_mutations(row, na_mutation_cols), axis=1)

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
    df['ScriptVersion'] = script_version
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
                        'Avg_Depth_A_H1_HA', 
                        'Avg_Depth_A_N1_NA', 
                        'Avg_Depth_A_XX_MP', 
                        'Avg_Depth_A_H1_NP',
                        'thymine_ratio_823_N1',
                        'cytosine_ratio_823_N1',
                        'Avg_Depth_A_H1_NS', 
                        'Avg_Depth_A_H1_PA', 
                        'Avg_Depth_A_H1_PB1', 
                        'Avg_Depth_A_H1_PB2',
                        'Avg_Depth_A_H3_HA', 
                        'Avg_Depth_A_N2_NA', 
                        'Avg_Depth_A_H3_NP',
                        'Avg_Depth_A_H3_NS', 
                        'Avg_Depth_A_H3_PA', 
                        'Avg_Depth_A_H3_PB1', 
                        'Avg_Depth_A_H3_PB2',
                        'Avg_Depth_A_H5_HA',  
                        'Avg_Depth_A_H9_HA',
                        'Avg_Depth_A_H5_NA',  
                        'Avg_Depth_A_H9_NA',
                        'Avg_Depth_A_H7_HA',
                        'Avg_Depth_A_H7_NA',   
                        'Avg_Depth_B_VIC_HA', 
                        'Avg_Depth_B_VIC_NA', 
                        'Differences_A_H1_HA', 
                        'Differences_A_N1_NA', 
                        'Differences_A_XX_MP', 
                        'Differences_A_H1_NP',
                        'Differences_A_H1_NS',
                        'Differences_A_H1_PA', 
                        'Differences_A_H1_PB1', 
                        'Differences_A_H1_PB2',
                        'Differences_A_H3_HA', 
                        'Differences_A_N2_NA', 
                        'Differences_A_H3_NP',
                        'Differences_A_H3_NS', 
                        'Differences_A_H3_PA', 
                        'Differences_A_H3_PB1', 
                        'Differences_A_H3_PB2',
                        'Differences_A_H5_HA',  
                        'Differences_A_H9_HA',
                        'Differences_A_H5_NA',  
                        'Differences_A_H9_NA',
                        'Differences_A_H7_HA',
                        'Differences_A_H7_NA',   
                        'Differences_B_VIC_HA', 
                        'Differences_B_VIC_NA',
                        'ScriptTimestamp',
                        'SampleApprovalStatus',
                        'Quality',
                        'ScriptVersion']
 
    column_names_map = {
        'SequenceID': 'SequenceID',
        'RunName': 'Run Name',
        'sample': 'Sample',
        'Subtype': 'Subtype',
        'FluserverV1': 'Fluserver Mutations',
        'Clade': 'Clade',
        'Thymine_ratio_823_N1': 'Neuraminidase 1 823 Thymine Ratio',
        'Cytosine_ratio_823_N1': 'Neuraminidase 1 823 Cytosine Ratio',
        'Avg_Depth_A_H1_HA': 'Average Depth A H1 HA',
        'Avg_Depth_A_N1_NA': 'Average Depth A N1 NA',
        'Avg_Depth_A_XX_MP': 'Average Depth A MP',
        'Avg_Depth_A_H1_NP': 'Average Depth A H1 NP',
        'Avg_Depth_A_H1_NS': 'Average Depth A H1 NS',
        'Avg_Depth_A_H1_PA': 'Average Depth A H1 PA',
        'Avg_Depth_A_H1_PB1': 'Average Depth A H1 PB1',
        'Avg_Depth_A_H1_PB2': 'Average Depth A H1 PB2',
        'Avg_Depth_A_H3_HA': 'Average Depth A H3 HA',
        'Avg_Depth_A_N2_NA': 'Average Depth A N2 NA',
        'Avg_Depth_A_H3_NP': 'Average Depth A H3 NP',
        'Avg_Depth_A_H3_NS': 'Average Depth A H3 NS',
        'Avg_Depth_A_H3_PA': 'Average Depth A H3 PA',
        'Avg_Depth_A_H3_PB1': 'Average Depth A H3 PB1',
        'Avg_Depth_A_H3_PB2': 'Average Depth A H3 PB2',
        'Avg_Depth_A_H5_HA': 'Average Depth A HA H5',  
        'Avg_Depth_A_H9_HA': 'Average Depth A HA H9',
        'Avg_Depth_A_H7_HA': 'Average Depth A HA H7',
        'Avg_Depth_A_H5_NA': 'Average Depth A NA H5',  
        'Avg_Depth_A_H9_NA': 'Average Depth A NA H9',
        'Avg_Depth_A_H7_NA': 'Average Depth A NA H7',   
        'Avg_Depth_B_VIC_HA': 'Average Depth B VIC HA',
        'Avg_Depth_B_VIC_NA': 'Average Depth B VIC NA',
        'Differences_A_H1_HA': 'Mutations A H1 HA',
        'Differences_A_N1_NA': 'Mutations A N1 NA',
        'Differences_A_XX_MP': 'Mutations A MP',
        'Differences_A_H1_NP': 'Mutations A H1 NP',
        'Differences_A_H1_NS': 'Mutations A H1 NS',
        'Differences_A_H1_PA': 'Mutations A H1 PA',
        'Differences_A_H1_PB1': 'Mutations A H1 PB1',
        'Differences_A_H1_PB2': 'Mutations A H1 PB2',
        'Differences_A_H3_HA': 'Mutations A H3 HA',
        'Differences_A_N2_NA': 'Mutations A N2 NA',
        'Differences_A_H3_NP': 'Mutations A H3 NP',
        'Differences_A_H3_NS': 'Mutations A H3 NS',
        'Differences_A_H3_PA': 'Mutations A H3 PA',
        'Differences_A_H3_PB1': 'Mutations A H3 PB1',
        'Differences_A_H3_PB2': 'Mutations A H3 PB2',
        'Differences_A_H5_HA': 'Mutations A H5 HA',  
        'Differences_A_H9_HA': 'Mutations A H9 HA',
        'Differences_A_H7_HA': 'Mutations A H7 HA',
        'Differences_A_H5_NA': 'Mutations A H5 NA',  
        'Differences_A_H9_NA': 'Mutations A H9 NA',
        'Differences_A_H7_NA': 'Mutations A H7 NA',    
        'Differences_B_VIC_HA': 'Mutations B VIC HA',
        'Differences_B_VIC_NA': 'Mutations B VIC NA',
        'ScriptTimestamp': 'Script Timestamp',
        'SampleApprovalStatus': 'Sample Approval Status',
        'Quality' : 'Quality',
        'ScriptVersion' : 'Script Version'
    }

    quality_columns = ['H1', 'H3', 'H5', 'H9', 'N1', 'N2', 'B']
    
    df = process_dataframe(df, runname, required_columns, column_names_map, quality_columns)

    df = add_average_depth_columns(df)

    df = add_mutation_columns(df)


    # Add vaccine mutations to summary file
    mutation_csv = pd.read_csv(mutation_file)

    mutation_vaccine_HA = mutation_csv[mutation_csv['Ref_Name'].str.contains('HA')]
    mutation_vaccine_NA = mutation_csv[mutation_csv['Ref_Name'].str.contains('NA')]

    pivoted_df1_HA = mutation_vaccine_HA.groupby('sample')['Differences'].apply(';'.join).reset_index().rename(columns={'Differences': 'Vaccine Mutations HA'})
    pivoted_df1_NA = mutation_vaccine_NA.groupby('sample')['Differences'].apply(';'.join).reset_index().rename(columns={'Differences': 'Vaccine Mutations NA'})

    merged_df = pd.merge(pivoted_df1_HA, pivoted_df1_NA, on='sample', how='outer')
    merged_df = merged_df.rename(columns={'sample': 'Sample'})

    final_merge = pd.merge(df, merged_df, on='Sample', how='outer')


    # Add NA resistance mutation interpetations to summary file

    def calculate_h1_na_resistance(row):
        if row['Subtype'] == 'H1N1':
            return 'S110;I117;E119Q136;R152;D199;I223;S247;H275;R293;N295;I427;I436;P458;I223' if row['Fluserver Mutations'] == 'NO MATCH' else row['Fluserver Mutations']
        return ''
    
    final_merge['H1 NA Resistance'] = final_merge.apply(calculate_h1_na_resistance, axis=1)
    
    def calculate_h3_na_resistance(row):
        if row['Subtype'] == 'H3N2':
            return 'E119;Q136;I222;R224;N245;K249;E276;R292;N294;N329;S334;R371' if row['Fluserver Mutations'] == 'NO MATCH' else row['Fluserver Mutations']
        return ''
    
    final_merge['H3 NA Resistance'] = final_merge.apply(calculate_h3_na_resistance, axis=1)
    
    def calculate_bvic_na_resistance(row):
        if row['Subtype'] == 'B_Victoria':
            return 'H101;G104;E105;G108;E117;H134;H134;Q138;P139;G140;Y142;G145;N151;K152;N169;D197;A200;I221;A245;S246;G247;H273;R292;N294;K360;I361;R374;A395;L396;G407;D432;H439;H439;M464' if row['Fluserver Mutations'] == 'NO MATCH' else row['Fluserver Mutations']
        return ''

    final_merge['BVIC NA Resistance'] = final_merge.apply(calculate_bvic_na_resistance, axis=1)


    def calculate_h1_na_resistance_status(row):
        if row['Subtype'] == 'H1N1':
            return 'AANI' if row['Fluserver Mutations'] == 'NO MATCH' else 'Review'
        return ''
    
    final_merge['H1 NA Resistance Status'] = final_merge.apply(calculate_h1_na_resistance_status, axis=1)

    def calculate_h3_na_resistance_status(row):
        if row['Subtype'] == 'H3N2':
            return 'AANI' if row['Fluserver Mutations'] == 'NO MATCH' else 'Review'
        return ''
    
    final_merge['H3 NA Resistance Status'] = final_merge.apply(calculate_h3_na_resistance_status, axis=1)

    def calculate_bvic_na_resistance_status(row):
        if row['Subtype'] == 'B_Victoria':
            return 'AANI' if row['Fluserver Mutations'] == 'NO MATCH' else 'Review'
        return ''
    
    final_merge['BVIC NA Resistance Status'] = final_merge.apply(calculate_bvic_na_resistance_status, axis=1)

    def add_NA_status(row):
        values = [row['H1 NA Resistance Status'], row['H3 NA Resistance Status'], row['BVIC NA Resistance Status']]
        non_empty_values = sum(1 for value in values if pd.notna(value) and value != '')
    
        return 'Review' if non_empty_values >= 2 else next((value for value in values if pd.notna(value) and value != ''), '')

    final_merge['NA Resistance Status'] = final_merge.apply(add_NA_status, axis=1)


    def add_NA_resistance_mutation_list(row):
        values = [row['H1 NA Resistance'], row['H3 NA Resistance'], row['BVIC NA Resistance']]
        non_empty_values = sum(1 for value in values if pd.notna(value) and value != '')
    
        return 'Review' if non_empty_values >= 2 else next((value for value in values if pd.notna(value) and value != ''), '')
    
    final_merge['NA Resistance Mutations'] = final_merge.apply(add_NA_resistance_mutation_list, axis=1)


    #final_merge['H1 NA Resistance Status'] = final_merge['Fluserver Mutations'].apply(lambda x: 'AANI' if x == 'NO MATCH' else 'Review')
    #final_merge['H3 NA Resistance'] = final_merge['Fluserver Mutations'].apply(lambda x: 'E119;Q136;I222;R224;N245;K249;E276;R292;N294;N329;S334;R371' if x == 'NO MATCH' else x)
    #final_merge['BVIC NA Resistance'] = final_merge['Fluserver Mutations'].apply(lambda x: 'H101;G104;E105;G108;E117;H134;H134;Q138;P139;G140;Y142;G145;N151;K152;N169;D197;A200;I221;A245;S246;G247;H273;R292;N294;K360;I361;R374;A395;L396;G407;D432;H439;H439;M464' if x == 'NO MATCH' else x)

    # Add PA resistance mutation interpetations to summary file
    pa_mutation_csv = pd.read_csv(mutation_pa)

    mutation_vaccine_PA = pa_mutation_csv[pa_mutation_csv['Ref_Name'].str.contains('PA')]
    pivoted_df1_PA = mutation_vaccine_PA.groupby('sample')['PA_mutations'].apply(';'.join).reset_index().rename(columns={'PA_mutations': 'PA Resistance Mutations'})
    pivoted_df1_PA = pivoted_df1_PA.rename(columns={'sample': 'Sample'})

    final_merge = pd.merge(final_merge, pivoted_df1_PA, on='Sample', how='outer')

    
    def calculate_h1_pa_resistance(row):
        if row['Subtype'] == 'H1N1':
            return 'E23;K34;A36;A37;I38;119;E198;E199' if row['PA Resistance Mutations'] == 'NO MATCH' else row['PA Resistance Mutations']
        return ''
    
    final_merge['H1N1 PA Resistance'] = final_merge.apply(calculate_h1_pa_resistance, axis=1)


    def calculate_h3_pa_resistance(row):
        if row['Subtype'] == 'H3N2':
            return 'E23;K34;A36;A37;I38;119;E198;E199' if row['PA Resistance Mutations'] == 'NO MATCH' else row['PA Resistance Mutations']
        return ''
    
    final_merge['H3N2 PA Resistance'] = final_merge.apply(calculate_h3_pa_resistance, axis=1)

    print(final_merge)

    def calculate_h1_pa_resistance_status(row):
        if row['Subtype'] == 'H1N1':
            return 'AARS' if row['PA Resistance Mutations'] == 'NO MATCH' else 'Review'
        return ''
    
    final_merge['H1 PA Resistance Status'] = final_merge.apply(calculate_h1_pa_resistance_status, axis=1)

    def calculate_h3_pa_resistance_status(row):
        if row['Subtype'] == 'H3N2':
            return 'AARS' if row['PA Resistance Mutations'] == 'NO MATCH' else 'Review'
        return ''
    
    final_merge['H3 PA Resistance Status'] = final_merge.apply(calculate_h3_pa_resistance_status, axis=1)


    def add_PA_status(row):
        values = [row['H1 PA Resistance Status'], row['H3 PA Resistance Status']]
        non_empty_values = sum(1 for value in values if pd.notna(value) and value != '')
    
        return 'Review' if non_empty_values >= 2 else next((value for value in values if pd.notna(value) and value != ''), '')

    final_merge['PA Resistance Status'] = final_merge.apply(add_PA_status, axis=1)


    def add_PA_resistance_mutation_list(row):
        values = [row['H1N1 PA Resistance'], row['H3N2 PA Resistance']]
        non_empty_values = sum(1 for value in values if pd.notna(value) and value != '')
    
        return 'Review' if non_empty_values >= 2 else next((value for value in values if pd.notna(value) and value != ''), '')
    
    final_merge['PA Resistance Mutations'] = final_merge.apply(add_PA_resistance_mutation_list, axis=1)

    #final_merge['HA PA Resistance'] = final_merge['PA Resistance Mutations'].apply(lambda x: 'E23;K34;A36;A37;I38;119;E198;E199' if x == 'NO MATCH' else x)

    

    # FINAL OUTPUT

    final_merge.to_csv(output_file, index=False)
    
    

if __name__ == "__main__":
    csv_file = '/Users/rasmuskopperudriis/Coding/work_folder/runfolder/INF055_results/stat/app_summary.csv'
    output_file = '/Users/rasmuskopperudriis/Coding/work_folder/runfolder/INF055_results/result.csv'
    runname = '55'
    script_version = 'V.1'
    mutation_file = '/Users/rasmuskopperudriis/Coding/work_folder/runfolder/INF055_results/mutation/app_merged_mutation_vaccine.csv'
    mutation_pa = '/Users/rasmuskopperudriis/Coding/work_folder/runfolder/INF055_results/mutation/app_pa_mutation copy.csv'
    main(csv_file, output_file, runname)



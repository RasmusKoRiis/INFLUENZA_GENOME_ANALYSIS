import pandas as pd

# Read the CSV file
data = pd.read_csv('/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results/stat/new_influenza_pipline_short_summary.csv')

# Extract the gene name from the Ref_Name
data['gene'] = data['Ref_Name'].apply(lambda x: x.split('_')[-1])

# Calculate the total mapped reads for each gene within each sample
sample_gene_total_mapped_reads = data.groupby(['sample', 'gene'])['Mapped_Reads'].sum().reset_index()
sample_gene_total_mapped_reads.columns = ['sample', 'gene', 'total_mapped_reads']

# Merge the sample_gene_total_mapped_reads with the original data
data = data.merge(sample_gene_total_mapped_reads, on=['sample', 'gene'])

# Calculate the percentage of mapped reads per entry within each sample and gene
data['percentage'] = (data['Mapped_Reads'] / data['total_mapped_reads']) * 100

# Filter the data to only include entries with percentages under 1
data_under_1 = data[data['percentage'] < 1]

# Output the result as a CSV
data_under_1[['sample', 'Ref_Name', 'gene', 'percentage']].to_csv('/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results/stat/output_1.csv', index=False)

print("Percentages under 1:")
print(data_under_1[['sample', 'Ref_Name', 'gene', 'percentage']])

# Filter the data to only include entries with percentages under 5 and greater than or equal to 1
data_under_5 = data[(data['percentage'] >= 5) & (data['percentage'] <= 85)]

# Output the result as a CSV
data_under_5[['sample', 'Ref_Name', 'gene', 'percentage']].to_csv('/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results/stat/output_5.csv', index=False)

print("\nPercentages under 5 and greater than or equal to 1:")
print(data_under_5[['sample', 'Ref_Name', 'gene', 'percentage']])

average_percentage = data['percentage'].mean()
print("Average percentage:", average_percentage)





import pysam
import pandas as pd
import numpy as np
import sys

bam_file = sys.argv[1]
stats_file = sys.argv[2]
output_file_sample = sys.argv[3]
output_file_summary = sys.argv[4]

# Open the BAM file and extract the read IDs
bam_file = pysam.AlignmentFile(bam_file, 'rb')

# Read the statistics file into a DataFrame
stats_file = pd.read_csv(stats_file, sep='\t')

# Add a new column for the reference
stats_file['reference'] = ''

# Loop through the reads in the BAM file and update the reference column with the reference name
for read in bam_file.fetch():
    if not read.is_unmapped:
        ref_name = read.reference_name
        read_id = read.query_name
        stats_file.loc[stats_file['read_id'] == read_id, 'reference'] = ref_name
# Close the BAM file
bam_file.close()

# Calculate the mean quality and read count for each reference
quality_summary = stats_file.groupby(['reference']).agg({'mean_quality': np.mean, 'sample_name': 'first', 'read_id': 'count'})

# Rename the reads column to count_reads
quality_summary = quality_summary.rename(columns={'read_id': 'count_reads'})

# Reset the index to turn the 'reference' column into a regular column
quality_summary = quality_summary.reset_index()

# Fill empty cells in the reference column with "not_mapped"
quality_summary['reference'] = quality_summary['reference'].fillna('not_mapped')

# Print the result
print(quality_summary)

# Save the updated statistics file
stats_file.to_csv(output_file_sample, sep='\t')
quality_summary.to_csv(output_file_summary, sep='\t')

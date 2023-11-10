import pysam
import pandas as pd
import numpy as np
import sys

bam_file_path = sys.argv[1]
stats_file_path = sys.argv[2]
output_file_sample = sys.argv[3]
output_file_summary = sys.argv[4]

# Open the BAM file
bam_file = pysam.AlignmentFile(bam_file_path, 'rb')

# Read the statistics file into a DataFrame
stats_file = pd.read_csv(stats_file_path, sep='\t')

# Initialize a dictionary to store reference names for each read ID
reference_dict = {}

# Iterate through the reads in the BAM file and store the reference names in the dictionary
for read in bam_file:
    if not read.is_unmapped:
        reference_dict[read.query_name] = read.reference_name

# Close the BAM file
bam_file.close()

# Map the reference names to the DataFrame using the dictionary
stats_file['reference'] = stats_file['read_id'].map(reference_dict).fillna('not_mapped')

# Calculate the mean quality and read count for each reference
quality_summary = stats_file.groupby('reference').agg({'mean_quality': np.mean, 'sample_name': 'first', 'read_id': 'count'}).rename(columns={'read_id': 'count_reads'}).reset_index()

# Save the updated statistics file and the summary
stats_file.to_csv(output_file_sample, sep='\t', index=False)
quality_summary.to_csv(output_file_summary, sep='\t', index=False)

print("Processing complete. Files saved.")
import csv
import sys
import os.path

csv_file = sys.argv[1]
output_file = sys.argv[2]

# Extract the sample name from the input file name
sample_name = os.path.basename(csv_file).split('.')[0:2]

# Join the sample name back together with an underscore
sample_name = '_'.join(sample_name)

# Add the Sample column to the input file and write the output to the output file
with open(csv_file, newline='') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile, delimiter='\t')
    writer = csv.writer(outfile, delimiter='\t')

    # Add a header row to the output file
    header = next(reader)
    header.append('Sample')
    writer.writerow(header)

    # Iterate over the input rows and write the new column to the output file
    for row in reader:
        row.append(sample_name)
        writer.writerow(row)
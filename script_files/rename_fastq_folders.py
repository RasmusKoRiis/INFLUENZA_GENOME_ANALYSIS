import csv
import os
import sys

fastq_path = sys.argv[1]
csv_path = sys.argv[2]

# path of the directory
directory_path = fastq_path
# path of the csv file
csv_file_path = csv_path

# read csv file as a list of dictionaries
with open(csv_file_path, 'r') as f:
    reader = csv.DictReader(f)
    rows = list(reader)

# create a mapping of barcodes to sampleids
barcode_to_sampleid = {row['barcode']: row['sampleid'] for row in rows}

# iterate over all directories in the directory
for dir_name in os.listdir(directory_path):
    if os.path.isdir(os.path.join(directory_path, dir_name)):
        # if the directory name (barcode) is in our mapping, rename it to the sampleid
        if dir_name in barcode_to_sampleid:
            new_name = barcode_to_sampleid[dir_name]
            os.rename(os.path.join(directory_path, dir_name), os.path.join(directory_path, new_name))

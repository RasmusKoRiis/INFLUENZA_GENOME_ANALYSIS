import pandas as pd
import sys
import numpy as np
from datetime import datetime

fasta_file = sys.argv[1]
output_file = sys.argv[2]

# Define the mapping of two-letter codes to numbers
code_to_number = {
    "HA": "01",
    "NA": "02",
    "MP": "03",
    "PB1": "04",
    "PB2": "05",
    "NP": "06",
    "PA": "07",
    "NS": "08"
}

# Function to convert the header without spaces
def convert_header_no_spaces(header):
    # Splitting the header to get the ID and other components
    parts = header.split("_")
    id_part = parts[0]
    region = parts[-1]  # The last part e.g. HA, NA, etc.
    
    # Mapping the region to its corresponding number
    region_number = code_to_number.get(region, "")
    
    # Constructing the new header without spaces
    new_header = f"{id_part}|{region_number}-{region}|{parts[-2]}"
    
    return new_header

# Reading the input FASTA file and converting the headers
input_file_path = fasta_file
with open(input_file_path, "r") as f:
    lines = f.readlines()

new_lines_no_spaces = []
for line in lines:
    if line.startswith(">"):
        new_lines_no_spaces.append(convert_header_no_spaces(line.strip()) + "\n")
    else:
        new_lines_no_spaces.append(line)

# Writing the sequences with new headers without spaces to a new file
output_file_no_spaces_path = output_file
with open(output_file_no_spaces_path, "w") as f:
    f.writelines(new_lines_no_spaces)

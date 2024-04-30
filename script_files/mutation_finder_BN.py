import pandas as pd
import sys
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


sequence_file = sys.argv[1]
reference_file = sys.argv[2]
segment = sys.argv[3]
output_file = sys.argv[4]
print("python segment: {}".format(segment))

def find_differences(reference, seq):
    # Align the sequences
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -2
    alignments = aligner.align(reference, seq)
    
    # Get the alignment with the highest score
    #best_alignment = max(alignments, key=lambda x: x.score)
    best_alignment = max(alignments)

    
    # Find the differences between the sequences
    differences = []
    #for i in range(len(best_alignment)-1):
    for i in range(len(best_alignment[1])):
        if best_alignment[0][i] != best_alignment[1][i]:
            differences.append(f"{best_alignment[0][i]}{i+1}{best_alignment[1][i]}")
    return ';'.join(differences)


def check_frameshift(seq):
    if 'X' in seq:
        return seq.index('X') + 1
    else:
        return 'N.A.'

def process_differences(row):
    if 'FRAMESHIFT' in str(row['Frameshift/Poor Seq']).upper().strip():
       return row['Frameshift/Poor Seq']
    else:
        return row['Differences']
    
    #xs_count = sum([1 for x in row['Differences'].split(';') if 'X' in x])
    #if xs_count > 20:
    #    return 'Pos. FS ' + str(row['Frameshift/Poor Seq'])
    #else:
    #    return row['Differences']

def find_frameshift(reference, seq):
    # Align the sequences
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    alignments = aligner.align(reference, seq)
    
    # Get the alignment with the highest score
    best_alignment = max(alignments)
    
    ref_aligned = best_alignment[0]
    seq_aligned = best_alignment[1]
    
        # Find frameshift
    consecutive_divergences = 0
    frameshift_detected = False
    start_divergence_position = None
    for i, (ref_aa, seq_aa) in enumerate(zip(ref_aligned, seq_aligned)):
        if ref_aa != seq_aa:
            if consecutive_divergences == 0:
                start_divergence_position = i
            consecutive_divergences += 1
            if consecutive_divergences > 5:
                frameshift_detected = True
                break
        else:
            consecutive_divergences = 0
    
    if frameshift_detected:
        return 'FRAMESHIFT POS: ' + str(int(start_divergence_position) + 1)
    else:
        return None


sequences = []

for ref in SeqIO.parse(reference_file, 'fasta'):
    reference = Seq(str(ref.seq))
    for record in SeqIO.parse(sequence_file, 'fasta'):
        sequence = Seq(str(record.seq))
        differences = find_differences(reference, sequence)
        frameshift = check_frameshift(str(sequence))
        frameshift2 = find_frameshift(reference, sequence)
        sequences.append({'ID': record.id, 'Differences': differences, 'Frameshift/Poor Seq': frameshift2})

# Create a DataFrame from the sequences
df = pd.DataFrame(sequences)
print(df)
    
df['Differences'] = df.apply(process_differences, axis=1)

# Split "ID" column into "sample" and "Ref_Name" columns
df[['Ref_Name', 'sample']] = df['ID'].str.split('|', n=1, expand=True)
print(df['Ref_Name'])

# Set the "Ref_Name" column to the 'segment' variable value
#df['Ref_Name'] = segment

# Combine the last character of the "sample" column with the first character of the "Ref_Name" column
#df['sample'] = df['sample'] + '_' + df['Ref_Name'].str[:1]


# Remove the first character of the "Ref_Name" column
df['Ref_Name'] = df['Ref_Name'].str[2:]
df['Ref_Name'] = segment

# Drop the original "ID" column
df.drop(columns=['ID'], inplace=True)

# Reorder the columns
df = df[['sample', 'Ref_Name', 'Differences', 'Frameshift/Poor Seq']]

# Save the final dataframe to a CSV file
df.to_csv(output_file, index=False)

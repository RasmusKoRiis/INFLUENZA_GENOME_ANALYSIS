import sys
from Bio import SeqIO

input_file = sys.argv[1]
output_file = sys.argv[2]
print(output_file)


with open(output_file, 'w') as output_handle:
    for record in SeqIO.parse(input_file, 'fasta'):
        seq = str(record.seq)
        n_count = seq.count('N') + seq.count('n')
        seq_length = len(seq)
        if n_count/seq_length < 0.15:
            SeqIO.write(record, output_handle, 'fasta')
           
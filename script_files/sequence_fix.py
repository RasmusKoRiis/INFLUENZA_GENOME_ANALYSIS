import sys
from Bio import SeqIO

input_file = sys.argv[1]
output_file = sys.argv[2]

def remove_dashes(fasta_in, fasta_out):
    sequences = list(SeqIO.parse(fasta_in, "fasta"))
    
    for seq in sequences:
        seq.seq = seq.seq.ungap("-")
    
    print(fasta_out)
        
    SeqIO.write(sequences, fasta_out, "fasta")

# Call the function with your input and output file paths
remove_dashes(input_file, output_file)
S
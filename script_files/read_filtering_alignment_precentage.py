import pysam
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

def proportion_aligned(cigar_tuples):
    """Calculate the proportion of the read that's aligned based on the CIGAR string."""
    if cigar_tuples is None:  # Handle unaligned reads
        return 0
    
    aligned_bases = sum(length for op, length in cigar_tuples if op in (0, 1, 2))  # M=0, I=1, D=2
    total_bases = sum(length for op, length in cigar_tuples if op != 5)  # H=5 (hard clip) is excluded
    return aligned_bases / total_bases

def filter_reads(input_file, output_file, threshold=0.85):
    """Filter reads from an input BAM file and write to an output BAM file."""
    with pysam.AlignmentFile(input_file, "rb") as infile, pysam.AlignmentFile(output_file, "wb", header=infile.header) as outfile:
        for read in infile:
            if proportion_aligned(read.cigartuples) >= threshold:
                outfile.write(read)

filter_reads(input_file, output_file)

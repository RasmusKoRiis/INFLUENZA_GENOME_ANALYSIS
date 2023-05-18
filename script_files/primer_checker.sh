#!/bin/bash

# PATHWAYS
startdir=$(pwd)
runname=$(basename $startdir)
script="$startdir/script_files"

result_folder="$startdir/results"
fasta_folder="$result_folder/fasta"
bam_folder="$result_folder/bam"
stat_folder="$result_folder/stat"
mutation_folder="$result_folder/mutation"
dataset_folder="$startdir/dataset"

arr=("A_H3_MP" "A_H1_MP" )


for i in "${arr[@]}"
do
    echo $i
    blastn -task blastn-short \
        -gapopen 3 \
        -penalty -1 \
        -outfmt 6 \
        -query $dataset_folder/primer_database_with_duplicate_names_v1.fasta \
        -subject $fasta_folder/${i}.fasta > $mutation_folder/blastn_${i}.csv


    python3 ${script}/primer_aligner.py $mutation_folder/blastn_${i}.csv $i "$mutation_folder"/"$i"

done


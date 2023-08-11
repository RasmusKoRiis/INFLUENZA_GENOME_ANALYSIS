#!/bin/bash

#Pipline for generation of Influenza A + B whole genome sequences, 
#including analysis producing mutation calling, vaccine-contamination protocols and more....


#TECHNICAL INFO OF PIPLINE
date=$(date +"%Y-%m-%d_%H-%M-%S")
startdir=$(pwd)
script="$startdir/script_files"
reference="$startdir/references/epi2me/reference_epi2me_FULL_NAMES.fasta"
fastq_in="$startdir/input_fastq_processed"

cd $startdir

mkdir output_fastq

mkdir output_fastq_cleaned

fastq_out="$startdir/output_fastq"
fastq_out_cleaned="$startdir/output_fastq_cleaned"

cd $fastq_in

for subdir in */; do
    subdir_name=$(basename $subdir)
    out_file_name="$subdir_name.fastq"
    cat $subdir/*.fastq > $out_file_name
    mv $out_file_name $cur_dir/$fastq_out/
done

cd $startdir

#PLANNED NEXTFLOW IMPLEMENTATION
NXF_VER=22.10.6 nextflow run "script_files/mainQA.nf" \
    --script_files "$script"\
    --in_fastq "$fastq_out" \
    --out_fastq "$fastq_out_cleaned"\
    --reference "$reference"\





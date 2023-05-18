#!/bin/bash

#Pipline for generation of Influenza A + B whole genome sequences, 
#including analysis producing mutation calling, vaccine-contamination protocols and more....


#TECHNICAL INFO OF PIPLINE
date=$(date +"%Y-%m-%d_%H-%M-%S")
startdir=$(pwd)
script="/app/script_files"
runname=$(basename $startdir)

cd $startdir

result_folder="/app/results"
fasta_folder="$result_folder/fasta"
bam_folder="$result_folder/bam"
stat_folder="$result_folder/stat"
mutation_folder="$result_folder/mutation"
dataset_folder="/app/dataset"
reference="/app/references"


#PLANNED NEXTFLOW IMPLEMENTATION
cd $startdir

NXF_VER=22.10.6 nextflow run "script_files/main.nf" \
    --script_files "$script"\
    --in_fasta "$result_folder/epi2me_wf_flu_output" \
    --out_fasta "$fasta_folder"\
    --in_bam "$result_folder/epi2me_wf_flu_output" \
    --out_bam "$bam_folder"\
    --in_stat "$result_folder/epi2me_wf_flu_output" \
    --out_stat "$stat_folder"\
    --out_mutation "$mutation_folder"\
    --in_dataset "$dataset_folder"\
    --reference "$reference"\
    --runname "$runname"\
    --results "$result_folder"\


#FINALIZIN§G SUMMARY
#python3 "$script/group_sample_summary.py" "$stat_folder/"$runname"_long_quality_mutation_summary.csv" "$result_folder"/"$runname"_final.csv


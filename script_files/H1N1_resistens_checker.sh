#!/bin/bash

# PATHWAYS
startdir=$(pwd)
runname=$(basename $startdir)
script="$startdir/script_files"

result_folder="$startdir/results"
fasta_folder="$result_folder/fasta"
bam_folder="$result_folder/bam"
stat_folder="$result_folder/stat"
dataset_folder="$startdir/dataset"
mutation_folder="$result_folder/mutation"

# MAKE RAW DEPTH FILES
for bam_file in "$bam_folder"/*.bam; do
    bam_name=$(basename $bam_file .bam)
    bam-readcount -w1 $bam_file > "$mutation_folder"/"$bam_name"_position_count_raw.tsv
done

# PROCESS RAW DEPTH FILES
for depth_file in "$mutation_folder"/*raw.tsv; do
    sample_name=$(basename $depth_file _position_count_raw.tsv)
    python3 "$script"/mutation_ratio_finder.py $depth_file "$mutation_folder"/"$sample_name"_processed_raw.csv $sample_name 823
done

# MERGE ALL DEPTH TABELS
python3 "$script/table_merger.py" "processed_raw.csv" "$mutation_folder/"$runname"_merged_depth_raw.csv" "$mutation_folder"

# APPEND TO LONG SUMMARY
python3 "$script/column_lookup_append.py"  \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \
        "$mutation_folder/"$runname"_merged_depth_raw.csv" \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "cytosine_ratio__823" \

# APPEND TO LONG SUMMARY
python3 "$script/column_lookup_append.py"  \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \
        "$mutation_folder/"$runname"_merged_depth_raw.csv" \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "thymine_ratio__823" \

# ADD RATIO FIX - REMOVING RATIO FROM SAMPLES WITH DEPT BELLOW 20

# CLEAN UP
#rm "$mutation_folder"/*raw.csv "$mutation_folder"/*raw.tsv
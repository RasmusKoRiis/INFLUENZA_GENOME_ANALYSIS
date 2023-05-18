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

# MOVE FASTA FILES
mv "$result_folder/epi2me_wf_flu_output"/*.fasta "$fasta_folder"

# REMOVE SEQUNECES WITHS ONLY Ns
for fasta_file in "$fasta_folder"/*.fasta; do
    output_file="${fasta_file}_filtered.fasta"
    python3 "$script/failed_sequence_removal.py" "$fasta_file" "$output_file"
    rm $fasta_file
done

# ADD SAMPLE NAME INSIDE FASTA
for fasta_file in "$fasta_folder"/*.fasta; do
    file_name=$(basename "$fasta_file" .fasta | cut -d '.' -f 1)
    part1=$(echo "$file_name" | cut -d '_' -f 1)
    part2=$(echo "$file_name" | cut -d '_' -f 2-)
    sed "s/>\(.*\)/>${part1}_${part2}_\1/" "$fasta_file" > "$fasta_file.tmp"
    mv "$fasta_file.tmp" "$fasta_file"
done

# MERGING FASTA FILES * fix merged.fasta files name
cat "$fasta_folder"/*.fasta > "$fasta_folder/"$runname"_merged.fasta"

# MOVE INDIVUAL FASTA FILES
mkdir "$fasta_folder/single_fasta_files"
mv "$fasta_folder"/*filtered.fasta "$fasta_folder/single_fasta_files"

# EXTRACT SEGMENTS INTO FILES
match_strings=("A_H1_NS" "A_H3_NS" "A_H3_MP" "A_H1_MP" "A_H3_NP" "A_H1_NP" "A_H3_PA" "A_H1_PA" "A_H3_PB1" "A_H1_PB1" "A_H3_PB2" "A_H1_PB2" "A_HA_H1" "A_HA_H10" "A_HA_H11" "A_HA_H13" "A_HA_H14" "A_HA_H15" "A_HA_H16" "A_HA_H2" "A_HA_H3" "A_HA_H4" "A_HA_H5" "A_HA_H6" "A_HA_H7" "A_HA_H8" "A_HA_H9" "A_NA_N1" "A_NA_N2" "A_NA_N3" "A_NA_N5" "A_NA_N6" "A_NA_N7" "A_NA_N8" "A_NA_N9" "B_MP" "B_VIC_HA" "B_VIC_NA")

for match_string in "${match_strings[@]}"; do
    seqkit grep -r -p "${match_string}$" "$fasta_folder"/*.fasta > "$fasta_folder/${match_string}.fasta"
done

# REMOVE EMPTY FASTA FILES
find "$fasta_folder" -name "*.fasta" -type f -size 0 | xargs rm


# MOVE BAM AND BAI FILES
mv "$result_folder/epi2me_wf_flu_output"/*.bam "$bam_folder"
mv "$result_folder/epi2me_wf_flu_output"/*.bai "$bam_folder"

# MOVE STAT FILES
mv "$result_folder/epi2me_wf_flu_output"/*.{stats,txt,csv} "$stat_folder"

# ADD SAMPLE NAME TO DEPTH FILE AND MERGE
for file in "$stat_folder"/*.depth.txt; do
    filename=$(basename "$file" .depth.txt)
    awk -v OFS="\t" -v filename="$filename" '{print filename, $0}' "$file" >> "$stat_folder"/"$runname"_heatmap_depth.csv
done

# GENERATE STAT SUMMARY FILE
for bam_file in "$bam_folder"/*.bam; do
    sample_name=$(basename $bam_file) 
    weeSAM --bam "$bam_file" --out "$stat_folder/$sample_name.csv"
done

# ADD SAMPLE COLUMN TO TABLE
for csv_file in "$stat_folder"/*.bam.csv; do
    output_file="${csv_file}_fixed.csv"
    python3 "$script/add_sample_column.py" "$csv_file" "$output_file"
done

# MERGE SUMMARY FILEs
cat "$stat_folder"/*fixed.csv | awk '!/^Ref_Name/' > "$stat_folder/"$runname"_merged_summary_temp.csv"
echo -e "Ref_Name\tRef_Len\tMapped_Reads\tBreadth\t%_Covered\tMin_Depth\tMax_Depth\tAvg_Depth\tStd_Dev\tAbove_0.2_Depth\tAbove_1_Depth\tAbove_1.8_Depth\tVariation_Coefficient\tSample" > "$stat_folder/header.csv"
cat "$stat_folder/header.csv" "$stat_folder/"$runname"_merged_summary_temp.csv" > "$stat_folder/"$runname"_merged_summary.csv"

# STAT CLEAN UP
rm "$stat_folder"/*bam.csv "$stat_folder/"$runname"_merged_summary_temp.csv" "$stat_folder"/*fixed.csv "$stat_folder"/header.csv "$stat_folder"/*depth.txt

# MAKE FINAL SUMMARIES
csv_file="$stat_folder/wf-flu-results.csv"
tab_file="$stat_folder/${runname}_merged_summary.csv"
output_file_long="$stat_folder/${runname}_long_summary.csv"
output_file_short="$stat_folder/${runname}_short_summary.csv"
python3 "$script/add_subtype_to_summary.py" "$csv_file" "$tab_file" "$output_file_long" "$output_file_short"

for bam_file in "$bam_folder"/*.bam; do
    sample_name=$(basename $bam_file .bam) 
    stats_file="$stat_folder/${sample_name}.stats"
    output_file_sample="$stat_folder/${sample_name}_processed.stats"
    output_file_summary="$stat_folder/${sample_name}_read_quality.csv"
    python3 "$script/read_quality_summaries.py" "$bam_file" "$stats_file" "$output_file_sample" "$output_file_summary"
done

# MERGE STATS SUMMARIES - EXT SAVE_NAME - FOLDER
python3 "$script/table_merger.py" "processed.stats" "$stat_folder/"$runname"_processed_summary.stats" "$stat_folder"
python3 "$script/table_merger.py" "read_quality.csv" "$stat_folder/"$runname"_read_quality_summary.csv" "$stat_folder"

# STAT CLEAN UP 2
rm "$stat_folder"/*_processed.stats "$stat_folder"/*_read_quality.csv "$stat_folder"/*merged_summary.csv "$stat_folder"/wf-flu-results.csv 

# ADD FRAGMENT QUALITY TO LONG SUMMARY
python3 "$script/column_lookup_append.py"  \
        "$stat_folder/${runname}_long_summary.csv" \
        "$stat_folder/${runname}_read_quality_summary.csv" \
        "$stat_folder/${runname}_long_quality_summary.csv" \
        "," \
        "\t" \
        "sample,Ref_Name" \
        "sample_name,reference" \
        "mean_quality" \

# STAT CLEAN UP 3
rm "$stat_folder"/*_long_summary.csv

# B-SUBTYPE FIX
python3 "$script/B_subtype_fix.py" "$stat_folder/${runname}_long_quality_summary.csv"  "$stat_folder/${runname}_long_quality_summary.csv"  


mv "$fasta_folder"/*merged.fasta "$fasta_folder/single_fasta_files"


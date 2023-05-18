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

# TRANSLATE FASTA FILES TO AMINOACID 

for fasta_file in "$fasta_folder"/*.fasta; do
    segment_name=$(basename $fasta_file .fasta)
    seqkit grep -r -p "^"$segment_name"$" $startdir/references/epi2me/reference_epi2me_FULL_NAMES.fasta > $mutation_folder/"$segment_name"_temp.fasta

    nextalign run \
        --input-ref=$mutation_folder/"$segment_name"_temp.fasta \
        --genemap=$startdir/references/epi2me/singel_files/"$segment_name"_genemap.gff \
        --output-all=$mutation_folder/$segment_name/ \
        $fasta_folder/"$segment_name".fasta

    mv "$mutation_folder"/"$segment_name"/*translation.fasta "$mutation_folder"
    rm -rf "$mutation_folder"/"$segment_name"
    rm "$mutation_folder"/*_temp.fasta 

    # Move the file without using segment_name in the last line
    for translation_file in "$mutation_folder"/nextalign_gene_*.translation.fasta; do
        new_translation_filename=$(echo "$(basename "$translation_file")" | sed 's/nextalign_gene_//')
        mv "$translation_file" "$mutation_folder"/"$new_translation_filename"
    done
done

# MUTATION LIST
for fasta_file in "$mutation_folder"/*.fasta; do
    fasta_name=$(basename $fasta_file .translation.fasta)
    segment=$(basename $fasta_file .translation.fasta)
    echo "this is the: $segment"
    reference="$startdir/references/epi2me/singel_files/${segment}_amino.fasta"
    output_file="$mutation_folder/${segment}.csv"
    python3 "$script/mutation_finder.py" $fasta_file $reference $segment $output_file
done

# MERGE ALL MUTATION TABELS
python3 "$script/table_merger.py" "csv" "$mutation_folder/"$runname"_merged_mutation.csv" "$mutation_folder"


# APPEND TO LONG SUMMARY
python3 "$script/column_lookup_append.py"  \
        "$stat_folder/${runname}_long_quality_summary.csv" \
        "$mutation_folder/"$runname"_merged_mutation.csv" \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "Differences" \

#ADD HA2 TO SUMMARY
python3 "$script/HA2_to_summary.py" \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \
        "$mutation_folder/"$runname"_merged_mutation.csv" \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \

# ADD FLUREVSER MUTATION MATCHES
python3 "$script/mutation_annotation.py" "$mutation_folder/"$runname"_merged_mutation.csv" \
        "$startdir/dataset/RESITENCE_MUTATION/NA_FLUSERVER.csv" \
        "$mutation_folder/"$runname"_merged_mutation.csv"  \


# APPEND TO LONG SUMMARY
python3 "$script/column_lookup_append.py"  \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \
        "$mutation_folder/"$runname"_merged_mutation.csv" \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "Matches" \

#ADD LOW HIG MEDIUM NA

#ADD LOW HIG MEDIUM PA

# DOWNLOAD DATASET?

# FIND CLADE
if [ -f "$fasta_folder/A_HA_H1.fasta" ]; then
    nextclade run \
      --input-dataset $dataset_folder/HA_NEXTCLADE/A_HA_H1 \
      --output-all=$mutation_folder/H1_HA_CLADE \
      $fasta_folder/A_HA_H1.fasta
    
    mv "$mutation_folder"/H1_HA_CLADE/nextclade.csv "$mutation_folder"/A_H1_HA_CLADE.csv
    python3 "$script/clade_table_fix.py" "$mutation_folder"/A_H1_HA_CLADE.csv "$mutation_folder"/A_H1_HA_CLADE.csv
    rm -rf "$mutation_folder"/H1_HA_CLADE

else
    echo "A_HA_H1.fasta file does not exist"
fi

if [ -f "$fasta_folder/A_HA_H3.fasta" ]; then
    nextclade run \
      --input-dataset $dataset_folder/HA_NEXTCLADE/A_HA_H3 \
      --output-all=$mutation_folder/H3_HA_CLADE \
      $fasta_folder/A_HA_H3.fasta

    mv "$mutation_folder"/H3_HA_CLADE/nextclade.csv "$mutation_folder"/A_H3_HA_CLADE.csv
    python3 "$script/clade_table_fix.py" "$mutation_folder"/A_H3_HA_CLADE.csv "$mutation_folder"/A_H3_HA_CLADE.csv
    rm -rf "$mutation_folder"/H3_HA_CLADE
else
    echo "A_HA_H3.fasta file does not exist"
fi

if [ -f "$fasta_folder/B_VIC_HA.fasta" ]; then
    nextclade run \
      --input-dataset $dataset_folder/HA_NEXTCLADE/B_VIC \
      --output-all=$mutation_folder/B_HA_CLADE \
      $fasta_folder/B_VIC_HA.fasta


    mv "$mutation_folder"/B_HA_CLADE/nextclade.csv "$mutation_folder"/B_HA_CLADE.csv
    python3 "$script/clade_table_fix.py" "$mutation_folder"/B_HA_CLADE.csv "$mutation_folder"/B_HA_CLADE.csv
    rm -rf "$mutation_folder"/B_HA_CLADE
else
    echo "B_VIC_HA.fasta file does not exist"
fi

# ADD CLADE TO SUMMARY
python3 "$script/table_merger.py" "CLADE.csv" "$mutation_folder/"$runname"_merged_clade.csv" "$mutation_folder"

python3 "$script/column_lookup_append.py"  \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \
        "$mutation_folder/"$runname"_merged_clade.csv" \
        "$stat_folder/${runname}_long_quality_mutation_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "clade" \













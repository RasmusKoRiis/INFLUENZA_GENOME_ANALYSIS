#!/usr/bin/env nextflow

params.in_fasta = ""
params.out_fasta = ""
params.in_bam = ""
params.out_bam = ""
params.in_stat = ""
params.out_stat = ""
params.out_mutation = ""
params.in_dataset = ""
params.script_files = ""
params.runname = ""
params.runname2 = ""
params.script_version = ""
params.reference = ""



log.info """\
    I N  H O U S E - N F   P I P E L I N E
    ===================================
    fasta: ${params.in_fasta}
    script: ${params.script_files}
    fasta_out: ${params.out_fasta}
    mutation_out: ${params.out_mutation}
    runname: ${params.runname2}
    """.stripIndent()

Channel
    .fromPath(params.in_fasta + "/*.fasta")
    .set { fasta_files_ch }
Channel
    .fromPath(params.in_bam + "/*.bam")
    .set { bam_files_ch }
Channel
    .fromPath(params.in_bam + "/*.bai")
    .set { bai_files_ch }
Channel
    .fromPath(params.in_stat + "/*.{stats,txt,csv}")
    .set { stat_files_ch }


process REMOVE_HIGH_N_SAMPLES {

    input:
    path in_fasta from fasta_files_ch

    output:
    path "n_removed_${in_fasta.name}" into n_removed_fasta_files_ch

    script:
    """
    python ${params.script_files}/failed_sequence_removal.py ${in_fasta} n_removed_${in_fasta.name}
    """
}

process CORRECT_FASTA {

    input:
    path in_fasta from n_removed_fasta_files_ch

    output:
    path "${in_fasta.name}" into coreccted_fasta_files_ch

    script:
    """
    python ${params.script_files}/sequence_fix.py ${in_fasta} ${in_fasta.name}
    """
}

process ADD_NAME_INSIDE_FASTA {

    input:
    path in_fasta from coreccted_fasta_files_ch

    output:
    path "names_added_${in_fasta.name}" into name_added_fasta_files_ch, name_added_fasta_files_ch2

    script:
    """
    file_name=\$(basename ${in_fasta} .fasta | sed 's/^n_removed_//' | cut -d '.' -f 1)
    part1=\$(echo "\$file_name" | cut -d '_' -f 1)
    part2=\$(echo "\$file_name" | cut -d '_' -f 2-)
    sed "s/>\\(.*\\)/>\${part1}_\${part2}_\\1/" ${in_fasta} > "${in_fasta}.tmp"
    mv "${in_fasta}.tmp" "names_added_${in_fasta.name}"
    """
}

process SAVE_FASTA_FILES {
    publishDir params.out_fasta, mode: 'copy'

    input:
    path in_fasta from name_added_fasta_files_ch.toList()

    output:
    path "single_fasta_files" into singel_fasta_ch, singel_fasta_ch2

    script:
    """
    mkdir -p single_fasta_files

    # MOVE INDIVUAL FASTA FILES
    for file in *names_added_*.fasta; do
        new_name=\$(echo "\$file" | sed 's/^names_added_n_removed_//')
        mv "\$file" single_fasta_files/"\$new_name"
    done
    """
}

process MERGE_AND_EXTRACT_FASTA {
    publishDir params.out_fasta, mode: 'copy'

    input:
    path in_fasta from singel_fasta_ch

    output:
    path "*.fasta" into merged_and_extracted_ch, merged_and_extracted_ch2, merged_and_extracted_ch3, merged_and_extracted_ch4
    path "*.fa" into merged_fasta_ch

    script:
    """
    runname=${params.runname2}

    # MERGING FASTA FILES
    cat single_fasta_files/*.fasta > "\${runname}_merged.fa"
    

    # EXTRACT SEGMENTS INTO FILES
    match_strings=("A_H1_NS" "A_H3_NS" "A_XX_MP" "A_H1_MP" "A_H3_NP" "A_H1_NP" "A_H3_PA" "A_H1_PA" "A_H3_PB1" "A_H1_PB1" "A_H3_PB2" "A_H1_PB2" "A_H1_HA" "A_H10_HA" "A_H11_HA" "A_H13_HA" "A_H14_HA" "A_H15_HA" "A_H16_HA" "A_H2_HA" "A_H3_HA" "A_H4_HA" "A_H5_HA" "A_H6_HA" "A_H7_HA" "A_H8_HA" "A_H9_HA" "A_N1_NA" "A_N2_NA" "A_N3_NA" "A_N5_NA" "A_N6_NA" "A_N7_NA" "A_N8_NA" "A_N9_NA" "B_MP" "B_VIC_HA" "B_VIC_NA")

    for match_string in "\${match_strings[@]}"; do
        seqkit grep -r -p "\${match_string}\$" "\${runname}_merged.fa" > "\${match_string}.fasta"
    done

    # REMOVE EMPTY FASTA FILES
    find . -name "*.fasta" -type f -size 0 -delete
    """
}

process SAVE_BAM_FILES {

    publishDir "${params.out_bam}", mode: 'copy'

    input:
    path in_bam from bam_files_ch

    output:
    path "${in_bam.name}" into moved_bam_ch, moved_bam_ch2

    script:
    """
    cp ${in_bam} ${params.out_bam}/${in_bam.name}
    """
}

process SAVE_BAI_FILES {

    publishDir "${params.out_bam}", mode: 'copy'

    input:
    path in_bai from bai_files_ch

    output:
    path "${in_bai.name}" into moved_bai_ch, moved_bai_ch2

    script:
    """
    cp ${in_bai} ${params.out_bam}/${in_bai.name}
    """
}

process SAVE_STAT_FILES {
    publishDir "${params.out_stat}/", mode: 'copy'

    input:
    path in_stat from stat_files_ch

    output:
    path "${in_stat.name}" into moved_stat_ch, moved_stat_ch2, moved_stat_ch3

    script:
    """
    cp ${in_stat} ${params.out_stat}/${in_stat.name}
    """
}

process ADD_SAMPLE_NAME_AND_MERGE_STAT {
    publishDir params.out_stat, mode: 'copy'

    input:
    path in_depth from moved_stat_ch.filter { it.name.endsWith('.depth.txt') }.collect()

    output:
    path("merged_heatmap_depth.csv") into merged_depth_ch

    script:
    """
    # ADD SAMPLE NAME TO DEPTH FILE AND MERGE
    for file in ${in_depth}; do
        filename=\$(basename "\$file" .depth.txt)
        awk -v OFS="\\t" -v filename="\$filename" '{print filename, \$0}' "\$file" >> "merged_heatmap_depth.csv"
    done
    """
}

process GENERATE_STAT_SUMMARY {

    input:
    path bam_file from moved_bam_ch.flatMap().filter { it.name.endsWith('.bam') }

    output:
    path "*.csv" into stat_summary_ch

    script:
    """
    sample_name=\$(basename "${bam_file}")
    weeSAM --bam "$bam_file" --out "\$sample_name".csv
    """

}

process ADD_SAMPLE_COLUMN_TO_SUMMARY {

    input:
    path csv_file from stat_summary_ch

    output:
    path "*.csv" into stat_summary_sample_ch

    script:
    """
    output_file="${csv_file}_fixed.csv"
    python3 ${params.script_files}/add_sample_column.py ${csv_file} \$output_file
    """
 
}

process MERGE_SUMMARY_FILE {

    input:
    path csv_file from stat_summary_sample_ch.collect()

    output:
    path "${params.runname}_merged_summary.csv" into stat_summary_merged_ch

    script:
    """
    awk '!/^Ref_Name/' ${csv_file} > ${params.runname}_merged_summary_temp.csv
    echo -e "Ref_Name\tRef_Len\tMapped_Reads\tBreadth\t%_Covered\tMin_Depth\tMax_Depth\tAvg_Depth\tStd_Dev\tAbove_0.2_Depth\tAbove_1_Depth\tAbove_1.8_Depth\tVariation_Coefficient\tSample" > "header.csv"
    cat "header.csv" ${params.runname}_merged_summary_temp.csv > ${params.runname}_merged_summary.csv
    """
}

process UPDATE_SUMMARY_FILE {

    input:
    path bam_folder from params.out_bam
    path csv_file from moved_stat_ch2.filter{ it.getName() == 'wf-flu-results.csv' }
    path stat_file from moved_stat_ch3.filter { it.name.endsWith('.depth.txt') }
    path tab_file from stat_summary_merged_ch

    output:
    path "*.csv" into stat_summary_updated_ch, stat_summary_updated_ch2

    script:
    """
    stat_folder=${params.out_stat}

    output_file_long="${params.runname}_long_summary.csv" 
    output_file_short="${params.runname}_short_summary.csv" 
    python3 ${params.script_files}/add_subtype_to_summary.py "${csv_file}" "${tab_file}" "\$output_file_long" "\$output_file_short"

    for bam_file in \$(ls ${bam_folder}/*.bam); do
        sample_name=\$(basename "\$bam_file" .bam) 
        stats_file="\$stat_folder/\${sample_name}.stats"
        output_file_sample="\$stat_folder/\${sample_name}_processed.stats"
        output_file_summary="\$stat_folder/\${sample_name}_read_quality.csv"
        python3 "${params.script_files}/read_quality_summaries.py" "\$bam_file" "\$stats_file" "\$output_file_sample" "\$output_file_summary"
    done
    """


}

process MERGE_STAT_SUMMARY {

    input:
    path csv_file from stat_summary_updated_ch

    output:
    path "*_read_quality_summary.csv" into merged_csv_summary_ch
    path "*_processed_summary.stats" into merged_stats_summary_ch

    script:
    """
    stat_folder=${params.out_stat}
    python3 "${params.script_files}/table_merger.py" "processed.stats" "\$stat_folder/${params.runname}_processed_summary.stats" "\$stat_folder"
    python3 "${params.script_files}/table_merger.py" "read_quality.csv" "\$stat_folder/${params.runname}_read_quality_summary.csv" "\$stat_folder"
    
    ln -s "\$stat_folder/${params.runname}_processed_summary.stats" .
    ln -s "\$stat_folder/${params.runname}_read_quality_summary.csv" .
    """
}

process ADD_FRAGMENT_QUALITY_TO_SUMMARY {

    input:
    path csv_file from stat_summary_updated_ch2.flatMap().filter { file -> file.name.endsWith('long_summary.csv') }
    path quality_file from merged_csv_summary_ch

    output:
    path "*csv" into quality_long_summary_ch
    

    script:
    """
    python3 "${params.script_files}/column_lookup_append.py"  \
        "${csv_file}"\
        "${quality_file}" \
        "${params.runname}_long_quality_summary_temp.csv" \
        "," \
        "\t" \
        "sample,Ref_Name" \
        "sample_name,reference" \
        "mean_quality" \
    """
}

process SUBTYPE_B_FIX {

    input:
    path quality_file from quality_long_summary_ch

    output:
    path "*_long_quality_summary.csv" into technical_summary_ch, technical_summary_ch2
    

    script:
    """
    
    python3 "${params.script_files}/B_subtype_fix.py" "${quality_file}"  "${params.runname}_long_quality_summary.csv" 

    """
}

process CLEANUP_OUT_STAT_DIR {

    input:
    val _ from technical_summary_ch
    val out_stat_dir from params.out_stat

    script:
    """
    for file in ${params.out_stat}/*; do
        case \$(basename "\${file}") in
            "new_influenza_pipline_short_summary.csv" | "new_influenza_pipline_read_quality_summary.csv" | "new_influenza_pipline_processed_summary.stats")
                # Do nothing, keep the file
                ;;
            *)
                # Remove the file
                rm "\${file}"
                ;;
        esac
    done
    """
}

process TRANSLATE_FAST_FILES_TO_AMINOACID {
    publishDir params.out_mutation, mode: 'copy'

    errorStrategy 'ignore'

    input:
    path fasta_file from merged_and_extracted_ch.flatten()

    output:
    path "*translation.fasta" into translated_fasta_ch, translated_fasta_ch2

    script:
    """
    fasta_folder="${params.out_fasta}"

    segment_name=\$(basename "${fasta_file}" .fasta)
    seqkit grep -r -p "\$segment_name" "${params.reference}/epi2me/reference_epi2me_FULL_NAMES_mutations.fasta" > "\${segment_name}_temp.fasta"

    nextalign run \
    --input-ref="\${segment_name}_temp.fasta" \
    --genemap=${params.reference}/epi2me/singel_files/"\${segment_name}_genemap.gff" \
    --output-all="\${segment_name}/" \
    "${fasta_file}"

    mv "\${segment_name}"/*translation.fasta "./"
    rm -rf "\${segment_name}"
    rm "\${segment_name}_temp.fasta"

    for translation_file in ./nextalign_gene_*.translation.fasta; do
        new_translation_filename=\$(echo "\$(basename "\$translation_file")" | sed 's/nextalign_gene_//')
        mv "\$translation_file" "./\$new_translation_filename"
    done

    echo "Generated file: ./\${segment_name}.translation.fasta"
    """
}

process GENERATE_MUTATION_LIST {

    publishDir params.out_stat, mode: 'copy'
    
    input:
    path fasta_file from translated_fasta_ch.collect().flatten()


    output:
    path "*.csv" into mutation_singel_summary_ch

    script:
    """
    fasta_name=\$(basename ${fasta_file} .translation.fasta)
    segment=\$(basename "${fasta_file}" .translation.fasta)
    reference="${params.reference}/epi2me/singel_files/"\$segment"_amino.fasta"
    output_file="\$segment.csv"

    python3 "${params.script_files}/mutation_finder.py" ${fasta_file} "\$reference" "\$segment" "\$output_file"
    """
}

process GENERATE_MUTATION_LIST_VACCINE {

    publishDir params.out_stat, mode: 'copy'
    
    input:
    path fasta_file from translated_fasta_ch2.collect().flatten()


    output:
    path "*.csv" into mutation_vaccine_singel_summary_ch

    errorStrategy 'ignore'

    script:
    """
    fasta_name=\$(basename ${fasta_file} .translation.fasta)
    segment=\$(basename "${fasta_file}" .translation.fasta)
    reference="${params.reference}/epi2me/vaccine/"\$segment"_amino.fasta"
    output_file=""\$segment"_vaccine.csv"

    python3 "${params.script_files}/mutation_finder.py" ${fasta_file} "\$reference" "\$segment" "\$output_file"
    """
}

process MERGE_MUTATION_LIST {

    publishDir params.out_mutation, mode: 'copy'
    
    input:
    path csv_file from mutation_singel_summary_ch.collect()
   

    output:
    path "*.csv" into mutation_merged_summary_ch, mutation_merged_summary_ch2, mutation_merged_summary_ch3, mutation_merged_summary_ch4, mutation_merged_summary_ch5, mutation_merged_summary_ch6, mutation_merged_summary_ch7

    script:
    """
    python3 "${params.script_files}/table_merger.py" "csv" "${params.runname}_merged_mutation.csv" "./"
    """
}

process MERGE_MUTATION_LIST_VACCINE {

    publishDir params.out_mutation, mode: 'copy'
    
    input:
    path csv_file from mutation_vaccine_singel_summary_ch.collect()
   

    output:
    path "*.csv" into mutation_merged_summary_vaccine_ch 

    script:
    """
    python3 "${params.script_files}/table_merger.py" "csv" "${params.runname}_merged_mutation_vaccine.csv" "./"
    """
}

process APPEND_MUTATION_LIST_TO_SUMMARY{
    
    input:
    path csv_file from mutation_merged_summary_ch
    path main_summary from technical_summary_ch2

    output:
    path "*.csv" into mutation_add_main_summary_ch

    script:
    """
    python3 "${params.script_files}/column_lookup_append.py"  \
        "${main_summary}" \
        "${csv_file}"  \
        "${params.runname}_long_quality_mutation_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "Differences"
    """
}

process ADD_HA2_MUTATION_LIST_TO_SUMMARY {

    input:
    path summary_file from mutation_add_main_summary_ch
    path csv_file from mutation_merged_summary_ch2

    output:
    path "${params.runname}_long_quality_mutation_summary.csv" into mutation_HA2_add_main_summary_ch

    script:
    """
    python3 "${params.script_files}/HA2_to_summary.py"  \
            "${summary_file}" \
            "${csv_file}" \
            "${params.runname}_long_quality_mutation_summary.csv"
    """
}

process FIND_FLUSERVER_MUTATIONS {

    input:
    path csv_file from mutation_merged_summary_ch3

    output:
    path "${params.runname}_fluserver_mutation.csv" into fluserver_mutation_ch

    script:
    """
    python3 "${params.script_files}/mutation_annotation.py" "${csv_file}" \
        "${params.in_dataset}/RESITENCE_MUTATION/NA_FLUSERVER.csv" \
        "${params.runname}_fluserver_mutation.csv"  \
        "FluserverV1"  \
    """
}

process FIND_PA_MUTATIONS {

    publishDir params.out_mutation, mode: 'copy'

    input:
    path csv_file from mutation_merged_summary_ch7

    output:
    path "${params.runname}_pa_mutation.csv" into pa_mutation_ch

    script:
    """
    python3 "${params.script_files}/mutation_annotation.py" "${csv_file}" \
        "${params.in_dataset}/RESITENCE_MUTATION/PA_RES.csv" \
        "${params.runname}_pa_mutation.csv"  \
        "PA_mutations"  \
    """
}

process FIND_NA_LOW_MUTATIONS {

    input:
    path csv_file from mutation_merged_summary_ch4

    output:
    path "${params.runname}_na_low_mutation.csv" into na_low_mutation_ch

    script:
    """
    python3 "${params.script_files}/mutation_annotation.py" "${csv_file}" \
        "${params.in_dataset}/RESITENCE_MUTATION/NA_LOW.csv" \
        "${params.runname}_na_low_mutation.csv"  \
        "NA_LOW"  \
    """
}

process FIND_NA_MEDIUM_MUTATIONS {

    input:
    path csv_file from mutation_merged_summary_ch5

    output:
    path "${params.runname}_na_medium_mutation.csv" into na_medium_mutation_ch

    script:
    """
    python3 "${params.script_files}/mutation_annotation.py" "${csv_file}" \
        "${params.in_dataset}/RESITENCE_MUTATION/NA_MEDIUM.csv" \
        "${params.runname}_na_medium_mutation.csv"  \
        "NA_MEDIUM"  \
    """
}

process FIND_NA_HIGH_MUTATIONS {

    input:
    path csv_file from mutation_merged_summary_ch6

    output:
    path "${params.runname}_na_high_mutation.csv" into na_high_mutation_ch

    script:
    """
    python3 "${params.script_files}/mutation_annotation.py" "${csv_file}" \
        "${params.in_dataset}/RESITENCE_MUTATION/NA_HIGH.csv" \
        "${params.runname}_na_high_mutation.csv"  \
        "NA_HIGH"  \
    """
}

process APPEND_FLUSERVER_LIST_TO_SUMMARY {

    input:
    path csv_file from fluserver_mutation_ch
    path main_summary from mutation_HA2_add_main_summary_ch

    output:
    path "${params.runname}_long_quality_mutation_summary.csv" into mutation_fluserver_add_main_summary_ch

    script:
    """
    python3 "${params.script_files}/column_lookup_append.py"  \
        "${main_summary}" \
        "${csv_file}"  \
        "${params.runname}_long_quality_mutation_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "FluserverV1" \
    """
}

process APPEND_NA_LOW_LIST_TO_SUMMARY {

    
    input:
    path csv_file from na_low_mutation_ch
    path main_summary from mutation_fluserver_add_main_summary_ch

    output:
    path "${params.runname}_long_quality_mutation_summary.csv" into mutation_na_low_add_main_summary_ch

    script:
    """
    python3 "${params.script_files}/column_lookup_append.py"  \
        "${main_summary}" \
        "${csv_file}"  \
        "${params.runname}_long_quality_mutation_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "NA_LOW" \
    """
}

process APPEND_NA_MEDIUM_LIST_TO_SUMMARY {

    
    input:
    path csv_file from na_medium_mutation_ch
    path main_summary from mutation_na_low_add_main_summary_ch

    output:
    path "${params.runname}_long_quality_mutation_summary.csv" into mutation_na_medium_add_main_summary_ch

    script:
    """
    python3 "${params.script_files}/column_lookup_append.py"  \
        "${main_summary}" \
        "${csv_file}"  \
        "${params.runname}_long_quality_mutation_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "NA_MEDIUM" \
    """
}

process APPEND_NA_HIGH_LIST_TO_SUMMARY {

    
    input:
    path csv_file from na_high_mutation_ch
    path main_summary from mutation_na_medium_add_main_summary_ch

    output:
    path "${params.runname}_long_quality_mutation_summary.csv" into mutation_na_high_add_main_summary_ch

    script:
    """
    python3 "${params.script_files}/column_lookup_append.py"  \
        "${main_summary}" \
        "${csv_file}"  \
        "${params.runname}_long_quality_mutation_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "NA_HIGH" \
    """
}

process FIND_CLADE {

    publishDir params.out_stat, mode: 'copy'
    
    errorStrategy 'ignore'
    
    input:
    path fasta_file from merged_and_extracted_ch.flatten()

    output:
    path "*_CLADE.csv" into clade_ch

    script:
    """
    echo "${fasta_file.name}"

    fasta_name="${fasta_file.name}"
    if [[ "\${fasta_name}" == *H1_HA.fasta ]]; then
        clade=A_H1_HA
        clade_2=A_H1_HA
        nextclade=flu_h1n1pdm_ha
    elif [[ "\${fasta_name}" == *VIC_HA.fasta ]]; then
        clade=B_VIC
        clade_2=B_VIC_HA
        nextclade=flu_vic_ha 
    elif [[ "\${fasta_name}" == A_H3_HA.fasta ]]; then
        clade=A_H3_HA
        clade_2=A_HA_H3
        nextclade=flu_h3n2_ha 
    else
        clade=""
        clade_2=""
    fi

    if [[ "\${clade}" != "" ]]; then
        mkdir -p nextclade_output
        nextclade dataset get --name "\${nextclade}" --output-dir \${nextclade}_nextclade_dataset
        nextclade run \
            --input-dataset \${nextclade}_nextclade_dataset \
            --output-all nextclade_output \
            ${fasta_file}



        python3 "${params.script_files}/clade_table_fix.py" nextclade_output/nextclade.csv \${clade_2}_CLADE.csv 

    fi

    echo "--input-dataset ${params.in_dataset}/HA_NEXTCLADE/\${clade}"

    """
}


process APPEND_CLADES_TO_SUMMARY {

    input:
    path clade_csv from clade_ch.collect()
    path "${params.runname}_long_quality_mutation_summary.csv" from mutation_na_high_add_main_summary_ch

    output:
    path "${params.runname}_clade_summary.csv" into merged_clade_ch

    script:
    """
    mkdir clade_files

    cp -r ${clade_csv} clade_files/

    python3 "${params.script_files}/table_merger.py" "CLADE.csv"  "${params.runname}_merged_clade.csv" "./"

    ls ./

    python3 "${params.script_files}/column_lookup_append.py"  \
        "${params.runname}_long_quality_mutation_summary.csv" \
        "${params.runname}_merged_clade.csv" \
        "${params.runname}_clade_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "clade"

    """
}

process GENERATE_DEPTH_FILES {

    errorStrategy 'ignore'

    input:
    path bam_file from moved_bam_ch2

    output:
    path "*_processed_raw.csv" into depth_files_ch

    script:
    """
    bam_name=\$(basename "${bam_file}" .bam)

    bam-readcount -w1 ${bam_file} > "\${bam_name}_position_count_raw.tsv"

    python3 "${params.script_files}/mutation_ratio_finder.py" "\${bam_name}_position_count_raw.tsv" "\${bam_name}_processed_raw.csv" "\${bam_name}" 823

    """

}

process MERGE_DEPTH_FILES {
    publishDir params.out_stat, mode: 'copy'

    input:
    path csv_file from depth_files_ch.collect()

    output:
    path "${params.runname}_final_merged_depth_raw.csv" into merged_depth_file_ch, merged_depth_file_ch2

    script:
    """
    echo "Creating folder for CSV files..."
    mkdir csv_files

    echo "Moving CSV files to folder..."
    cp -r *csv csv_files/

    echo "Files to be merged:"
    ls -l csv_files/

    echo "Merging process starting..."
    python3 "${params.script_files}/table_merger.py" "_processed_raw.csv" "${params.runname}_merged_depth_raw.csv" "./csv_files/"

    echo "Moving output file to working directory..."
    mv "./csv_files/${params.runname}_merged_depth_raw.csv" "./${params.runname}_final_merged_depth_raw.csv"
    """

}

process APPEND_CYTOSIN_RATIO_823_LIST_TO_SUMMARY {

    input:
    path csv_file from merged_depth_file_ch
    path main_summary from merged_clade_ch

    output:
    path "${params.runname}_cytosine_summary.csv" into summary_cytosin_ch
   

    script:
    """
    python3 "${params.script_files}/column_lookup_append.py"  \
        "${main_summary}" \
        "${csv_file}"  \
        "${params.runname}_cytosine_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "cytosine_ratio__823"
    """
}

process APPEND_THYMINE_RATIO_823_LIST_TO_SUMMARY {
    
    input:
    path csv_file from merged_depth_file_ch2
    path main_summary from summary_cytosin_ch

    output:
    path "${params.runname}_thymine_summary.csv" into summary_thymine_ch

    script:
    """
    python3 "${params.script_files}/column_lookup_append.py"  \
        "${main_summary}" \
        "${csv_file}"  \
        "${params.runname}_thymine_summary.csv" \
        "," \
        "," \
        "sample,Ref_Name" \
        "sample,Ref_Name" \
        "thymine_ratio__823" \
    """
}

process SUMMARY_FILTERING {

    publishDir params.out_stat, mode: 'copy'
        
    input:
    path csv_file from summary_thymine_ch
   
    output:
    path "${params.runname}_summary.csv" into summary_tmep_ch

    script:
    """
    python3 "${params.script_files}/summary_filtering.py"  \
        "${csv_file}" \
        "${params.runname}_summary.csv"  \
    """
}

process FINALIZING_SUMMARY {

    errorStrategy 'ignore'
    
    publishDir params.out_stat, mode: 'copy'
    
    input:
    path csv_file from summary_tmep_ch
    path mutation_file from mutation_merged_summary_vaccine_ch
    path pa_mutations from pa_mutation_ch
   
    output:
    path "${params.runname2}_summary.csv" into summary_ch

    script:
    """
    python3 "${params.script_files}/group_sample_summary.py"  \
        "${csv_file}" \
        "${params.runname2}_summary.csv"  \
        "${params.runname2}"  \
        "${params.script_version}"  \
        "${mutation_file}" \
        "${pa_mutations}" \
    """
}

process FASTA_FILE {

    errorStrategy 'ignore'
    
    publishDir params.out_fasta, mode: 'copy'
    
    input:
    path fasta_file from merged_fasta_ch
   
    output:
    path "${params.runname2}.fasta" into final_fasta_ch

    script:
    """
    python3 "${params.script_files}/final_fasta.py"  \
        "${fasta_file}" \
        "${params.runname2}.fa"
    
    sed '/^>/ s/MP/M/g' "${params.runname2}.fa" > "${params.runname2}.fasta"

    
    """
}

process FIND_COVERAGE {

    
    publishDir params.out_mutation, mode: 'copy'
    
    input:
    path fasta_file from name_added_fasta_files_ch2

    output: 
    path "*.csv" into coverage_ch

    script:
    """
    fasta_file_name=\$(basename ${fasta_file} .fasta)

    python3 "${params.script_files}/coverage_finder.py"  \
            ${fasta_file} \
            "\${fasta_file_name}_coverage.csv"  \

    
    cat "\${fasta_file_name}_coverage.csv"

    """
    
}

process MERGE_COVERAGE {

    publishDir params.out_mutation, mode: 'copy'
    
    input:
    path fasta_file from ncoverage_ch.collect()

    output: 
    path "merged.csv" into merged_coverage_ch

    script:
    """
    # Define the header
    echo "Sample,HA,NA,PA,coverage status,average coverage" > merged.csv

    # Concatenate the content of all files
    for file in ${fasta_file}; do
        cat $file >> merged.csv
    done
    """
}







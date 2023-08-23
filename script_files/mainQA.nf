#!/usr/bin/env nextflow

log.info """\
    I N  H O U S E - N F   P I P E L I N E
    ===================================
    fastq: ${params.in_fastq}
    script: ${params.script_files}
    fastq_out: ${params.out_fastq}
    """.stripIndent()

// Create a channel for each subdirectory that contains *fastq files
fastq_files_ch = Channel.fromPath("${params.in_fastq}/*.fastq")

process REMOVE_DOUBLE_ALIGNES_READS {
    publishDir params.out_fastq, mode: 'copy'

    input:
    path fastq from fastq_files_ch
    path reference from params.reference

    output:
    path "*.fastq"

    script:
    """
    output_file=${fastq.baseName}

    minimap2 -a $reference $fastq | samtools view -bS - | samtools sort -o ${fastq.baseName}.bam -
    samtools index ${fastq.baseName}.bam

    samtools view ${fastq.baseName}.bam | awk '{print \$1}' > ${fastq.baseName}_raw.txt

    python3 "${params.script_files}/read_filtering_alignment_precentage.py" ${fastq.baseName}.bam ${fastq.baseName}_portion_filtered.bam

    # Convert BAM to SAM
    samtools view -h ${fastq.baseName}_portion_filtered.bam > ${fastq.baseName}_portion_filtered.sam

    # Extract read names from \$output_file.sam and sort them
    awk '{print \$1}' ${fastq.baseName}_portion_filtered.sam | sort > ${fastq.baseName}_portion_filtered.txt

    awk '{print \$1 "\t" \$3}' ${fastq.baseName}_portion_filtered.sam > ${fastq.baseName}_reads_refs.txt
    awk '{print \$1}' ${fastq.baseName}_reads_refs.txt | sort | uniq -d > ${fastq.baseName}_multi_align_reads.txt

    sort ${fastq.baseName}_raw.txt > ${fastq.baseName}_raw_sorted.txt
    sort ${fastq.baseName}_portion_filtered.txt > ${fastq.baseName}_portion_filtered_sorted.txt

    # Use comm to get the reads that are unique to fastq (i.e., not in sam)
    comm -23 ${fastq.baseName}_raw_sorted.txt ${fastq.baseName}_portion_filtered_sorted.txt >> ${fastq.baseName}_multi_align_reads.txt

    filterbyname.sh in=$fastq out=${fastq.baseName}_cleaned.fastq names=${fastq.baseName}_multi_align_reads.txt exclude=t
    """
}
















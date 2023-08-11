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
    minimap2 -a $reference $fastq > ${fastq.baseName}.sam
    awk '{print \$1 "\t" \$3}' ${fastq.baseName}.sam > ${fastq.baseName}_reads_refs.txt
    awk '{print \$1}' ${fastq.baseName}_reads_refs.txt | sort | uniq -d > ${fastq.baseName}_multi_align_reads.txt
    filterbyname.sh in=$fastq out=${fastq.baseName}_cleaned.fastq names=${fastq.baseName}_multi_align_reads.txt exclude=t
    """
}














#!/usr/bin/env nextflow

params.outdir = "results"
params.fastq = "$baseDir/data/reads/*.fastq.gz"
params.ref = "$baseDir/data/genome.fna.gz"
params.bowtie2_index = "$baseDir/data/bowtie2_index"
params.script = "$baseDir/bin/deseq2_test.R"
params.data = "$baseDir/data"
params.isomir_sea = "$baseDir/bin/isomiR-SEA_1.6_webpacket/**/*"

process FASTQC {
    executor 'local'

    memory '2 GB'
    cpus 1
    tag "FASTQC on $read"
    publishDir "$params.outdir/fastqc", mode: 'move'

    input:
    path read

    output:
    path "fastqc_${read.baseName}_logs", emit: logs

    script:
    """
    mkdir fastqc_${read.baseName}_logs
    fastqc -o fastqc_${read.baseName}_logs -f fastq -q ${read}
    """
}

process TRIMMING {
    executor 'pbspro'

    memory '1 GB'
    cpus 5
    
    tag "Trimming $read"
    publishDir "$params.outdir/trimming", mode: 'copy'

    input:
    path read

    output:
    path "trimmed_${read.baseName}.gz", emit: trimmed

    script:
    """
    cutadapt --compression-level=2 \
    -a AACTGTAGGCACCATCAAT \
    -g GTTCAGAGTTCTACAGTCCGACGATC \
    -q 20 \
    -m 15 \
    -M 28 \
    -j 5 \
    -o trimmed_${read.baseName}.gz $read
    """
}

process FASTQC2 {
    executor 'local'

    memory '2 GB'
    cpus 1
    tag "FASTQC on $read"
    publishDir "$params.outdir/fastqc2", mode: 'move'

    input:
    path read

    output:
    path "fastqc_${read.baseName}_logs", emit: logs

    script:
    """
    mkdir fastqc_${read.baseName}_logs
    fastqc -o fastqc_${read.baseName}_logs -f fastq -q ${read}
    """
}

process MAPPING {
    executor 'pbspro'

    memory '4 GB'
    cpus 7
    tag "Mapping $read"
    publishDir "$params.outdir/mapped", mode: 'copy'

    input:
    path read

    output:
    path "mapped_${read.baseName}.sam", emit: mapped_sam
    path "mapped_${read.baseName}.fasta.gz", emit: mapped_fa

    script:
    """
    bowtie2 -p 5 \
    --very-sensitive \
    -x $params.bowtie2_index/hsa \
    -U ${read} \
    -S mapped_${read.baseName}.sam \
    --al-gz mapped_${read.baseName}.fasta.gz \
    """
}

process QUANTIFICATION {
    executor 'pbspro'

    memory '1 GB' 
    
    tag "Quantifying $sample_id"
    conda 'bioconda::htseq=2.0.5'
    publishDir "$params.outdir/quantification", mode: 'copy'

    input:
    path gtf
    tuple val(sample_id), path(bam)

    output:
    path "counts_${sample_id}.txt", emit: counts

    script:
    """
    samtools sort -n -o mapped_${sample_id}_sorted.bam ${bam} 
    htseq-count -f bam -r name -m union -s reverse -t gene -i gene_id mapped_${sample_id}_sorted.bam ${gtf} > counts_${sample_id}.txt
    
    rm ${bam} \$(readlink -f ${bam})
    rm mapped_${sample_id}_sorted.bam
    """
}

process DIFFERENTIAL_EXPRESSION {
    executor 'local'

    tag "Differential Expression"
    conda 'bioconda::deseq2=1.42.0'
    publishDir "$params.outdir/deseq2", mode: 'copy'

    input:
    path htseq_outputs 

    output:
    path "deseq2.pdf", emit: deseq2_results

    script:
    """
    mkdir htseq_outputs
    mv counts_*.txt ./htseq_outputs
    Rscript ${params.script} './htseq_outputs' 
    """
}

workflow {
    read_ch = Channel.fromPath(params.fastq, checkIfExists: true)
    FASTQC(read_ch)
    TRIMMING(read_ch)
    FASTQC2(read_ch)
    MAPPING(TRIMMING.out)    

    // For differential expression only 
    // QUANTIFICATION.out
    // .collect()
    // .set{ htseq_outputs }

    // DIFFERENTIAL_EXPRESSION(htseq_outputs)
}
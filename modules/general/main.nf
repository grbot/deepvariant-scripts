#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Get params
ref = file(params.ref, type: 'file')
nr_cores= params.nr_cores

process align {
    tag { "${sample_id}.align" }
    publishDir "${params.out_dir}", mode: 'copy', overwrite: false
    label 'bwa'

    input:
        tuple val(sample_id), path(fastq_r1_file), path(fastq_r2_file)

    output:
        tuple val("$sample_id"), path("${sample_id}.cram"), path("${sample_id}.cram.crai"), emit: cram 

    script:
        flowcell = new Random().nextInt(10000)
        lane = new Random().nextInt(10000) 
        readgroup_info="@RG\\tID:${flowcell}.${lane}\\tLB:LIBA\\tSM:$sample_id\\tPL:Illumina"

        """
        bwa mem \
        -R \"${readgroup_info}\" \
        -t ${nr_cores}  \
        -K 100000000 \
        -Y \
        ${ref} \
        ${fastq_r1_file} \
        ${fastq_r2_file} | \
        samtools sort \
        -@ ${nr_cores} \
        -m ${params.memory_per_thread} \
        -O cram \
        - > ${sample_id}.cram	
	samtools index ${sample_id}.cram
        """ 
}

process call_variants {
    tag { "${sample_id}.call_variants" }
    publishDir "${params.out_dir}/", mode: 'copy', overwrite: false
    label 'deepvariant'

    input:
        tuple val(sample_id), path(cram_file), path(cram_file_index)

    output:
        tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.g.vcf.gz"), emit: vcf

    script:
        """
       /opt/deepvariant/bin/run_deepvariant \
       --model_type=WGS \
       --ref=${ref} \
       --reads=${cram_file} \
       --output_vcf=${sample_id}.vcf.gz \
       --output_gvcf=${sample_id}.g.vcf.gz \
       --num_shards=${nr_cores} \
        """
}


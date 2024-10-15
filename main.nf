#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// For general use
include { align; call_variants } from './modules/general/main.nf'

// Get params

Channel.fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row -> [ "${row.SampleID}",
                       "${row.FastqR1}",
                       "${row.FastqR2}"] }
       .set { samples }


workflow {
  
  bam = align(samples)
  vcf = call_variants(bam)
       
}

workflow.onComplete {
    println ( workflow.success ? """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
     Duration    : ${workflow.duration}
     Success     : ${workflow.success}
     workDir     : ${workflow.workDir}
     exit status : ${workflow.exitStatus}
     """ : """
     Failed: ${workflow.errorReport}
     exit status : ${workflow.exitStatus}
     """
     )
}


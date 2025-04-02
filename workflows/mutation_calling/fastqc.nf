#!/usr/bin/env nextflow

include { FASTQC } from "../../modules/variant_calling/fastqc/main.nf"

// Default parameters
params.samplesheet = "${launchDir}/data/csv/samplesheet.csv"
params.outdir = "${launchDir}/results/fastqc"

workflow {
    // Process samplesheet to create channel with [meta, read] format for each file
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .view {"CSV row: ${it}"}
        .flatMap { row -> 
            def meta = [
                id: row.sample,
                patient: row.patient,
                sample: row.sample
            ]
            
            // Return each file separately with its meta
            return [
                [ meta, file(row.fastq_1) ],
                [ meta, file(row.fastq_2) ]
            ]
        }
        .view {"Channel item: meta=${it[0]}, file=${it[1].name}"}
        .set { reads_ch }
    
    // Run FastQC
    FASTQC(reads_ch)
}
#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Default parameters
params.samplesheet = "/home/zhonggr/projects/250224_DFSP_WES/wes/csv/samplesheet.csv"
params.outdir = "${launchDir}/test/cnv_facets2"
params.dbsnp = params.genomes[params.genome]?.dbsnp
params.dbsnp_tbi = params.genomes[params.genome]?.dbsnp_tbi

// Process: Run CNV-FACETS on paired tumor-normal
process CNV_FACETS {
    
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(input_tumour), path(input_index_tumor), path(input_normal), path(input_index_normal)
    path dbsnp
    path dbsnp_tbi
    
    output:
    tuple val(meta), path("*.cnv.png"),     emit: cnv
    tuple val(meta), path("*.cov.pdf"),     emit: cov
    tuple val(meta), path("*.csv.gz"),      emit: csv
    tuple val(meta), path("*.spider.pdf"),  emit: spider
    tuple val(meta), path("*.vcf.gz"),      emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"),  emit: vcf_tbi
    
    script:
    def prefix = meta.tumour_id

    """
    # Run CNV-FACETS for tumor-normal pair
    cnv_facets.R \\
        --snp-normal ${input_normal} \\
        --snp-tumour ${input_tumour} \\
        --snp-vcf ${dbsnp} \\
        --snp-nprocs ${task.cpus} \\
        --out ${prefix}
    """
}

// Define workflow
workflow {
    // Get reference files
    dbsnp = file(params.dbsnp)
    dbsnp_tbi = file(params.dbsnp_tbi)
    
    // Read and parse samplesheet
    // Create a channel for BAM files
    input_samples = Channel.fromPath(params.samplesheet)
        .ifEmpty { exit(1, "Samplesheet not found: ${params.samplesheet}") }
        .splitCsv(header: true)
        .map { row ->
            // Extract sample info
            def patient_id = row.patient ? row.patient.trim() : null
            def sample_id = row.sample ? row.sample.trim() : null
            def status = row.status ? row.status.trim() : null

            // Validate required fields
            if (!patient_id) {
                error("Missing or empty 'patient' field in row: ${row}")
            }
            if (!sample_id) {
                error("Missing or empty 'sample' field in row: ${row}")
            }
            if (!status) {
                error("Missing or empty 'status' field in row: ${row}")
            }
            if (!row.bam) {
                error("Missing 'bam' field in row: ${row}")
            }

            // Process BAM and BAI paths
            def bam = file(row.bam)
            def bai = row.bai ? file(row.bai) : file("${row.bam}.bai")

            // Check if files exist
            if (!bam.exists()) {
                error("BAM file not found: ${bam}")
            }
            if (!bai.exists()) {
                error("BAI file not found: ${bai}")
            }

            // Return a tuple with patient_id, sample_id, status, bam, and bai
            [
                patient_id: patient_id,
                sample_id: sample_id,
                status: status.toInteger(),
                bam: bam,
                bai: bai,
            ]
        }

        // Split samples into tumour and normal
        tumour_samples = input_samples.filter { it.status == 1 }
        normal_samples = input_samples.filter { it.status == 0 }
    
        // Create paired samples channel
        paired_samples = tumour_samples
            .map { tumour -> [ tumour.patient_id, tumour] }
            .combine(
                normal_samples.map { normal -> [normal.patient_id, normal] },
                by: 0
            )
            .map {
                patient_id, tumour, normal -> 
                def meta = [
                    id: "${tumour.sample_id}_vs_${normal.sample_id}",
                    patient_id: patient_id,
                    tumour_id: tumour.sample_id,
                    normal_id: normal.sample_id,
                    is_paired: true
                ]
                [meta, tumour.bam, tumour.bai, normal.bam, normal.bai]
            }
    
    // Run CNV-FACETS on paired tumor-normal samples
    CNV_FACETS(paired_samples, dbsnp, dbsnp_tbi)
}
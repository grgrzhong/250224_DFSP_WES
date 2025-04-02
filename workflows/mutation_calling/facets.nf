#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Default parameters
// params.samplesheet = "/home/zhonggr/projects/250224_DFSP_WES/data/WES/DFSP/samplesheet.csv"
params.samplesheet = "/home/zhonggr/projects/250224_DFSP_WES/test/SARC/samplesheet.csv"
params.outdir = "${launchDir}/results"
params.dbsnp = params.genomes[params.genome]?.dbsnp
params.dbsnp_tbi = params.genomes[params.genome]?.dbsnp_tbi

// Process: Run CNV-FACETS on paired tumor-normal
process CNV_FACETS {
    
    tag "${meta.id}"
    
    publishDir "${params.outdir}/variant_calling/cnv/facets/${meta.patient_id}/${meta.tumor_id}_vs_${meta.normal_id}", mode: "copy"
    
    // container "/home/zhonggr/projects/250224_DFSP_WES/containers/singularity/cnv_facets-0.16.1.sif"
    
    input:
    tuple val(meta), path(input_normal), path(input_index_normal), path(input_tumor), path(input_index_tumor)
    path dbsnp
    path dbsnp_tbi
    
    output:
    tuple val(meta), path("*.cnv.png"), emit: cnv
    tuple val(meta), path("*.cov.pdf"), emit: cov
    tuple val(meta), path("*.cnv.facets.vcf.gz"), emit: vcf
    tuple val(meta), path("*.cnv.facets.csv.gz"), emit: csv
    
    script:
    def prefix = meta.id

    """
    # Run CNV-FACETS for tumor-normal pair
    cnv_facets.R \\
        --snp-normal ${input_normal} \\
        --snp-tumour ${input_tumor} \\
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
    samples = Channel.fromPath(params.samplesheet)
        .ifEmpty { exit 1, "Samplesheet not found: ${params.samplesheet}" }
        .splitCsv(header: true)
        .map { row ->
            // Extract sample info
            def patient_id = row.patient ? row.patient.trim() : null
            def sample_id = row.sample ? row.sample.trim() : null
            def status = row.status ? row.status.trim() : null
            
            // Validate required fields
            if (!patient_id) error "Missing or empty 'patient' field in row: ${row}"
            if (!sample_id) error "Missing or empty 'sample' field in row: ${row}"
            if (!status) error "Missing or empty 'status' field in row: ${row}"
            if (!row.bam) error "Missing 'bam' field in row: ${row}"
            
            // Process BAM and BAI paths
            def bam = file(row.bam)
            def bai = row.bai ? file(row.bai) : file("${row.bam}.bai")
            
            // Check if files exist
            if (!bam.exists()) error "BAM file not found: ${bam}"
            if (!bai.exists()) error "BAI file not found: ${bai}"
            
            [patient_id, sample_id, status, bam, bai]
        }
        // .view()
    
    // Split samples into tumor and normal
    tumor_samples = samples.filter { it[2] == '1' } // Tumor samples (status == 1)
    normal_samples = samples.filter { it[2] == '0' } // Normal samples (status == 0)
    
    // Match tumor samples with normal samples by patient ID
    tumor_normal_pairs = tumor_samples
        .map { tumor -> [tumor[0], tumor] } // Format as [patient_id, tumor_data]
        .combine(
            normal_samples.map { normal -> [normal[0], normal] }, // Format as [patient_id, normal_data]
            by: 0 // Match by patient_id
        )
        .map { patient_id, tumor, normal ->
            def tumor_id = tumor[1]
            def normal_id = normal[1]
            def tumor_bam = tumor[3]
            def tumor_bai = tumor[4]
            def normal_bam = normal[3]
            def normal_bai = normal[4]
            
            // Create meta map required by the process
            def meta = [
                id: "${tumor_id}_vs_${normal_id}",
                patient_id: patient_id,
                tumor_id: tumor_id,
                normal_id: normal_id
            ]
            
            // log.info "Matched tumor sample ${tumor_id} with normal sample ${normal_id} for patient ${patient_id}"
            
            // Return in format expected by CNV_FACETS process
            [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
        }
    
    log.info "Created ${tumor_normal_pairs.count().val} tumor-normal pairs for analysis"
    
    // Run CNV-FACETS on paired samples
    CNV_FACETS(tumor_normal_pairs, dbsnp, dbsnp_tbi)
}
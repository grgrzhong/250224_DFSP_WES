#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Default parameters
params.samplesheet = "/home/zhonggr/projects/250224_DFSP_WES/test/sarc/samplesheet.csv"
params.outdir = "${launchDir}/results"
params.dbsnp = params.genomes[params.genome]?.dbsnp
params.dbsnp_tbi = params.genomes[params.genome]?.dbsnp_tbi

// Process: Run CNV-FACETS on paired tumor-normal
process CNV_FACETS {
    
    tag "${meta.id}"
    
    input:
    tuple val(meta), path(input_normal), path(input_index_normal), path(input_tumour), path(input_index_tumor)
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
            
            // Return a tuple with patient_id, sample_id, status, bam, and bai
            [
                patient_id: patient_id, 
                sample_id: sample_id, 
                status: status, 
                bam: bam, 
                bai:bai
            ]
        }
        // .view()
    
    // Split samples into tumor and normal
    tumor_samples = samples.filter { it.status == '1' } // Tumor samples (status == 1)
    normal_samples = samples.filter { it.status == '0' } // Normal samples (status == 0)
    
    // tumor_samples.view()
    // normal_samples.view()
     // Extract patient IDs to identify samples without matching normal
    // tumor_patient_ids = tumor_samples.map { sample -> sample.patient_id }.unique()
    // normal_patient_ids = normal_samples.map { sample -> sample.patient_id }.unique()
    
    // Find and report tumor samples without matching normal samples
    // tumor_patient_ids
    //     .combine(normal_patient_ids.collect())
    //     .map { tumor_patient_id, all_normal_patient_ids ->
    //         def has_matching_normal = all_normal_patient_ids.contains(tumor_patient_id)
    //         [tumor_patient_id, has_matching_normal]
    //     }
    //     .filter { tumor_patient_id, has_matching_normal -> !has_matching_normal }
    //     .map { tumor_patient_id, _ -> 
    //         log.warn "WARNING: Tumor patient ${tumor_patient_id} has no matching normal sample"
    //         return tumor_patient_id 
    //     }
    //     .count()
    //     .view { unmatched_count -> 
    //         if (unmatched_count > 0) {
    //             return "WARNING: Found ${unmatched_count} tumor patients without matching normal samples"
    //         } else {
    //             return "All tumor patients have matching normal samples"
    //         }
    //     }

    // Match tumor samples with normal samples by patient ID
    tumor_normal_pairs = tumor_samples
        .map { tumor -> [tumor.patient_id, tumor] } // Format as [patient_id, tumor_data]
        .combine(
            normal_samples.map { normal -> [normal.patient_id, normal] }, // Format as [patient_id, normal_data]
            by: 0 // Match by patient_id
        )
        .map { patient_id, tumor, normal ->
            def tumour_id = tumor.sample_id
            def normal_id = normal.sample_id
            def tumor_bam = tumor.bam
            def tumor_bai = tumor.bai
            def normal_bam = normal.bam
            def normal_bai = normal.bai
            
            // Create meta map required by the process
            def meta = [
                id: "${tumour_id}_vs_${normal_id}",
                patient_id: patient_id,
                tumour_id: tumour_id,
                normal_id: normal_id
            ]
            
            // log.info "Matched tumor sample ${tumour_id} with normal sample ${normal_id} for patient ${patient_id}"
            
            // Return in format expected by CNV_FACETS process
            [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
        }
    
    // log.info "Created ${tumor_normal_pairs.count().val} tumor-normal pairs for analysis"
    // tumor_normal_pairs.view()
    // Run CNV-FACETS on paired samples
    CNV_FACETS(tumor_normal_pairs, dbsnp, dbsnp_tbi)
}
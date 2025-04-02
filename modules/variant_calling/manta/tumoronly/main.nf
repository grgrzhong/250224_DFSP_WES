#!/usr/bin/env nextflow

// Parameters
params.samplesheet = "${launchDir}/data/WES/DFSP/samplesheet_DFSP.csv"
params.outdir = "${launchDir}/results/variant_calling/manta"
params.reference = "/path/to/your/reference/genome.fa" // Replace with actual reference path


// Define process to run Manta for tumor-normal pairs
process MANTA_TUMOR_NORMAL {
    tag "${patient}_T${tumor_id}_N${normal_id}"
    
    publishDir "${params.outdir}/manta/${tumor_id}", mode: 'copy'
    
    container "${params.singularity_container_dir}/manta.1.6.0.img"
    
    input:
    tuple val(patient), val(tumor_id), path(tumor_bam), path(tumor_bai), val(normal_id), path(normal_bam), path(normal_bai)
    path reference
    
    output:
    tuple val(patient), val(tumor_id), val(normal_id), path("${patient}_T${tumor_id}_N${normal_id}/results/variants/diploidSV.vcf.gz"), emit: diploid_sv
    tuple val(patient), val(tumor_id), val(normal_id), path("${patient}_T${tumor_id}_N${normal_id}/results/variants/somaticSV.vcf.gz"), emit: somatic_sv
    tuple val(patient), val(tumor_id), val(normal_id), path("${patient}_T${tumor_id}_N${normal_id}/results/variants/candidateSV.vcf.gz"), emit: candidate_sv
    tuple val(patient), val(tumor_id), val(normal_id), path("${patient}_T${tumor_id}_N${normal_id}/results/variants/*"), emit: all_variants
    
    script:
    """
    # Configure Manta for tumor-normal analysis
    configManta.py \
        --normalBam ${normal_bam} \
        --tumorBam ${tumor_bam} \
        --referenceFasta ${reference} \
        --runDir ${patient}_T${tumor_id}_N${normal_id}
    
    # Run Manta using available threads
    python ${patient}_T${tumor_id}_N${normal_id}/runWorkflow.py \
        --mode local \
        --jobs ${task.cpus}
    """
}

// Define process to run Manta for tumor-only samples
process MANTA_TUMOR_ONLY {
    tag "${patient}_T${tumor_id}"
    
    publishDir "${params.outdir}/manta/${patient}", mode: 'copy'
    
    container "${params.singularity_container_dir}/manta.sif"
    
    cpus params.threads
    
    input:
    tuple val(patient), val(tumor_id), path(tumor_bam), path(tumor_bai)
    path reference
    
    output:
    tuple val(patient), val(tumor_id), path("${patient}_T${tumor_id}/results/variants/tumorSV.vcf.gz"), emit: tumor_sv
    tuple val(patient), val(tumor_id), path("${patient}_T${tumor_id}/results/variants/candidateSV.vcf.gz"), emit: candidate_sv
    tuple val(patient), val(tumor_id), path("${patient}_T${tumor_id}/results/variants/*"), emit: all_variants
    
    script:
    """
    # Configure Manta for tumor-only analysis
    configManta.py \
        --tumorBam ${tumor_bam} \
        --referenceFasta ${reference} \
        --runDir ${patient}_T${tumor_id}
    
    # Run Manta using available threads
    python ${patient}_T${tumor_id}/runWorkflow.py \
        --mode local \
        --jobs ${task.cpus}
    """
}

workflow {
    // Read sample information from the CSV file
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row -> 
            def patient = row.patient
            def sample = row.sample
            def status = row.status.toInteger()  // 0 for normal, 1 for tumor
            def bam = file(row.bam)
            def bai = file(row.bai)
            
            // Return [patient, sample, status, bam, bai]
            return [patient, sample, status, bam, bai]
        }
        .branch {
            normals: it[2] == 0  // Normal samples
            tumors: it[2] == 1   // Tumor samples
        }
        .set { samples }

    // Create channel for reference genome
    reference_ch = Channel.value(file(params.reference))
    
    // Group normal samples by patient
    normals_by_patient = samples.normals
        .map { patient, sample, status, bam, bai -> [patient, sample, bam, bai] }
        .groupTuple(by: 0)
    
    // Group tumor samples by patient
    tumors_by_patient = samples.tumors
        .map { patient, sample, status, bam, bai -> [patient, sample, bam, bai] }
        .groupTuple(by: 0)
    
    // Combine tumor and normal samples by patient for paired analysis
    paired_samples = tumors_by_patient
        .join(normals_by_patient, by: 0, remainder: true)
        .map { patient, tumor_info, normal_info ->
            if (normal_info) {
                // For each tumor sample in this patient
                def result = []
                tumor_info[0].eachWithIndex { tumor_sample, t_idx ->
                    // Use the first normal sample for this patient
                    def normal_sample = normal_info[0][0]
                    def normal_bam = normal_info[2][0]
                    def normal_bai = normal_info[3][0]
                    
                    // Create a tumor-normal pair
                    result << [patient, tumor_sample, tumor_info[2][t_idx], tumor_info[3][t_idx], normal_sample, normal_bam, normal_bai]
                }
                return result
            } else {
                // For tumor-only analysis
                def result = []
                tumor_info[0].eachWithIndex { tumor_sample, t_idx ->
                    result << [patient, tumor_sample, tumor_info[2][t_idx], tumor_info[3][t_idx]]
                }
                return result
            }
        }
        .flatMap { it }
    
    // Split into tumor-normal pairs and tumor-only samples
    paired_samples
        .branch {
            tumor_normal: it.size() == 7
            tumor_only: it.size() == 4
        }
        .set { analysis_samples }
    
    // Run Manta on tumor-normal pairs
    MANTA_TUMOR_NORMAL(
        analysis_samples.tumor_normal,
        reference_ch
    )
    
    // Run Manta on tumor-only samples
    MANTA_TUMOR_ONLY(
        analysis_samples.tumor_only,
        reference_ch
    )
}
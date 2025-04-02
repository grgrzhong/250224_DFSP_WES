#!/usr/bin/env nextflow

// Default parameter values
params.samplesheet = "${launchDir}/samplesheet.csv"
params.outdir = "${launchDir}/results"
params.facets_vcf = "/path/to/reference/SNP_vcf.vcf.gz"
params.publishDirMode = "copy"
params.genome = "hg38"

// FACETS parameters
params.cval = 100
params.purity_cval = 200
params.snp_nbhd = 250
params.ndepth = 35
params.min_nhet = 15
params.purity_min_nhet = 15
params.seed = 1234

// Read samplesheet
Channel
    .fromPath(params.samplesheet)
    .ifEmpty { exit 1, "Samplesheet not found: ${params.samplesheet}" }
    .splitCsv(header: true)
    .map { row -> 
        def sample_id = row.sample
        def patient_id = row.patient
        def status = row.status
        def bam_file = file(row.bam)
        def bai_file = file(row.bai)

        return [patient_id, sample_id, status, bam_file, bai_file]
    }
    .branch {
        tumor: it[2] == '1'
        normal: it[2] == '0'
    }
    .set { samples_ch }

// Group tumor-normal pairs by patient
samples_ch.tumor
    .combine(samples_ch.normal, by: 0)  // Join by patient_id
    .map { patient_id, tumor_id, tumor_status, tumor_bam, tumor_bai, normal_id, normal_status, normal_bam, normal_bai ->
        [tumor_id, normal_id, patient_id, tumor_bam, tumor_bai, normal_bam, normal_bai]
    }
    .set { tumor_normal_pairs }

// Process: Run FACETS
process FACETS {
    tag "${tumor_id}_vs_${normal_id}"
    publishDir "${params.outdir}/${patient_id}/${tumor_id}_vs_${normal_id}/facets", 
               mode: params.publishDirMode,
               pattern: "*.{png,txt,seg,out,Rdata,gz}"

    input:
    tuple val(tumor_id), val(normal_id), val(patient_id), 
          file(tumor_bam), file(tumor_bai), file(normal_bam), file(normal_bai)
    path facets_vcf from Channel.value(file(params.facets_vcf))

    output:
    tuple val(tumor_id), val(normal_id), val(patient_id), 
          file("${pair_id}.snp_pileup.gz"),
          file("${output_dir}/*.png"),
          file("${output_dir}/*.seg"),
          file("${output_dir}/*.out"),
          file("${output_dir}/*cncf.txt")
    file("${pair_id}_facets_qc.txt")

    script:
    pair_id = "${tumor_id}_vs_${normal_id}"
    output_dir = "facets_c${params.cval}_pc${params.purity_cval}"
    """
    # Run SNP pileup
    export SNP_PILEUP=/usr/bin/snp-pileup
    Rscript /usr/bin/facets-suite/snp-pileup-wrapper.R \
        --pseudo-snps 50 \
        --vcf-file ${facets_vcf} \
        --output-prefix ${pair_id} \
        --normal-bam ${normal_bam} \
        --tumor-bam ${tumor_bam}

    # Create output directory
    mkdir -p ${output_dir}

    # Run FACETS
    Rscript /usr/bin/facets-suite/run-facets-wrapper.R \
        --cval ${params.cval} \
        --snp-window-size ${params.snp_nbhd} \
        --normal-depth ${params.ndepth} \
        --min-nhet ${params.min_nhet} \
        --purity-cval ${params.purity_cval} \
        --purity-min-nhet ${params.purity_min_nhet} \
        --genome ${params.genome} \
        --counts-file ${pair_id}.snp_pileup.gz \
        --sample-id ${pair_id} \
        --directory ${output_dir} \
        --facets-lib-path /usr/local/lib/R/site-library \
        --seed ${params.seed} \
        --everything \
        --legacy-output T

    # Generate summary
    python3 /usr/bin/summarize_project.py \
        -p ${pair_id} \
        -c ${output_dir}/*cncf.txt \
        -o ${output_dir}/*out \
        -s ${output_dir}/*seg

    # Generate QC report
    mkdir -p refit_watcher/bin/ refit_watcher/refit_jobs/
    R -e "facetsPreview::generate_genomic_annotations('${pair_id}', '\$(pwd)/', '/usr/bin/facets-preview/tempo_config.json')"
    cp facets_qc.txt ${pair_id}_facets_qc.txt
    """
}

workflow {
    FACETS(tumor_normal_pairs, params.facets_vcf)
}
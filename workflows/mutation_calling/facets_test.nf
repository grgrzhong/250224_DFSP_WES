#!/usr/bin/env nextflow

// Default parameters
params.normal_bam = "/home/zhonggr/projects/250224_DFSP_WES/test/SARC/bam/SARC-004-N/SARC-004-N_recalibrated.bam"
params.normal_bai = "/home/zhonggr/projects/250224_DFSP_WES/test/SARC/bam/SARC-004-N/SARC-004-N_recalibrated.bai"
params.tumour_bam = "/home/zhonggr/projects/250224_DFSP_WES/test/SARC/bam/SARC-004-T/SARC-004-T_recalibrated.bam"
params.tumour_bai = "/home/zhonggr/projects/250224_DFSP_WES/test/SARC/bam/SARC-004-T/SARC-004-T_recalibrated.bai"
params.outdir = "${launchDir}/results/variant_calling/cnv/facets"
params.snp_vcf = "/home/zhonggr/projects/250224_DFSP_WES/data/Reference/dbSNP.vcf.gz"
params.snp_vcf_tbi = "/home/zhonggr/projects/250224_DFSP_WES/data/Reference/dbSNP.vcf.gz.tbi"

// Process: Run CNV-FACETS on paired tumor-normal
process CNV_FACETS {
    
    publishDir params.outdir, mode: "copy", overwrite: true
    container "/home/zhonggr/projects/250224_DFSP_WES/containers/singularity/cnv_facets-0.16.1.sif"
    
    input:
    path tumor_bam
    path tumor_bai
    path normal_bam
    path normal_bai
    path snp_vcf
    path snp_vcf_tbi
    
    output:
        path("*.{cnv.png,cov.pdf,vcf.gz,csv.gz}")
    
    script:
    """
    # Run CNV-FACETS
    cnv_facets.R \\
        --snp-tumour ${tumor_bam} \\
        --snp-normal ${normal_bam} \\
        --snp-vcf ${snp_vcf} \\
        --snp-nprocs $task.cpus\\
        --out "${tumor_bam.baseName}_vs_${normal_bam.baseName}"
        
    """
}

// Define workflow
workflow {
    // Get reference files
    snp_vcf = file(params.snp_vcf)
    snp_vcf_tbi = file(params.snp_vcf_tbi)

    tumour_bam_ch = file(params.tumour_bam)
    tumour_bai_ch = file(params.tumour_bai)
    normal_bam_ch = file(params.normal_bam)
    normal_bai_ch = file(params.normal_bai)

    // Run CNV-FACETS on valid tumor-normal pairs
    CNV_FACETS(tumour_bam_ch, tumour_bai_ch, normal_bam_ch, normal_bai_ch, snp_vcf, snp_vcf_tbi)
}
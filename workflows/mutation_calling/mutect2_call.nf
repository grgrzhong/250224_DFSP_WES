#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Load modules
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_TUMOR  } from '../../modules/variant_calling/gatk4/getpileupsummaries/main.nf'
include { GATK4_GETPILEUPSUMMARIES as GETPILEUPSUMMARIES_NORMAL } from '../../modules/variant_calling/gatk4/getpileupsummaries/main.nf'
include { GATK4_CALCULATECONTAMINATION                          } from '../../modules/variant_calling/gatk4/calculatecontamination/main.nf'
include { GATK4_MUTECT2                                         } from '../../modules/variant_calling/gatk4/mutect2/main.nf'
include { GATK4_LEARNREADORIENTATIONMODEL                       } from '../../modules/variant_calling/gatk4/learnreadorientationmodel/main.nf'
include { GATK4_FILTERMUTECTCALLS                               } from '../../modules/variant_calling/gatk4/filtermutectcalls/main.nf'
include { BCFTOOLS_NORM                                         } from '../../modules/variant_calling/bcftools/norm/main.nf'
include { BCFTOOLS_VIEW                                         } from '../../modules/variant_calling/bcftools/view/main.nf'
include { GATK4_FUNCOTATOR                                      } from '../../modules/variant_calling/gatk4/funcotator/main.nf'
include { ANNOVAR                                               } from '../../modules/variant_calling/annovar/main.nf'

workflow SOMATIC_VARIANT_CALLING {
    take:
    input_sample          // Channel: [val(meta), path(bam), path(bai), val(normal_meta), path(normal_bam), path(normal_bai)]
    fasta                 // Channel: path(fasta)
    fai                   // Channel: path(fai)
    dict                  // Channel: path(dict)
    germline_resource     // Channel: path(germline_resource_vcf)
    germline_resource_tbi // Channel: path(germline_resource_tbi)
    panel_of_normals      // Channel: path(panel_of_normals)
    panel_of_normals_tbi  // Channel: path(panel_of_normals_tbi)
    interval_list         // Channel: path(interval_list)
    small_exac            // Channel: path(small_exac_common)
    small_exac_tbi        // Channel: path(small_exac_common_tbi)
    funcotator_sources    // Channel: path(funcotator_sources)
    annovar_db            // Channel: path(annovar_db)

    main:
    ch_versions = Channel.empty()

    // Get pileup for tumor samples
    GETPILEUPSUMMARIES_TUMOR (
        input_sample.map { meta, bam, bai, normal_meta, normal_bam, normal_bai ->
            [meta, bam]
        },
        small_exac,
        small_exac_tbi,
        small_exac // Use as intervals
    )
    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_TUMOR.out.versions)

    // [... rest of the pipeline same as before until FilterMutectCalls ...]

    // Normalize VCFs using bcftools norm
    BCFTOOLS_NORM (
        GATK4_FILTERMUTECTCALLS.out.vcf.map { meta, vcf -> [meta, vcf] },
        fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    // Filter normalized VCFs to keep only PASS variants
    BCFTOOLS_VIEW (
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi, by: 0)
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    // Funcotator annotation - now working with filtered VCFs from BCFTOOLS_VIEW
    GATK4_FUNCOTATOR (
        BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi, by: 0),
        fasta,
        fai,
        dict,
        funcotator_sources,
        interval_list,
        "MAF"
    )
    ch_versions = ch_versions.mix(GATK4_FUNCOTATOR.out.versions)

    // Annovar annotation - also working with filtered VCFs
    ANNOVAR (
        BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi, by: 0),
        annovar_db,
        "hg38",
        "refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70",
        "gx,r,f,f,f,f,f"
    )
    ch_versions = ch_versions.mix(ANNOVAR.out.versions)

    emit:
    filtered_vcf     = GATK4_FILTERMUTECTCALLS.out.vcf       // [meta, vcf]
    normalized_vcf   = BCFTOOLS_NORM.out.vcf                 // [meta, vcf]
    filtered_norm_vcf = BCFTOOLS_VIEW.out.vcf                // [meta, vcf]
    funcotator_maf   = GATK4_FUNCOTATOR.out.annotated_variants // [meta, maf]
    annovar_txt      = ANNOVAR.out.annovar_txt               // [meta, txt]
    versions         = ch_versions                           // path: versions.yml
}
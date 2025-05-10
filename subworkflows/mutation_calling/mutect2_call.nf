
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                GATK4 Mutect2 SNV/Indels somatic variant calling workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load required modules
include { GATK4_MUTECT2_TUMOUR as MUTECT2_TUMOUR                 } from '../../modules/variant_calling/gatk4/mutect2/tumour'
include { GATK4_MUTECT2_PAIRED as MUTECT2_PAIRED                 } from '../../modules/variant_calling/gatk4/mutect2/paired'
include { GATK4_GETPILEUPSUMMARIES as PILEUP_PAIRED_TUMOUR       } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as PILEUP_PAIRED_NORMAL       } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as PILEUPS_UNPAIRED_TUMOUR    } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_CALCULATECONTAMINATION as CONTAMINATION_PAIRED   } from '../../modules/variant_calling/gatk4/calculatecontamination'
include { GATK4_CALCULATECONTAMINATION as CONTAMINATION_UNPAIRED } from '../../modules/variant_calling/gatk4/calculatecontamination'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNMODEL_PAIRED   } from '../../modules/variant_calling/gatk4/learnreadorientationmodel'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNMODEL_UNPAIRED } from '../../modules/variant_calling/gatk4/learnreadorientationmodel'
include { GATK4_FILTERMUTECTCALLS as FILTERMUTECTCALLS           } from '../../modules/variant_calling/gatk4/filtermutectcalls'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_MUTECT2                 } from '../../modules/variant_calling/bcftools/norm'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_MUTECT2                 } from '../../modules/variant_calling/bcftools/view'
include { TABIX as TABIX_MUTECT2                                 } from '../../modules/variant_calling/tabix/tabix'
include { GATK4_FUNCOTATOR as FUNCOTATOR                         } from '../../modules/variant_calling/gatk4/funcotator'

workflow MUTECT2_CALL {
    take:
    fasta
    fai
    dict
    pileup_variants
    pileup_variants_tbi
    germline_resource
    germline_resource_tbi
    panel_of_normals
    panel_of_normals_tbi
    intervals
    bam_tumour_normal
    bam_tumour_only

    main:
    
    // Initialize empty channels for results
    tumour_normal_vcf = Channel.empty()
    tumour_normal_tbi = Channel.empty()
    tumour_normal_stats = Channel.empty()
    tumour_normal_orientation = Channel.empty()
    tumour_normal_contamination = Channel.empty()
    tumour_normal_segmentation = Channel.empty()
    
    tumour_only_vcf = Channel.empty()
    tumour_only_tbi = Channel.empty()
    tumour_only_stats = Channel.empty()
    tumour_only_orientation = Channel.empty()
    tumour_only_contamination = Channel.empty()
    tumour_only_segmentation = Channel.empty()

    /*
     * =================== TUMOR-NORMAL PAIRED ANALYSIS ================
     */
    
    // Extract tumour samples for pileup
    paired_tumour_samples = bam_tumour_normal
        .map { meta, tumour_bam, tumour_bai, _normal_bam, _normal_bai ->
            [meta, tumour_bam, tumour_bai]
        }
    
    // Extract normal samples for pileup
    paired_normal_samples = bam_tumour_normal
        .map { meta, _tumour_bam, _tumour_bai, normal_bam, normal_bai ->
                [meta, normal_bam, normal_bai]
        }
    
    // Get pileup summaries for tumour samples
    paired_tumour_pileup = PILEUP_PAIRED_TUMOUR(
        paired_tumour_samples, 
        pileup_variants,
        pileup_variants_tbi
    )
    
    // Get pileup summaries for normal samples
    paired_normal_pileup = PILEUP_PAIRED_NORMAL(
        paired_normal_samples,
        pileup_variants,
        pileup_variants_tbi
    )

    // Calculate contamination
    paired_tumour_normal_pileup = paired_tumour_pileup
        .table
        .join(paired_normal_pileup.table, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, t_table, n_table -> [meta, t_table, n_table] }
        
    paired_contamination = CONTAMINATION_PAIRED(
        paired_tumour_normal_pileup
    )
    
    // Run Mutect2
    paired_mutect2 = MUTECT2_PAIRED(
        paired_tumour_samples,
        paired_normal_samples,
        fasta,
        fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi,
        intervals
    )
    
    // Learn read orientation model
    paired_orientation = LEARNMODEL_PAIRED(paired_mutect2.f1r2)
    
    // Store output channels
    tumour_normal_vcf = paired_mutect2.vcf
    tumour_normal_tbi = paired_mutect2.tbi
    tumour_normal_stats = paired_mutect2.stats
    tumour_normal_orientation = paired_orientation.orientation
    tumour_normal_contamination = paired_contamination.contamination
    tumour_normal_segmentation = paired_contamination.segmentation

    /*
     * =================== TUMOR-ONLY UNPAIRED ANALYSIS ====================
     */
    
    unpaired_tumour_samples = bam_tumour_only
        .map { meta, tumour_bam, tumour_bai, _normal_bam, _normal_bai ->
            [meta, tumour_bam, tumour_bai]
        }

    // Get pileup summaries for tumor-only samples
    unpaired_tumour_pileup = PILEUPS_UNPAIRED_TUMOUR(
        unpaired_tumour_samples,
        pileup_variants,
        pileup_variants_tbi
    )
    
    // Calculate contamination (no matched normal)
    tumour_only_input = unpaired_tumour_pileup.table
        .map { meta, table -> [meta, table, []] }

    unpaired_contamination = CONTAMINATION_UNPAIRED(tumour_only_input)

    
    // Run Mutect2 for tumor-only samples
    unpaired_mutect2 = MUTECT2_TUMOUR(
        unpaired_tumour_samples,
        fasta,
        fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi,
        intervals
    )
    
    // Learn read orientation model
    unpaired_orientation = LEARNMODEL_UNPAIRED(unpaired_mutect2.f1r2)
    
    // Store output channels
    tumour_only_vcf = unpaired_mutect2.vcf
    tumour_only_tbi = unpaired_mutect2.tbi
    tumour_only_stats = unpaired_mutect2.stats
    tumour_only_orientation = unpaired_orientation.orientation
    tumour_only_contamination = unpaired_contamination.contamination
    tumour_only_segmentation = unpaired_contamination.segmentation
    
    // Combine all Mutect2 outputs for further processing
    mutect2_vcf = tumour_normal_vcf.mix(tumour_only_vcf)
    mutect2_tbi = tumour_normal_tbi.mix(tumour_only_tbi)
    mutect2_stats = tumour_normal_stats.mix(tumour_only_stats)
    mutect2_orientation = tumour_normal_orientation.mix(tumour_only_orientation)
    mutect2_contamination = tumour_normal_contamination.mix(tumour_only_contamination)
    mutect2_segmentation = tumour_normal_segmentation.mix(tumour_only_segmentation)
    
    // Prepare inputs for FilterMutectCalls
    filter_input = mutect2_vcf
        .join(mutect2_tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(mutect2_stats, failOnDuplicate: true, failOnMismatch: true)
        .join(mutect2_orientation, failOnDuplicate: true, failOnMismatch: true)
        .join(mutect2_contamination, failOnDuplicate: true, failOnMismatch: true)
        .join(mutect2_segmentation, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi, stats, orientation, contamination, segmentation ->
            [meta, vcf, tbi, stats, orientation, contamination, segmentation]
        }
    
    // Filter Mutect2 calls
    filtered_calls = FILTERMUTECTCALLS(
        filter_input,
        fasta,
        fai,
        dict
    )
    
    // Normalize variants with bcftools
    normalized_input = filtered_calls
        .vcf
        .join(filtered_calls.tbi, failOnDuplicate: true, failOnMismatch: true)
    
    normalized_vcf = BCFTOOLS_NORM_MUTECT2(
        normalized_input,
        fasta,
        fai,
        dict
    )
    
    // Filter for PASS variants with bcftools
    pass_vcf = BCFTOOLS_VIEW_MUTECT2(normalized_vcf.vcf)
    
    // Index the filtered VCF
    indexed_vcf = TABIX_MUTECT2(pass_vcf.vcf)

    emit:
    vcf = pass_vcf.vcf                // channel: [ meta, vcf ]
    vcf_tbi = indexed_vcf.tbi             // channel: [ meta, tbi ]
}    

include { PREPARE_GENOME } from '../../subworkflows/mutation_calling/prepare_genome'
include { PREPARE_SAMPLE } from '../../subworkflows/mutation_calling/prepare_sample'


workflow {
    
    PREPARE_GENOME(params.genome)

    PREPARE_SAMPLE(params.input)
    
    // Run Mutect2 somatic variant calling
    MUTECT2_CALL(
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.dict,
        PREPARE_GENOME.out.pileup_variants,
        PREPARE_GENOME.out.pileup_variants_tbi,
        PREPARE_GENOME.out.germline_resource,
        PREPARE_GENOME.out.germline_resource_tbi,
        PREPARE_GENOME.out.panel_of_normals,
        PREPARE_GENOME.out.panel_of_normals_tbi,
        PREPARE_GENOME.out.intervals,
        PREPARE_SAMPLE.out.bam_tumour_normal,
        PREPARE_SAMPLE.out.bam_tumour_only
    )
}
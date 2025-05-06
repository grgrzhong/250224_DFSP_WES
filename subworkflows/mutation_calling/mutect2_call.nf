/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                GATK4 Mutect2 SNV/Indels somatic variant calling workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load required modules
include { GATK4_GETPILEUPSUMMARIES as  GATK4_GETPILEUPSUMMARIES_PAIRED_TUMOUR               } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as  GATK4_GETPILEUPSUMMARIES_PAIRED_NORMAL               } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as  GATK4_GETPILEUPSUMMARIES_TUMOUR_ONLY                 } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_CALCULATECONTAMINATION as GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL        } from '../../modules/variant_calling/gatk4/calculatecontamination'
include { GATK4_CALCULATECONTAMINATION as GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY          } from '../../modules/variant_calling/gatk4/calculatecontamination'
include { GATK4_MUTECT2 as GATK4_MUTECT2_TUMOR_NORMAL                                       } from '../../modules/variant_calling/gatk4/mutect2'
include { GATK4_MUTECT2 as GATK4_MUTECT2_TUMOR_ONLY                                         } from '../../modules/variant_calling/gatk4/mutect2'
include { GATK4_LEARNREADORIENTATIONMODEL as GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_NORMAL  } from '../../modules/variant_calling/gatk4/learnreadorientationmodel'
include { GATK4_LEARNREADORIENTATIONMODEL as GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_ONLY    } from '../../modules/variant_calling/gatk4/learnreadorientationmodel'
include { GATK4_FILTERMUTECTCALLS                                                           } from '../../modules/variant_calling/gatk4/filtermutectcalls'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_MUTECT2                                            } from '../../modules/variant_calling/bcftools/norm'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_MUTECT2                                            } from '../../modules/variant_calling/bcftools/view'
include { TABIX as TABIX_MUTECT2                                                            } from "../../modules/variant_calling/tabix/tabix"

workflow MUTECT2_CALL{

    take:
        // Reference genome files
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

        // Sample data
        paired_samples
        unpaired_samples
    
    main:
        /*
            ===================== TUMOR-NORMAL PAIRED ANALYSIS ====================
        */
        if (paired_samples) {

            // Extract tumor samples from paired samples for pileup
            paired_tumour_samples = paired_samples
                .map { 
                    meta, tumour_bam, tumour_bai, _normal_bam, _normal_bai ->
                    [meta, tumour_bam, tumour_bai]
            }
            paired_tumour_samples.view()

            // Get pileup summaries for tumor samples in paired mode
            GATK4_GETPILEUPSUMMARIES_PAIRED_TUMOUR(
                paired_tumour_samples,
                pileup_variants,
                pileup_variants_tbi
            )

            // Extract normal samples from paired samples for pileup
            paired_normal_samples = paired_samples
                .map { meta, _tumour_bam, _tumour_bai, normal_bam, normal_bai ->
                    def normal_meta = meta.clone()
                    normal_meta.tumour_id = normal_meta.normal_id
                    [normal_meta, normal_bam, normal_bai]
            }
            paired_normal_samples.view()

            GATK4_GETPILEUPSUMMARIES_PAIRED_NORMAL(
                paired_normal_samples,
                pileup_variants,
                pileup_variants_tbi
            )

            // Calculate contamination for paired samples
            GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL(
                GATK4_GETPILEUPSUMMARIES_PAIRED_TUMOUR.out.table
                .join(GATK4_GETPILEUPSUMMARIES_PAIRED_NORMAL.out.table, by: [0])
                .map { meta, tumor_table, normal_table ->
                    [meta, tumor_table, normal_table]
                }
            )

            // Run Mutect2 for paired samples
            GATK4_MUTECT2_TUMOR_NORMAL(
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

            GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_NORMAL(
                GATK4_MUTECT2_TUMOR_NORMAL.out.f1r2
            )

        }

        /*
            =================== TUMOR-ONLY UNPAIRED ANALYSIS ====================
        */
        if (unpaired_samples) {
        
            // Extract tumor samples for unpaired analysis
            unpaired_tumour_samples = unpaired_samples
                .map { meta, tumor_bam, tumor_bai, _normal_bam, _normal_bai ->
                [meta, tumor_bam, tumor_bai]
            }
            unpaired_tumour_samples.view()

            // Get pileup summaries for tumor-only samples
            GATK4_GETPILEUPSUMMARIES_TUMOUR_ONLY(
                unpaired_tumour_samples,
                pileup_variants,
                pileup_variants_tbi
            )

            // Calculate contamination for tumor-only samples
            GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY(
                GATK4_GETPILEUPSUMMARIES_TUMOUR_ONLY.out.table
                .map { meta, table ->
                    [meta, table, []]  // Empty file for matched normal
                }
            )

            // Run Mutect2 for tumor-only samples
            GATK4_MUTECT2_TUMOR_ONLY(
                unpaired_tumour_samples,
                unpaired_samples.map { 
                    meta, _tumor_bam, _tumor_bai, _normal_bam, _normal_bai ->
                    [meta, [], []]
                },
                fasta,
                fai,
                dict,
                germline_resource,
                germline_resource_tbi,
                panel_of_normals,
                panel_of_normals_tbi,
                intervals
            )

            GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_ONLY(
                GATK4_MUTECT2_TUMOR_ONLY.out.f1r2
            )
        }
        
        // Combine all Mutect2 outputs for further processing
        mutect2_read_orientation_models = GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_NORMAL.out.artifactprior
            .mix(GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_ONLY.out.artifactprior)

        mutect2_vcf = GATK4_MUTECT2_TUMOR_NORMAL.out.vcf.mix(GATK4_MUTECT2_TUMOR_ONLY.out.vcf)

        mutect2_tbi = GATK4_MUTECT2_TUMOR_NORMAL.out.tbi.mix(GATK4_MUTECT2_TUMOR_ONLY.out.tbi)

        mutect2_stats = GATK4_MUTECT2_TUMOR_NORMAL.out.stats.mix(GATK4_MUTECT2_TUMOR_ONLY.out.stats)

        // Combine contamination tables
        mutect2_contamination_tables = GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL.out.contamination
            .mix(GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY.out.contamination)
        
        // Combine segmentation tables
        mutect2_segmentation_tables = GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL.out.segmentation
            .mix(GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY.out.segmentation)
        
        // Prepare inputs for FilterMutectCalls
        filter_input = mutect2_vcf
            .join(mutect2_tbi, by: [0])
            .join(mutect2_stats, by: [0])
            .join(mutect2_read_orientation_models, by: [0])
            .join(mutect2_contamination_tables, by: [0])
            .join(mutect2_segmentation_tables, by: [0])
            .map { meta, vcf, tbi, stats, orientation_model, contamination, segmentation ->
                [meta, vcf, tbi, stats, orientation_model, contamination, segmentation]
            }
        
        filter_input.view()

        // Filter Mutect2 calls
        GATK4_FILTERMUTECTCALLS(
            filter_input,
            fasta,
            fai,
            dict
        )
        
        // Normalize variants with bcftools
        BCFTOOLS_NORM_MUTECT2(
            GATK4_FILTERMUTECTCALLS.out.vcf.join(GATK4_FILTERMUTECTCALLS.out.tbi, by: [0]),
            fasta,
            fai,
            dict
        )
        
        // Filter for PASS variants with bcftools
        BCFTOOLS_VIEW_MUTECT2(BCFTOOLS_NORM_MUTECT2.out.vcf)
        
        TABIX_MUTECT2(BCFTOOLS_VIEW_MUTECT2.out.vcf)

    // Define workflow outputs
    emit:
        vcf = BCFTOOLS_VIEW_MUTECT2.out.vcf
        tbi = TABIX_MUTECT2.out.tbi
}
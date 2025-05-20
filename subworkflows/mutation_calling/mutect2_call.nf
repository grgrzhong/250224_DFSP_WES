/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                GATK4 Mutect2 SNV/Indels somatic variant calling workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load required modules
include { GATK4_MUTECT2_TUMOUR as MUTECT2_TUMOUR                     } from '../../modules/local/gatk4/mutect2/tumour'
include { GATK4_MUTECT2_PAIRED as MUTECT2_PAIRED                     } from '../../modules/local/gatk4/mutect2/paired'
include { GATK4_GETPILEUPSUMMARIES_TUMOUR as PILEUP_PAIRED_TUMOUR    } from '../../modules/local/gatk4/getpileupsummaries/tumour'
include { GATK4_GETPILEUPSUMMARIES_NORMAL as PILEUP_PAIRED_NORMAL    } from '../../modules/local/gatk4/getpileupsummaries/normal'
include { GATK4_GETPILEUPSUMMARIES_TUMOUR as PILEUPS_UNPAIRED_TUMOUR } from '../../modules/local/gatk4/getpileupsummaries/tumour'
include { GATK4_CALCULATECONTAMINATION as CONTAMINATION_PAIRED       } from '../../modules/local/gatk4/calculatecontamination'
include { GATK4_CALCULATECONTAMINATION as CONTAMINATION_UNPAIRED     } from '../../modules/local/gatk4/calculatecontamination'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNMODEL_PAIRED       } from '../../modules/local/gatk4/learnreadorientationmodel'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNMODEL_UNPAIRED     } from '../../modules/local/gatk4/learnreadorientationmodel'
include { GATK4_FILTERMUTECTCALLS as FILTERMUTECTCALLS               } from '../../modules/local/gatk4/filtermutectcalls'
include { BCFTOOLS_NORM                                              } from '../../modules/local/bcftools/norm'
include { BCFTOOLS_VIEW                                              } from '../../modules/local/bcftools/view'
include { BCFTOOLS_ANNOTATE_REPEATMASKER as ANNOTATE_REPEATMASKER    } from '../../modules/local/bcftools/annotate/repeatmasker'
include { BCFTOOLS_ANNOTATE_BLACKLIST as ANNOTATE_BLACKLIST          } from '../../modules/local/bcftools/annotate/blacklist'
include { BCFTOOLS_FILTER as FILTER_REPEATMASKER_BLACKLIST           } from '../../modules/local/bcftools/filter'
include { ANNOVAR                                                    } from "../../modules/local/annovar/main.nf"
include { GATK4_FUNCOTATOR as FUNCOTATOR                             } from "../../modules/local/gatk4/funcotator/main.nf"

workflow MUTECT2_CALL {
    take:
        bam_tumour_normal
        bam_tumour_only
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
        repeatmasker
        blacklist
        funcotator_resources
        funcotator_ref_version
        annovar_db
        annovar_buildver
        annovar_protocol
        annovar_operation
        annovar_xreffile

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
     * =================== TUMOUR-NORMAL PAIRED ANALYSIS ================
     */
    

        // Extract tumour samples for pileup
        paired_tumour_samples = bam_tumour_normal
            .map { meta, tumour_bam, tumour_bai, _normal_bam, _normal_bai ->
                def new_meta = meta.clone()
                new_meta.id = meta.tumour_id
                [new_meta, tumour_bam, tumour_bai]
            }

        // Extract normal samples for pileup
        paired_normal_samples = bam_tumour_normal
            .map { meta, _tumour_bam, _tumour_bai, normal_bam, normal_bai ->
                def new_meta = meta.clone()
                new_meta.id = meta.normal_id
                    [new_meta, normal_bam, normal_bai]
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
        paired_tumour_normal_pileup = paired_tumour_pileup.table
            .map { meta, table -> [meta.patient_id, meta, table]}
            .join(
                paired_normal_pileup.table
                .map{
                    meta, table -> [meta.patient_id, meta, table]
                }, 
                remainder: true
            )
            .map {
                patient_id, tumour_meta, tumour_table, normal_meta, normal_table -> 
                [tumour_meta, tumour_table, normal_table]
            }
            
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
     * =================== TUMOUR-ONLY UNPAIRED ANALYSIS ====================
     */
    // [[id,patient_id,t_id,n_id,is_paired], t_bam, t_bai, n_bam, n_bai]

        unpaired_tumour_samples = bam_tumour_only
            .map { meta, tumour_bam, tumour_bai, _normal_bam, _normal_bai ->
                def new_meta = meta.clone()
                new_meta.id = meta.tumour_id
                [new_meta, tumour_bam, tumour_bai]
            }
        
        // Get pileup summaries for tumour-only samples
        unpaired_tumour_pileup = PILEUPS_UNPAIRED_TUMOUR(
            unpaired_tumour_samples,
            pileup_variants,
            pileup_variants_tbi
        )
        
        // Calculate contamination (no matched normal)
        unpaired_tumour_only_pileup = unpaired_tumour_pileup.table
            .map { meta, table -> [meta, table, []] }

        unpaired_contamination = CONTAMINATION_UNPAIRED(unpaired_tumour_only_pileup)

        
        // Run Mutect2 for tumour-only samples
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
        .join(mutect2_tbi)
        .join(mutect2_stats)
        .join(mutect2_orientation)
        .join(mutect2_contamination)
        .join(mutect2_segmentation)
        .map { meta, vcf, tbi, stats, orientation, contamination, segmentation ->
            [meta, vcf, tbi, stats, orientation, contamination, segmentation]
        }
    
    // filter_input.view()

    // Filter Mutect2 calls
    filtered_calls = FILTERMUTECTCALLS(
        filter_input,
        fasta,
        fai,
        dict
    )
    
    // Normalize variants with bcftools
    normalized_vcf = BCFTOOLS_NORM(
        filtered_calls.vcf.join(filtered_calls.tbi),
        fasta,
        fai,
        dict
    )
    
    // Filter for PASS variants with bcftools
    pass_vcf = BCFTOOLS_VIEW(
        normalized_vcf.vcf.join(normalized_vcf.tbi),
    )

    // Annotate with RepeatMasker and blacklist
    repeatmasker_vcf = ANNOTATE_REPEATMASKER(
        pass_vcf.vcf.join(pass_vcf.tbi),
        repeatmasker
    )

    blacklist_vcf = ANNOTATE_BLACKLIST(
        repeatmasker_vcf.vcf.join(repeatmasker_vcf.tbi),
        blacklist
    )

    // Filter out variants in RepeatMasker or Mapability
    final_vcf = FILTER_REPEATMASKER_BLACKLIST(
        blacklist_vcf.vcf.join(blacklist_vcf.tbi)
    )

    // Annotate variants with GATK Funcotator
    FUNCOTATOR(
        final_vcf.vcf.join(final_vcf.tbi),
        fasta,
        fai,
        dict,
        intervals,
        funcotator_resources,
        funcotator_ref_version
    )

    // Annotate variants with ANNOVAR
    ANNOVAR(
        final_vcf.vcf.join(final_vcf.tbi),
        annovar_db,
        annovar_buildver,
        annovar_protocol,
        annovar_operation,
        annovar_xreffile
    )

    emit:
    vcf = final_vcf.vcf              // channel: [ meta, vcf ]
    tbi = final_vcf.tbi             // channel: [ meta, tbi ]
}

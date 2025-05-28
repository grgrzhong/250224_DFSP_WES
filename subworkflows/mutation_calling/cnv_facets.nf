/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                CNV analysiws using FACETS subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CNV_FACETS_PAIRED as FACETS_PAIRED     } from '../../modules/local/facets/paired/main.nf'
include { CNV_FACETS_UNPAIRED as FACETS_UNPAIRED } from '../../modules/local/facets/unpaired/main.nf'

workflow CNV_FACETS {

    take:
        bam_tumour_normal
        bam_tumour_only
        dbsnp
        dbsnp_tbi
        defined_normal
        defined_normal_index

    main:

    // Initialize empty channels for results
    paired_cnv = Channel.empty()
    paired_cov = Channel.empty()
    paired_csv = Channel.empty()
    paired_spider = Channel.empty()
    paired_vcf = Channel.empty()
    paired_tbi = Channel.empty()
    
    unpaired_cnv = Channel.empty()
    unpaired_cov = Channel.empty()
    unpaired_csv = Channel.empty()
    unpaired_spider = Channel.empty()
    unpaired_vcf = Channel.empty()
    unpaired_tbi = Channel.empty()
    
    // Run FACETS for paired samples
    bam_tumour_normal
        .ifEmpty { 
            log.info "No paired samples found for FACETS analysis"
            Channel.empty()
        }
        .set { paired_input }

    // Only run FACETS_PAIRED if there are paired samples
    FACETS_PAIRED(
        paired_input, 
        dbsnp, 
        dbsnp_tbi
    )
    
    // Set paired results if process ran
    paired_cnv = FACETS_PAIRED.out.cnv
    paired_cov = FACETS_PAIRED.out.cov
    paired_csv = FACETS_PAIRED.out.csv
    paired_spider = FACETS_PAIRED.out.spider
    paired_vcf = FACETS_PAIRED.out.vcf
    paired_tbi = FACETS_PAIRED.out.vcf_tbi

    // Run FACETS for unpaired samples
    unpaired_ch = bam_tumour_only
        .ifEmpty { 
            log.info "No tumour-only samples found for FACETS analysis"
            Channel.empty()
        }
        .map { meta, tumour_bam, tumour_bai, normal_bam, normal_bai ->
            // Create proper tuple structure expected by FACETS_UNPAIRED
            [meta, tumour_bam, tumour_bai]
        }

    // Only run FACETS_UNPAIRED if there are unpaired samples
    FACETS_UNPAIRED(
        unpaired_ch, 
        defined_normal,
        defined_normal_index,
        dbsnp, 
        dbsnp_tbi
    )
    
    // Set unpaired results if process ran
    unpaired_cnv = FACETS_UNPAIRED.out.cnv
    unpaired_cov = FACETS_UNPAIRED.out.cov
    unpaired_csv = FACETS_UNPAIRED.out.csv
    unpaired_spider = FACETS_UNPAIRED.out.spider
    unpaired_vcf = FACETS_UNPAIRED.out.vcf
    unpaired_tbi = FACETS_UNPAIRED.out.vcf_tbi
    
    emit:
    // Mix results from both processes (will be empty if not run)
    cnv    = paired_cnv.mix(unpaired_cnv)
    cov    = paired_cov.mix(unpaired_cov)
    csv    = paired_csv.mix(unpaired_csv)
    spider = paired_spider.mix(unpaired_spider)
    vcf    = paired_vcf.mix(unpaired_vcf)
    tbi    = paired_tbi.mix(unpaired_tbi)

}
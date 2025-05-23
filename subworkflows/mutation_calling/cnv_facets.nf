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

    // Initialize the empty results channels
    paired_results = Channel.empty()
    unpaired_results = Channel.empty()
    
    // Run FACETS for paired samples if available
    FACETS_PAIRED(
        bam_tumour_normal, 
        dbsnp, 
        dbsnp_tbi
    )
    paired_results = FACETS_PAIRED.out


    // Run FACETS for unpaired samples if all requirements are met
    unpaired_ch = bam_tumour_only
        .map { meta, tumour_bam, tumour_bai, normal_bam, normal_bai ->
            // Create proper tuple structure expected by FACETS_UNPAIRED
            [meta, tumour_bam, tumour_bai]
        }

    FACETS_UNPAIRED(
        unpaired_ch, 
        defined_normal,
        defined_normal_index,
        dbsnp, 
        dbsnp_tbi
    )
    
    unpaired_results = FACETS_UNPAIRED.out

    
    // Create empty channels for unpaired results if not run
    unpaired_cnv    = unpaired_results.cnv ?: Channel.empty()
    unpaired_cov    = unpaired_results.cov ?: Channel.empty()
    unpaired_csv    = unpaired_results.csv ?: Channel.empty()
    unpaired_spider = unpaired_results.spider ?: Channel.empty()
    unpaired_vcf    = unpaired_results.vcf ?: Channel.empty()
    unpaired_tbi    = unpaired_results.vcf_tbi ?: Channel.empty()
    
    // Create empty channels for paired results if not run
    paired_cnv    = paired_results.cnv ?: Channel.empty()
    paired_cov    = paired_results.cov ?: Channel.empty()
    paired_csv    = paired_results.csv ?: Channel.empty()
    paired_spider = paired_results.spider ?: Channel.empty()
    paired_vcf    = paired_results.vcf ?: Channel.empty()
    paired_tbi    = paired_results.vcf_tbi ?: Channel.empty()
    
    emit:
    // Mix results from both processes (will be empty if not run)
    cnv    = paired_cnv.mix(unpaired_cnv)
    cov    = paired_cov.mix(unpaired_cov)
    csv    = paired_csv.mix(unpaired_csv)
    spider = paired_spider.mix(unpaired_spider)
    vcf    = paired_vcf.mix(unpaired_vcf)
    tbi    = paired_tbi.mix(unpaired_tbi)

}
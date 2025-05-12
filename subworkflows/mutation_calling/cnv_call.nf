/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                Copy number variation (CNV) calling subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Impport required modules
include { CNV_FACETS } from '../../modules/variant_calling/facets'

workflow CNV_CALL {
    
    take:
    bam_tumour_normal
    dbsnp
    dbsnp_tbi

    main:
    
    CNV_FACETS(
        bam_tumour_normal,
        dbsnp,
        dbsnp_tbi
    )

    emit:
    cnv_vcf = CNV_FACETS.out.vcf
    cnv_vcf_tbi = CNV_FACETS.out.vcf_tbi

}
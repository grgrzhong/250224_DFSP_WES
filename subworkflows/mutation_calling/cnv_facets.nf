/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                CNV analysiws using FACETS subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CNV_FACETS } from '../../modules/local/cnvfacets//main.nf'
include { VCF2TSVPY  } from '../../modules/local/vcf2tsvpy/main.nf'

workflow CNV_CALL_FACETS {

    take:
    bam_tumour_normal
    dbsnp
    dbsnp_tbi
    
    main:

    // Run FACETS for paired samples
    CNV_FACETS(
        bam_tumour_normal, 
        dbsnp, 
        dbsnp_tbi
    )

    // Convert the VCF output to TSV format
    VCF2TSVPY(CNV_FACETS.vcf)
    // Convert the vcf to tsv format
    emit:

    cnv    = CNV_FACETS.cnv
    cov    = CNV_FACETS.cov
    csv    = CNV_FACETS.csv
    tsv    = VCF2TSVPY.tsv
    spider = CNV_FACETS.spider
    vcf    = CNV_FACETS.vcf
    tbi    = CNV_FACETS.tbi

}
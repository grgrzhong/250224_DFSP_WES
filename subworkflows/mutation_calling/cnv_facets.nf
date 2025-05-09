
// Load required modules
include { CNV_FACETS } from '../../modules/variant_calling/facets/main.nf'

// Define workflow
workflow VARIANT_CALLING_CNV {
    
    take:
    paired_samples
    dbsnp
    dbsnp_tbi
    
    main:
    // Run CNV-FACETS on paired tumor-normal samples
    CNV_FACETS(
        paired_samples, 
        dbsnp, 
        dbsnp_tbi
    )
}
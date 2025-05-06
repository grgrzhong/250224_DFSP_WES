/*
    ====================================================================
    CNVkit subworkflow for variant calling
    ====================================================================
*/

// Load required modules
include { CNVKIT_REFERENCE } from '../../modules/variant_calling/cnvkit/reference/main'
include { CNVKIT_BATCH       } from '../../modules/variant_calling/cnvkit/batch/main'
include { CNVKIT_CALL        } from '../../modules/variant_calling/cnvkit/call/main'
include { CNVKIT_EXPORT      } from '../../modules/variant_calling/cnvkit/export/main'
include { CNVKIT_GENEMETRICS } from '../../modules/variant_calling/cnvkit/genemetrics/main'

workflow VARIANT_CALLING_CNVKIT {

    take:
    normal_bams         // channel: [optional]  meta, bam
    normal_bais        // channel: [optional]  meta, bai
    fasta               // channel: [mandatory] meta, fasta
    fai           // channel: [optional]  meta, fasta_fai
    dict
    targets             // channel: [mandatory] meta, bed

    main:
    // Create reference using all normal samples
    CNVKIT_REFERENCE(
        normal_bams,
        normal_bais,
        fasta,
        fai,
        dict,
        targets
    )

    emit:
    cnv_reference    = CNVKIT_REFERENCE.out.cnn
}

/*
    ====================================================================
    CNVkit subworkflow for variant calling
    ====================================================================
*/

// Load required modules
include { CNVKIT_REFERENCE   } from '../../modules/local/cnvkit/reference/main'
include { CNVKIT_BATCH       } from '../../modules/local/cnvkit/batch/main'
include { CNVKIT_CALL        } from '../../modules/local/cnvkit/call/main'
include { CNVKIT_EXPORT      } from '../../modules/local/cnvkit/export/main'
include { CNVKIT_GENEMETRICS } from '../../modules/local/cnvkit/genemetrics/main'

workflow CNV_CNVKIT {

    take:
    normal_bams
    normal_bais
    fasta
    fai
    dict
    targets

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

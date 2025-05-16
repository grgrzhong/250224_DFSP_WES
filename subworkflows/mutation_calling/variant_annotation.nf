
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                Variants annotation subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GATK4_FUNCOTATOR as FUNCOTATOR                } from "../../modules/variant_calling/gatk4/funcotator/main.nf"
include { GATK4_FUNCOTATOR_EXPORT as FUNCOTATOR_EXPORT  } from "../../modules/variant_calling/gatk4/funcotator/export/main.nf"
include { ANNOVAR                                       } from "../../modules/variant_calling/annovar/main.nf"
include { ANNOVAR_EXPORT                                } from "../../modules/variant_calling/annovar/export/main.nf"

workflow VARIANT_ANNOTATION {

    take:
    input_vcf
    fasta
    fai
    dict
    intervals
    funcotator_resources
    funcotator_ref_version
    annovar_db
    annovar_buildver
    annovar_protocol
    annovar_operation
    annovar_xreffile
    
    main:
    
    // Annotate variants with GATK Funcotator

    FUNCOTATOR(
        input_vcf,
        fasta,
        fai,
        dict,
        intervals,
        funcotator_resources,
        funcotator_ref_version
    )

    // Annotate variants with ANNOVAR
    ANNOVAR(
        input_vcf,
        annovar_db,
        annovar_buildver,
        annovar_protocol,
        annovar_operation,
        annovar_xreffile
    )

    emit:
    funcotator_maf          = FUNCOTATOR.out.maf
    funcotator_tsv          = FUNCOTATOR.out.tsv
    annovar_txt             = ANNOVAR.out.annovar
}
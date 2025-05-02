
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                Variants annotation workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GATK4_FUNCOTATOR } from '../../modules/variant_calling/gatk4/funcotator'
include { ANNOVAR          } from '../../modules/variant_calling/annovar'

workflow VARIANT_ANNOTATION {

    take:
        vcf             // channel: [val(meta), path(vcf)]
        vcf_tbi         // channel: [val(meta), path(tbi)]
        fasta           // path: Reference genome fasta
        fai             // path: Reference genome index 
        dict            // path: Reference genome dictionary
        funcotator_resources // path: GATK Funcotator resources
        intervals       // path: Intervals file (optional)
        annovar_db      // path: ANNOVAR database
    
    main:
        // Annotate variants with GATK Funcotator
        GATK4_FUNCOTATOR(
            BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi, by: [0]),
            fasta,
            fai,
            dict,
            funcotator_resources,
            intervals
        )
        
        // Annotate variants with ANNOVAR
        ANNOVAR(
            BCFTOOLS_VIEW.out.vcf,
            annovar_db
        )
    
    emit:
        funcotator_vcf = GATK4_FUNCOTATOR.out.vcf    // channel: [val(meta), path(vcf)]
        funcotator_html = GATK4_FUNCOTATOR.out.html  // channel: [val(meta), path(html)] 
        annovar_txt = ANNOVAR.out.txt                // channel: [val(meta), path(txt)]
}
    
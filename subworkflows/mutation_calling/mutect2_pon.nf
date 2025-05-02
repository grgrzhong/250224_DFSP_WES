
include { GATK4_MUTECT2 as GATK4_MUTECT2_NORMAL } from '../../modules/variant_calling/gatk4/mutect2/'
include { GATK4_GENOMICSDBIMPORT } from '../../modules/variant_calling/gatk4/genomicsdbimport/'
include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../modules/variant_calling/gatk4/createsomaticpanelofnormals/'

workflow MUTECT2_SOMATIC_PON {
    

    main:
       // Prepare reference files
        reference_fai = Channel.fromPath("${reference}.fai")
        reference_dict = Channel.fromPath("${reference.parent}/${reference.baseName}.dict")
        
        // Run Mutect2 on each normal sample
        GATK4_MUTECT2_NORMAL(
            normal_samples,
            reference,
            reference_fai,
            reference_dict,
            interval
        )
        
        // Collect all VCFs for GenomicsDB Import
        vcfs = GATK4_MUTECT2_NORMAL.out.vcfs
            .map { sample_id, vcf, tbi -> vcf }
            .collect()
            
        tbis = GATK4_MUTECT2_NORMAL.out.vcfs
            .map { sample_id, vcf, tbi -> tbi }
            .collect()
            
        // Create GenomicsDB
        GATK4_GENOMICSDBIMPORT(
            vcfs,
            tbis,
            reference,
            reference_fai,
            reference_dict,
            interval
        )
        
        // Create Panel of Normals
        GATK4_CREATESOMATICPANELOFNORMALS(
            GATK4_GENOMICSDBIMPORT.out.pon_db,
            reference,
            reference_fai,
            reference_dict,
            interval
        )
        
    emit:
        pon_vcf = CREATE_SOMATIC_PON.out.pon_vcf
        pon_vcf_tbi = CREATE_SOMATIC_PON.out.pon_vcf_tbi
    
}
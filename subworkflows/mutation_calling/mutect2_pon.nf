
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                GATK4 Mutect2 create PON subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load the required modules
include { GATK4_MUTECT2_NORMAL              } from '../../modules/local/gatk4/mutect2/normal'
include { GATK4_GENOMICSDBIMPORT            } from '../../modules/local/gatk4/genomicsdbimport'
include { GATK4_CREATESOMATICPANELOFNORMALS } from '../../modules/local/gatk4/createsomaticpanelofnormals'

workflow MUTECT2_PON {
    
    take:
        normal_samples
        fasta
        fai
        dict
        interval

    main: 
        // Run Mutect2 on each normal sample
        GATK4_MUTECT2_NORMAL(
            normal_samples,
            fasta,
            fai,
            dict,
            interval
        )
        
        // Collect all VCFs for GenomicsDB Import
        vcfs = GATK4_MUTECT2_NORMAL.out.vcf.collect()
        tbis = GATK4_MUTECT2_NORMAL.out.tbi.collect()
        
            
        // Create GenomicsDB
        GATK4_GENOMICSDBIMPORT(
            vcfs,
            tbis,
            fasta,
            fai,
            dict,
            interval
        )
        
        // Create Panel of Normals
        GATK4_CREATESOMATICPANELOFNORMALS(
            GATK4_GENOMICSDBIMPORT.out.pon_db,
            fasta,
            fai,
            dict,
            interval
        )
        
    emit:
        pon_vcf = GATK4_CREATESOMATICPANELOFNORMALS.out.vcf
        pon_vcf_tbi = GATK4_CREATESOMATICPANELOFNORMALS.out.tbi
    
}
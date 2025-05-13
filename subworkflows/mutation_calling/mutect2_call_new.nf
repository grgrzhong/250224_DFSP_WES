/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                GATK4 Mutect2 SNV/Indels somatic variant calling workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load required modules
include { GATK4_MUTECT2_TUMOUR as MUTECT2_TUMOUR                     } from '../../modules/variant_calling/gatk4/mutect2/tumour'
include { GATK4_MUTECT2_PAIRED as MUTECT2_PAIRED                     } from '../../modules/variant_calling/gatk4/mutect2/paired'

params.fasta                    = params.genomes[params.genome]?.fasta
params.fai                      = params.genomes[params.genome]?.fai
params.dict                     = params.genomes[params.genome]?.dict
params.dbsnp                    = params.genomes[params.genome]?.dbsnp
params.dbsnp_tbi                = params.genomes[params.genome]?.dbsnp_tbi
params.germline_resource        = params.genomes[params.genome]?.germline_resource
params.germline_resource_tbi    = params.genomes[params.genome]?.germline_resource_tbi
params.panel_of_normals         = params.genomes[params.genome]?.pon
params.panel_of_normals_tbi     = params.genomes[params.genome]?.pon_tbi
params.pileup_variants          = params.genomes[params.genome]?.contamination_variants
params.pileup_variants_tbi      = params.genomes[params.genome]?.contamination_variants_tbi
params.intervals                = params.genomes[params.genome]?.intervals
params.bait_intervals           = params.genomes[params.genome]?.bait_intervals
params.target_intervals         = params.genomes[params.genome]?.target_intervals
params.targets                  = params.genomes[params.genome]?.targets
params.funcotator_resources     = params.genomes[params.genome]?.funcotator
params.annovar_db               = params.genomes[params.genome]?.annovar_db

params.repeatmasker           = params.genomes[params.genome]?.repeatmasker
params.repeatmasker_tbi       = params.genomes[params.genome]?.repeatmasker_tbi
params.blacklist              = params.genomes[params.genome]?.blacklist
params.blacklist_tbi          = params.genomes[params.genome]?.blacklist_tbi

workflow MUTECT2_CALL {
    take:
        fasta
        fai
        dict
        germline_resource
        germline_resource_tbi
        panel_of_normals
        panel_of_normals_tbi
        intervals
        bam_tumour_normal
        bam_tumour_only

    main:
 // Initialize empty channels for results
    tumour_normal_vcf = Channel.empty()
    tumour_normal_tbi = Channel.empty()
    tumour_only_vcf = Channel.empty()
    tumour_only_tbi = Channel.empty()

    /*
     * =================== TUMOUR-NORMAL PAIRED ANALYSIS ================
     */
    
    bam_tumour_normal_count = bam_tumour_normal.count()
    
    
    // Extract tumour samples for pileup
    paired_tumour_samples = bam_tumour_normal
        .map { meta, tumour_bam, tumour_bai, _normal_bam, _normal_bai ->
            def new_meta = meta.clone()
            new_meta.id = meta.tumour_id
            [new_meta, tumour_bam, tumour_bai]
        }
    
    // Extract normal samples for pileup
    paired_normal_samples = bam_tumour_normal
        .map { meta, _tumour_bam, _tumour_bai, normal_bam, normal_bai ->
            def new_meta = meta.clone()
            new_meta.id = meta.normal_id
                [new_meta, normal_bam, normal_bai]
        }
    
    // Run Mutect2
    paired_mutect2 = MUTECT2_PAIRED(
        paired_tumour_samples,
        paired_normal_samples,
        fasta,
        fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi,
        intervals
    )

    // Store output channels
    tumour_normal_vcf = paired_mutect2.vcf
    tumour_normal_tbi = paired_mutect2.tbi

    /*
     * =================== TUMOUR-ONLY UNPAIRED ANALYSIS ====================
     */
    // [[id,patient_id,t_id,n_id,is_paired], t_bam, t_bai, n_bam, n_bai]
    bam_tumour_only_count = bam_tumour_only.count()

    
    // Extract tumour samples for pileup
    // [[id,patient_id,t_id,n_id,is_paired], t_bam, t_bai, n_bam, n_bai]
    unpaired_tumour_samples = bam_tumour_only
        .map { meta, tumour_bam, tumour_bai, _normal_bam, _normal_bai ->
            def new_meta = meta.clone()
            new_meta.id = meta.tumour_id
            [new_meta, tumour_bam, tumour_bai]
        }

    // Run Mutect2 for tumour-only samples
    unpaired_mutect2 = MUTECT2_TUMOUR(
        unpaired_tumour_samples,
        fasta,
        fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi,
        intervals
    )

    // Store output channels
    tumour_only_vcf = unpaired_mutect2.vcf
    tumour_only_tbi = unpaired_mutect2.tbi

    

    // Combine all Mutect2 outputs for further processing
    mutect2_vcf = tumour_normal_vcf.mix(tumour_only_vcf)
    mutect2_tbi = tumour_normal_tbi.mix(tumour_only_tbi)
    
    // Prepare inputs for FilterMutectCalls
    output = mutect2_vcf
        .join(mutect2_tbi)
        .map { meta, vcf, tbi ->
            [meta, vcf, tbi]
        }

    emit:
    vcf = output                // channel: [ meta, vcf, tbi ]
}

include { PREPARE_SAMPLE } from '../../subworkflows/mutation_calling/prepare_sample.nf'

workflow {
    
    fasta                   = params.fasta ? file(params.fasta) : null
    fai                     = params.fai ? file(params.fai) : null
    dict                    = params.dict ? file(params.dict) : null
    germline_resource       = params.germline_resource ? file(params.germline_resource) : null
    germline_resource_tbi   = params.germline_resource_tbi ? file(params.germline_resource_tbi) : null
    panel_of_normals        = params.panel_of_normals ? file(params.panel_of_normals) : null
    panel_of_normals_tbi    = params.panel_of_normals_tbi ? file(params.panel_of_normals_tbi) : null
    pileup_variants         = params.pileup_variants ? file(params.pileup_variants) : null
    pileup_variants_tbi     = params.pileup_variants_tbi ? file(params.pileup_variants_tbi) : null
    intervals               = params.intervals ? file(params.intervals) : null
    bait_intervals          = params.bait_intervals ? file(params.bait_intervals) : null
    target_intervals        = params.target_intervals ? file(params.target_intervals) : null
    targets                 = params.targets ? file(params.targets) : null
    funcotator_resources    = params.funcotator_resources ? file(params.funcotator_resources) : null
    annovar_db              = params.annovar_db ? file(params.annovar_db) : null
    repeatmasker            = params.repeatmasker ? file(params.repeatmasker) : null
    repeatmasker_tbi        = params.repeatmasker_tbi ? file(params.repeatmasker_tbi) : null
    blacklist               = params.blacklist ? file(params.blacklist) : null
    blacklist_tbi           = params.blacklist_tbi ? file(params.blacklist_tbi) : null

    samples = PREPARE_SAMPLE(params.input)
    
    MUTECT2_CALL(
        fasta,
        fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi,
        intervals,
        samples.bam_tumour_normal,
        samples.bam_tumour_only
    )
}
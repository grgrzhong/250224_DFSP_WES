#!/usr/bin/env nextflow

// Import modules and subworkflows
include { PREPARE_SAMPLE     } from "../subworkflows/mutation_calling/prepare_sample.nf"
include { PREPROCESSING      } from "../subworkflows/mutation_calling/preprocessing.nf"
include { MUTECT2_CALL       } from "../subworkflows/mutation_calling/mutect2_call.nf"
include { CNV_SEQUENZA       } from "../subworkflows/mutation_calling/cnv_sequenza.nf"
include { CNV_FACETS         } from "../subworkflows/mutation_calling/cnv_facets.nf"

// Input csv file
params.input = "/home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test2.csv"

//  Reference genome and resources
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
params.repeatmasker             = params.genomes[params.genome]?.repeatmasker
params.blacklist                = params.genomes[params.genome]?.blacklist
params.funcotator_resources     = params.genomes[params.genome]?.funcotator_resources
params.annovar_db               = params.genomes[params.genome]?.annovar_db
params.annovar_buildver         = params.genomes[params.genome]?.annovar_buildver
params.annovar_protocol         = params.genomes[params.genome]?.annovar_protocol
params.annovar_operation        = params.genomes[params.genome]?.annovar_operation
params.annovar_xreffile         = params.genomes[params.genome]?.annovar_xreffile
params.funcotator_ref_version   = params.genomes[params.genome]?.funcotator_ref_version
params.wigfile                  = params.genomes[params.genome]?.wigfile
params.window_size              = params.genomes[params.genome]?.window_size

// Main workflow
workflow {

    /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Prepare the genome and resource files
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    fasta                   = params.fasta ? file(params.fasta) : null
    fai                     = params.fai ? file(params.fai) : null
    dict                    = params.dict ? file(params.dict) : null
    dbsnp                   = params.dbsnp ? file(params.dbsnp) : null
    dbsnp_tbi               = params.dbsnp_tbi ? file(params.dbsnp_tbi) : null
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
    repeatmasker            = params.repeatmasker ? file(params.repeatmasker) : null
    blacklist               = params.blacklist ? file(params.blacklist) : null

    funcotator_resources    = params.funcotator_resources ? file(params.funcotator_resources) : null
    funcotator_ref_version  = params.funcotator_ref_version

    annovar_db              = params.annovar_db ? file(params.annovar_db) : null
    annovar_buildver        = params.annovar_buildver
    annovar_protocol        = params.annovar_protocol
    annovar_operation       = params.annovar_operation
    annovar_xreffile        = params.annovar_xreffile ? file(params.annovar_xreffile) : null
    
    wigfile                 = params.wigfile ? file(params.wigfile) : null
    window_size             = params.window_size

    // Manuall defined normal sample for the toumour sample without normal
    defined_normal          = params.defined_normal ? file(params.defined_normal) : null
    defined_normal_index    = params.defined_normal_index ? file(params.defined_normal_index) : null

    /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        Print Genome and Resource File Information
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    log.info "================== Genome & Resources =================="
    log.info "Reference genome           = ${fasta}"
    log.info "FASTA index                = ${fai}"
    log.info "Dictionary                 = ${dict}"
    log.info "dbSNP                      = ${dbsnp}"
    log.info "dbSNP index                = ${dbsnp_tbi}"
    log.info "Germline resource          = ${germline_resource}"
    log.info "Germline resource index    = ${germline_resource_tbi}"
    log.info "Panel of normals           = ${panel_of_normals}"
    log.info "Panel of normals index     = ${panel_of_normals_tbi}"
    log.info "Pileup variants            = ${pileup_variants}"
    log.info "Pileup variants index      = ${pileup_variants_tbi}"
    log.info "Intervals file             = ${intervals}"
    log.info "Bait intervals             = ${bait_intervals}"
    log.info "Target intervals           = ${target_intervals}"
    log.info "Targets                    = ${targets}"
    log.info "RepeatMasker               = ${repeatmasker}"
    log.info "Blacklist                  = ${blacklist}"
    log.info "Funcotator resources       = ${funcotator_resources}"
    log.info "ANNOVAR resources          = ${annovar_db}"
    log.info "Wigfile                    = ${wigfile}"
    log.info "Definied normal            = ${defined_normal}"

    /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        Prepare the samples and metadata
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    samples = PREPARE_SAMPLE(params.input)
    
    input_samples       = samples.input_samples
    tumour_samples      = samples.tumour_samples
    normal_samples      = samples.normal_samples
    fastq               = samples.fastq
    bam_tumour_normal   = samples.bam_tumour_normal
    bam_tumour_only     = samples.bam_tumour_only
    bam_all_samples     = samples.bam_all_samples

    /*
        ======================================================================
                        Run the full workflow
        ======================================================================
    */
    // // Preprocessing fastq files and bwa mapping
    // PREPROCESSING(
    //     fastq,
    //     fasta,
    //     fai,
    //     dict,
    //     dbsnp,
    //     dbsnp_tbi,
    //     bait_intervals,
    //     target_intervals,
    //     intervals
    // )

    // Run Mutect2 SNV/indels variant calling and annotation
    // MUTECT2_CALL(
    //     bam_tumour_normal,
    //     bam_tumour_only,
    //     fasta,
    //     fai,
    //     dict,
    //     pileup_variants,
    //     pileup_variants_tbi,
    //     germline_resource,
    //     germline_resource_tbi,
    //     panel_of_normals,
    //     panel_of_normals_tbi,
    //     intervals,
    //     repeatmasker,
    //     blacklist,
    //     funcotator_resources,
    //     funcotator_ref_version,
    //     annovar_db,
    //     annovar_buildver,
    //     annovar_protocol,
    //     annovar_operation,
    //     annovar_xreffile
    // )

    // annotation_input = MUTECT2_CALL.out.vcf.join(MUTECT2_CALL.out.tbi)

    // Run CNV analysis using FACETS
    CNV_FACETS(
        bam_tumour_normal,
        bam_tumour_only,
        dbsnp,
        dbsnp_tbi,
        defined_normal,
        defined_normal_index
    )

    // Run CNV analysis using Sequenza
    // CNV_SEQUENZA(
    //     bam_tumour_normal,
    //     fasta,
    //     wigfile,
    //     window_size
    // )

}

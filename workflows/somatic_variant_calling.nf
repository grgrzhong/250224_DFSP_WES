#!/usr/bin/env nextflow

// Import modules and subworkflows
include { CNV_FACETS         } from '../modules/variant_calling/facets/main.nf'
include { PREPARE_SAMPLE     } from "../subworkflows/mutation_calling/prepare_sample.nf"
include { PREPROCESSING      } from "../subworkflows/mutation_calling/preprocessing.nf"
include { MUTECT2_CALL       } from "../subworkflows/mutation_calling/mutect2_call.nf"
include { VARIANT_ANNOTATION } from "../subworkflows/mutation_calling/variant_annotation.nf"

// Define workflow steps and their dependencies as parameters
params.step_options = [
    "preprocessing",
    "mutect2_call",
    "variant_annotation",
    "cnv_facets"
]

// Default step (run full workflow)
params.step = null

// Define valid starting and ending points
params.valid_steps = [
    "preprocessing",
    "mutect2_call",
    "variant_annotation"
]

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
params.funcotator_ref_version   =params.genomes[params.genome]?.funcotator_ref_version

// Parameter validation for steps
def validate_params() {
    if (params.step) {
        if (!params.valid_steps.contains(params.step)) {
            log.error("Invalid starting step '${params.step}'. Valid options are: ${params.valid_steps.join(', ')}")
            workflow.exit(1)
        }
    }
}

// Main workflow
workflow {
    
    validate_params()
    
    if (params.step) {
        log.info "Starting workflow from step: ${params.step}"
    } else {
        log.info "Running full workflow"
    }

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
    funcotator_ref_version   =params.funcotator_ref_version

    annovar_db              = params.annovar_db ? file(params.annovar_db) : null
    annovar_buildver        = params.annovar_buildver
    annovar_protocol        = params.annovar_protocol
    annovar_operation       = params.annovar_operation
    annovar_xreffile        = params.annovar_xreffile ? file(params.annovar_xreffile) : null
    
    /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        Print Genome and Resource File Information
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    log.info "======= GENOME AND RESOURCES ======="
    log.info "Reference genome:          ${fasta}"
    log.info "FASTA index:               ${fai}"
    log.info "Dictionary:                ${dict}"
    log.info "dbSNP:                     ${dbsnp}"
    log.info "dbSNP index:               ${dbsnp_tbi}"
    log.info "Germline resource:         ${germline_resource}"
    log.info "Germline resource index:   ${germline_resource_tbi}"
    log.info "Panel of normals:          ${panel_of_normals}"
    log.info "Panel of normals index:    ${panel_of_normals_tbi}"
    log.info "Pileup variants:           ${pileup_variants}"
    log.info "Pileup variants index:     ${pileup_variants_tbi}"
    log.info "Intervals file:            ${intervals}"
    log.info "Bait intervals:            ${bait_intervals}"
    log.info "Target intervals:          ${target_intervals}"
    log.info "Targets:                   ${targets}"
    log.info "RepeatMasker:              ${repeatmasker}"
    log.info "Blacklist:                 ${blacklist}"
    log.info "Funcotator resources:      ${funcotator_resources}"
    log.info "ANNOVAR database:          ${annovar_db}"
    log.info "ANNOVAR build version:     ${annovar_buildver}"
    log.info "ANNOVAR protocol:          ${annovar_protocol}"
    log.info "ANNOVAR operation:         ${annovar_operation}"
    log.info "ANNOVAR xref file:         ${annovar_xreffile}"
    log.info "===================================="

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
    // Preprocessing fastq files and bwa mapping
    PREPROCESSING(
        fastq,
        fasta,
        fai,
        dict,
        dbsnp,
        dbsnp_tbi,
        bait_intervals,
        target_intervals,
        intervals
    )

    // Run Mutect2 variant calling
    MUTECT2_CALL(
        bam_tumour_normal,
        bam_tumour_only,
        fasta,
        fai,
        dict,
        pileup_variants,
        pileup_variants_tbi,
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi,
        intervals,
        repeatmasker,
        blacklist
    )

    annotation_input = MUTECT2_CALL.out.vcf.join(MUTECT2_CALL.out.tbi)

    // Run variant annotation
    VARIANT_ANNOTATION(
        annotation_input,
        fasta,
        fai,
        dict,
        intervals,
        funcotator_resources,
        funcotator_ref_version,
        annovar_db,
        annovar_buildver,
        annovar_protocol,
        annovar_operation,
        annovar_xreffile
    )
    
    // Run CNV analysis using FACETS
    CNV_FACETS(
        bam_tumour_normal,
        dbsnp,
        dbsnp_tbi
    )

    /*
        ======================================================================
                        Define steps to execute 
        ======================================================================
    */
    // validate_params()

    // // Define step execution boolean flags - simplified for your workflow
    // run_preprocessing = !params.step || params.step == 'preprocessing'
    // run_mutect2 = !params.step || params.step == 'mutect2'
    // run_cnv_facets = !params.step || params.step == 'cnv_facets'

    // // Variables to hold intermediate results
    // def preprocessed_results = null

    // // Execute workflow based on the step parameter
    // if (run_preprocessing) {

    //     log.info("Running preprocessing step")

    //     PREPROCESSING(
    //         fasta,
    //         fai,
    //         dict,
    //         dbsnp,
    //         dbsnp_tbi,
    //         bait_intervals,
    //         target_intervals,
    //         intervals,
    //         reads,
    //     )

    //     // Store preprocessing results for potential use in next step
    //     preprocessed_results = PREPROCESSING.out
    // }

    // if (run_mutect2) {

    //     log.info("Running Mutect2 variant calling step")

    //     // Determine the input for Mutect2
    //     def mutect2_input

    //     if (params.step == 'mutect2') {
    //         // Starting from mutect2 step directly - use sample sheet info
    //         log.info("Using sample sheet data for Mutect2 input")
    //     }
    //     else {
    //         // Coming from preprocessing step - use its output
    //         log.info("Using preprocessing output for Mutect2 input")
    //     }

    //     // Run Mutect2 somatic variant calling
    //     MUTECT2_CALL(
    //         fasta,
    //         fai,
    //         dict,
    //         pileup_variants,
    //         pileup_variants_tbi,
    //         germline_resource,
    //         germline_resource_tbi,
    //         panel_of_normals,
    //         panel_of_normals_tbi,
    //         intervals,
    //         paired_samples,
    //         unpaired_samples,
    //     )

    //     // Capture Mutect2 outputs for potential downstream use
    //     mutect2_vcf = MUTECT2_CALL.out.vcf
    //     mutect2_tbi = MUTECT2_CALL.out.tbi
    // }

    // Future steps could be added here
    // For example, VCF annotation, reporting, etc.

    // if (run_cnvkit) {

    //     log.info("Running CNVkit variant calling step")

    //     // Run CNVkit variant calling
    //     normal_bams = paired_normal_samples.map { _meta, normal_bam, _normal_bai ->
    //         normal_bam
    //     }

    //     normal_bais = paired_normal_samples.map { _meta, _normal_bam, normal_bai ->
    //         normal_bai
    //     }

    //     VARIANT_CALLING_CNVKIT(
    //         normal_bams,
    //         normal_bais,
    //         fasta,
    //         fai,
    //         dict,
    //         targets,
    //     )

    //     // Capture CNVkit outputs for potential downstream use
    //     cnv_reference = VARIANT_CALLING_CNVKIT.out.cnv_reference
    // }


    workflow.onComplete = {
        log.info("Pipeline completed at ${workflow.complete}")
        // Consider adding more detailed summary if needed, e.g., workflow.success
        if (workflow.success) {
            log.info("Pipeline completed successfully.")
        } else {
            log.error("Pipeline completed with errors.")
        }
    }
}

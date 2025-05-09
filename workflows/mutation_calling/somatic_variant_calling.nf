#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// input
params.input        = "/home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test2.csv"
// params.input = "/home/zhonggr/projects/250224_DFSP_WES/data/test_data/csv/test.csv"

// test the preprocesing
// params.input            = params.input ?: "/home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test.csv"


// Define workflow steps and their dependencies as parameters
params.step_options = [
    'preprocessing',
    'mutect2',
    "cnvkit",
    "cnv_facets"
]

// Default step (run full workflow)
params.step = null

// Define valid starting and ending points
params.valid_steps  = [
    'preprocessing',
    'mutect2',
    "cnvkit",
    "cnv_facets"
]

// Import modules and subworkflows
include { PREPARE_GENOME            } from "../../subworkflows/mutation_calling/prepare_genome.nf"
include { PREPARE_SAMPLE            } from "../../subworkflows/mutation_calling/prepare_sample.nf"
include { PREPROCESSING             } from "../../subworkflows/mutation_calling/preprocessing.nf"
include { MUTECT2_SOMATIC_PON       } from '../../subworkflows/mutation_calling/mutect2_pon.nf'
include { MUTECT2_CALL              } from "../../subworkflows/mutation_calling/mutect2_call.nf"
include { VARIANT_CALLING_CNVKIT    } from "../../subworkflows/mutation_calling/cnvkit.nf"
include { VARIANT_CALLING_CNV } from '../../subworkflows/mutation_calling/cnv_facets.nf'

// Parameter validation for steps
def validate_params() {
    // Check if step is valid
    if (params.step) {
        if (!params.valid_steps.contains(params.step)) {
            log.error("Invalid step '${params.step}'. Valid options are: ${params.valid_steps.join(', ')}")
            System.exit(1)
        }
    }
}

// Main workflow logic
workflow {

    /*
        ======================================================================
                    Prepare the genome and resource files
        ======================================================================
    */

    PREPARE_GENOME(params.genome)

    // Extract reference channels from PREPARE_GENOME
    fasta = PREPARE_GENOME.out.fasta
    fai = PREPARE_GENOME.out.fai
    dict = PREPARE_GENOME.out.dict
    germline_resource = PREPARE_GENOME.out.germline_resource
    germline_resource_tbi = PREPARE_GENOME.out.germline_resource_tbi
    panel_of_normals = PREPARE_GENOME.out.panel_of_normals
    panel_of_normals_tbi = PREPARE_GENOME.out.panel_of_normals_tbi
    pileup_variants = PREPARE_GENOME.out.pileup_variants
    pileup_variants_tbi = PREPARE_GENOME.out.pileup_variants_tbi
    dbsnp = PREPARE_GENOME.out.dbsnp
    dbsnp_tbi = PREPARE_GENOME.out.dbsnp_tbi
    intervals = PREPARE_GENOME.out.intervals
    bait_intervals = PREPARE_GENOME.out.bait_intervals
    target_intervals = PREPARE_GENOME.out.target_intervals
    targets = PREPARE_GENOME.out.targets

    /*
        ======================================================================
                    Prepare the input samples and their metadata
        ======================================================================
    */
    if (params.input == null) {
        error("Please provide an input CSV file with --input")
    }

    PREPARE_SAMPLE(params.input)
    reads                   = PREPARE_SAMPLE.out.reads
    paired_samples          = PREPARE_SAMPLE.out.paired_samples
    unpaired_samples        = PREPARE_SAMPLE.out.unpaired_samples
    paired_tumour_samples   = PREPARE_SAMPLE.out.paired_tumour_samples
    paired_normal_samples   = PREPARE_SAMPLE.out.paired_normal_samples
    unpaired_tumour_samples = PREPARE_SAMPLE.out.unpaired_tumour_samples
    all_samples             = PREPARE_SAMPLE.out.all_samples

    /*
        ======================================================================
                    Define steps to execute based on the workflow 
        ======================================================================
    */
    validate_params()

    // Define step execution boolean flags - simplified for your workflow
    run_preprocessing = !params.step || params.step == 'preprocessing'
    run_mutect2 = !params.step || params.step == 'mutect2'
    run_cnvkit = !params.step || params.step == 'cnvkit'
    run_cnv_facets = !params.step || params.step == 'cnv_facets'

    // Variables to hold intermediate results
    def preprocessed_results = null

    // Execute workflow based on the step parameter
    if (run_preprocessing) {

        log.info("Running preprocessing step")

        PREPROCESSING(
            fasta,
            fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            bait_intervals,
            target_intervals,
            intervals,
            reads,
        )

        // Store preprocessing results for potential use in next step
        preprocessed_results = PREPROCESSING.out
    }

    if (run_mutect2) {

        log.info("Running Mutect2 variant calling step")

        // Determine the input for Mutect2
        def mutect2_input

        if (params.step == 'mutect2') {
            // Starting from mutect2 step directly - use sample sheet info
            log.info("Using sample sheet data for Mutect2 input")
        }
        else {
            // Coming from preprocessing step - use its output
            log.info("Using preprocessing output for Mutect2 input")
        }

        // Run Mutect2 somatic variant calling
        MUTECT2_CALL(
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
            paired_samples,
            unpaired_samples,
        )

        // Capture Mutect2 outputs for potential downstream use
        mutect2_vcf = MUTECT2_CALL.out.vcf
        mutect2_tbi = MUTECT2_CALL.out.tbi
    }

    // Future steps could be added here
    // For example, VCF annotation, reporting, etc.

    if (run_cnvkit) {

        log.info("Running CNVkit variant calling step")

        // Run CNVkit variant calling
        normal_bams = paired_normal_samples.map(
            { _meta, normal_bam, _normal_bai -> 
                normal_bam
                }
        )
        
        normal_bais = paired_normal_samples.map(
            { _meta, _normal_bam, normal_bai -> 
                normal_bai
                }
        )
        
        VARIANT_CALLING_CNVKIT(
            normal_bams,
            normal_bais,
            fasta,
            fai,
            dict,
            targets
        )

        // Capture CNVkit outputs for potential downstream use
        cnv_reference = VARIANT_CALLING_CNVKIT.out.cnv_reference
    }

    if(run_cnv_facets) {

        log.info("Running CNV-FACETS variant calling")

        VARIANT_CALLING_CNV(
            paired_samples,
            dbsnp,
            dbsnp_tbi
        )
    }

    workflow.onComplete = {

        log.info("Pipeline completed at ${workflow.complete}")
    }
}

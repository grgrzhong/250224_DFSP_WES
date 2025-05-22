#!/usr/bin/env nextflow

// Import modules and subworkflows
include { PREPARE_SAMPLE     } from "../subworkflows/mutation_calling/prepare_sample.nf"
include { PREPROCESSING      } from "../subworkflows/mutation_calling/preprocessing.nf"
include { MUTECT2_CALL       } from "../subworkflows/mutation_calling/mutect2_call.nf"
include { CNV_SEQUENZA       } from "../subworkflows/mutation_calling/cnv_sequenza.nf"
include { CNV_FACETS         } from "../subworkflows/mutation_calling/cnv_facets.nf"

// Input csv file
params.input = "/home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test2.csv"

// Workflow step parameters - control which steps to run
params.step = null // Define step parameter, similar to nf-core/sarek

// Define valid workflow steps
params.validSteps = ['preprocessing', 'mutect2', 'facets', 'sequenza']

// Set default values for individual step parameters
params.run_preprocessing = false
params.run_mutect2 = false 
params.run_facets = false
params.run_sequenza = false

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

    // Process the step parameter to determine which steps to run
    // Create local variables for step control - these are not params
    def run_preprocessing = false
    def run_mutect2 = false
    def run_facets = false
    def run_sequenza = false
    
    if (params.step) {
        if (!params.validSteps.contains(params.step)) {
            exit 1, "Invalid step: '${params.step}'. Valid steps are: ${params.validSteps.join(', ')}"
        }
        
        // Enable only the selected step
        run_preprocessing = params.step == 'preprocessing'
        run_mutect2 = params.step == 'mutect2'
        run_facets = params.step == 'facets'
        run_sequenza = params.step == 'sequenza'
    } else {
        // If no step parameter is provided, run all steps (full workflow)
        run_preprocessing = true
        run_mutect2 = true
        run_facets = true
        run_sequenza = true
    }

    // Print workflow step settings
    log.info "================== Workflow Steps =================="
    log.info "Selected workflow step      = ${params.step ?: 'FULL WORKFLOW'}"
    log.info "Run preprocessing           = ${run_preprocessing}"
    log.info "Run Mutect2 variant calling = ${run_mutect2}"
    log.info "Run CNV analysis (FACETS)   = ${run_facets}"
    log.info "Run CNV analysis (Sequenza) = ${run_sequenza}"

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

    fastq               = samples.fastq
    bam_tumour_normal   = samples.bam_tumour_normal
    bam_tumour_only     = samples.bam_tumour_only

    /*
        ======================================================================
                        Run the full workflow
        ======================================================================
    */
    // Initialize channels for connecting workflow steps
    bam_for_mutect2_tn = Channel.empty()
    bam_for_mutect2_to = Channel.empty()
    
    // Store the preprocessing output to avoid "Parameter was not used" warning
    preprocessed_bam = Channel.empty()
    preprocessed_bai = Channel.empty()

    // Preprocessing fastq files and bwa mapping
    if (run_preprocessing) {
        log.info "Running preprocessing step..."
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
        
        // Store preprocessing output
        preprocessed_bam = PREPROCESSING.out.bam
        preprocessed_bai = PREPROCESSING.out.bai
        
        // Connect preprocessing to later steps in full workflow mode
        if (!params.step) {
            // Process paired tumor-normal samples
            bam_for_mutect2_tn = preprocessed_bam
                .join(preprocessed_bai)
                .map { meta, bam, bai -> [meta.patient_id, [meta, bam, bai]] }
                .groupTuple(by: 0)
                .map { patient_id, sample_entries ->
                    def tumor = sample_entries.find { it[0].status == 1 }
                    def normal = sample_entries.find { it[0].status == 0 }
                    
                    if (tumor && normal) {
                        def meta = [
                            id: "${tumor[0].sample_id}_vs_${normal[0].sample_id}",
                            patient_id: patient_id,
                            tumour_id: tumor[0].sample_id,
                            normal_id: normal[0].sample_id,
                            is_paired: true
                        ]
                        
                        return [meta, tumor[1], tumor[2], normal[1], normal[2]]
                    } else {
                        return null
                    }
                }
                .filter { it != null }

            // Process tumor-only samples
            bam_for_mutect2_to = preprocessed_bam
                .join(preprocessed_bai)
                .filter { meta, _bam, _bai -> meta.status == 1 }
                .map { meta, bam, bai ->
                    def tumors_with_normal = bam_for_mutect2_tn
                        .map { tn_meta, _t_bam, _t_bai, _n_bam, _n_bai -> tn_meta.tumour_id }
                        .toList()
                        .map { ids -> ids.contains(meta.sample_id) }
                        .ifEmpty { false }
                        
                    if (!tumors_with_normal.val) {
                        def meta_new = [
                            id: meta.sample_id,
                            patient_id: meta.patient_id,
                            tumour_id: meta.sample_id,
                            normal_id: null,
                            is_paired: false
                        ]
                        return [meta_new, bam, bai, null, null]
                    } else {
                        return null
                    }
                }
                .filter { it != null }
        }
    }

    // Run Mutect2 SNV/indels variant calling and annotation
    if (run_mutect2) {
        log.info "Running Mutect2 variant calling..."
        
        // Use different inputs depending on whether we're running the full workflow or just Mutect2
        def mutect2_tn = params.step ? bam_tumour_normal : bam_for_mutect2_tn.mix(bam_tumour_normal).unique { it[0].id }
        def mutect2_to = params.step ? bam_tumour_only : bam_for_mutect2_to.mix(bam_tumour_only).unique { it[0].id }
        
        MUTECT2_CALL(
            mutect2_tn,
            mutect2_to,
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
            blacklist,
            funcotator_resources,
            funcotator_ref_version,
            annovar_db,
            annovar_buildver,
            annovar_protocol,
            annovar_operation,
            annovar_xreffile
        )
    }

    // Run CNV analysis using FACETS
    if (run_facets) {
        log.info "Running CNV analysis using FACETS..."
        
        // Use preprocessed BAMs when in full workflow mode
        def facets_tn = params.step ? bam_tumour_normal : bam_for_mutect2_tn.mix(bam_tumour_normal).unique { it[0].id }
        def facets_to = params.step ? bam_tumour_only : bam_for_mutect2_to.mix(bam_tumour_only).unique { it[0].id }
        
        CNV_FACETS(
            facets_tn,
            facets_to,
            dbsnp,
            dbsnp_tbi,
            defined_normal,
            defined_normal_index
        )
    }

    // Run CNV analysis using Sequenza
    if (run_sequenza) {
        log.info "Running CNV analysis using Sequenza..."
        
        // Use preprocessed BAMs when in full workflow mode, maintain consistent naming
        def sequenza_tn = params.step ? bam_tumour_normal : bam_for_mutect2_tn.mix(bam_tumour_normal).unique { it[0].id }
        
        CNV_SEQUENZA(
            sequenza_tn,
            fasta,
            wigfile,
            window_size
        )
    }

}

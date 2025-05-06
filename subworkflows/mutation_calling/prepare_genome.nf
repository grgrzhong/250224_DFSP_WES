/*
 * Prepare reference genome and related resources for variant calling
 */

workflow PREPARE_GENOME {

    take:
        genome_id  // The genome identifier (e.g., 'GRCh38', 'GRCh19', etc.)

    main:
        // First, check if the genome_id is provided and exists in params.genomes
        if (genome_id == null) {
            error "No genome ID provided to PREPARE_GENOME workflow"
        }
        
        if (!params.genomes.containsKey(genome_id)) {
            error "Genome ID '${genome_id}' not found in configuration. Available genomes: ${params.genomes.keySet().join(', ')}"
        }

        // Setup reference channels - map to variables in genome.config
        fasta                = Channel.fromPath(params.genomes[genome_id].fasta, checkIfExists: true)
        fai                  = Channel.fromPath(params.genomes[genome_id].fai, checkIfExists: true)
        dict                 = Channel.fromPath(params.genomes[genome_id].dict, checkIfExists: true)
        
        // Map dbSNP and tbi
        dbsnp               = Channel.fromPath(params.genomes[genome_id].dbsnp, checkIfExists: true)
        dbsnp_tbi           = Channel.fromPath(params.genomes[genome_id].dbsnp_tbi, checkIfExists: true)
        
        // Map germline_resource and tbi
        germline_resource    = Channel.fromPath(params.genomes[genome_id].germline_resource, checkIfExists: true)
        germline_resource_tbi= Channel.fromPath(params.genomes[genome_id].germline_resource_tbi, checkIfExists: true)
        
        // Map panel_of_normals (pon) and tbi
        panel_of_normals     = params.genomes[genome_id].pon ? 
                               Channel.fromPath(params.genomes[genome_id].pon, checkIfExists: true) : 
                               Channel.empty()
        panel_of_normals_tbi = params.genomes[genome_id].pon_tbi ? 
                               Channel.fromPath(params.genomes[genome_id].pon_tbi, checkIfExists: true) : 
                               Channel.empty()
        
        // Map pileup variants (contamination_variants)
        pileup_variants      = Channel.fromPath(params.genomes[genome_id].contamination_variants, checkIfExists: true)
        pileup_variants_tbi  = Channel.fromPath(params.genomes[genome_id].contamination_variants_tbi, checkIfExists: true)
        
        // Map intervals
        intervals            = Channel.fromPath(params.genomes[genome_id].intervals, checkIfExists: true)
        bait_intervals       = Channel.fromPath(params.genomes[genome_id].bait_intervals, checkIfExists: true)
        target_intervals     = Channel.fromPath(params.genomes[genome_id].target_intervals, checkIfExists: true)

        targets            = Channel.fromPath(params.genomes[genome_id].targets, checkIfExists: true)
        // Additional resources
        funcotator_resources = params.genomes[genome_id].funcotator ? 
                               Channel.fromPath(params.genomes[genome_id].funcotator, checkIfExists: true) :
                               Channel.empty()
        annovar_db           = params.genomes[genome_id].annovar_db ? 
                               Channel.fromPath(params.genomes[genome_id].annovar_db, checkIfExists: true) :
                               Channel.empty()

         // Print genome resource information similar to nf-core/sarek
        log.info"""
        
        Fasta                : ${params.genomes[genome_id].fasta}
        Fasta index          : ${params.genomes[genome_id].fai}
        Dict                 : ${params.genomes[genome_id].dict}
        Germline resource    : ${params.genomes[genome_id].germline_resource}
        Panel of normals     : ${params.genomes[genome_id].containsKey('pon') ? params.genomes[genome_id].pon : 'Not provided'}
        Pileup variants      : ${params.genomes[genome_id].contamination_variants}
        Intervals            : ${params.genomes[genome_id].intervals}
        Targets              : ${params.genomes[genome_id].targets}
        Bait intervals       : ${params.genomes[genome_id].bait_intervals}
        Target intervals     : ${params.genomes[genome_id].target_intervals}
        dbSNP                : ${params.genomes[genome_id].dbsnp}
        dbSNP index          : ${params.genomes[genome_id].dbsnp_tbi}
        Funcotator resources : ${params.genomes[genome_id].containsKey('funcotator') ? params.genomes[genome_id].funcotator : 'Not provided'}
        Annovar database     : ${params.genomes[genome_id].containsKey('annovar_db') ? params.genomes[genome_id].annovar_db : 'Not provided'}
        
        """

    emit:
        fasta                = fasta
        fai                  = fai
        dict                 = dict
        dbsnp                = dbsnp
        dbsnp_tbi            = dbsnp_tbi
        germline_resource    = germline_resource
        germline_resource_tbi= germline_resource_tbi
        panel_of_normals     = panel_of_normals
        panel_of_normals_tbi = panel_of_normals_tbi
        pileup_variants      = pileup_variants
        pileup_variants_tbi  = pileup_variants_tbi
        intervals            = intervals
        bait_intervals       = bait_intervals
        target_intervals     = target_intervals
        targets              = targets
        funcotator_resources = funcotator_resources
        annovar_db           = annovar_db
}
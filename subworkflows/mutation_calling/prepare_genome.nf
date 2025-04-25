
/*
========================================================================================
This subworkflow prepares the reference genome and its associated resources.
========================================================================================
*/

// Load required modules
include { SAMTOOLS_FAIDX                 } from '../../modules/variant_calling/samtools/faidx'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/variant_calling/gatk4/createsequencedictionary'
include { TABIX_TABIX                    } from '../../modules/variant_calling/tabix/tabix'

workflow PREPARE_GENOME {
    take:
    genome_name        // val: genome name (e.g. 'GRCh38')

    main:
    ch_versions = Channel.empty()
    log.info "Preparing reference genome: ${genome_name}"
    
    // Check if genome exists in config
    if (!params.genomes.containsKey(genome_name)) {
        error "Genome ${genome_name} not found in params.genomes configuration"
    }
    
    // Reference FASTA - required
    def fasta_file = file(params.genomes[genome_name].fasta, checkIfExists: true)
    def fasta_ch = Channel.value(fasta_file)
    
    // FASTA index - generate if missing
    def fai_ch
    if (params.genomes[genome_name].containsKey('fai')) {
        def fai_file = file(params.genomes[genome_name].fai, checkIfExists: false)
        if (fai_file.exists()) {
            log.info "Using existing FASTA index: ${fai_file}"
            fai_ch = Channel.value(fai_file)
        } else {
            log.info "FASTA index file not found: ${fai_file}. Will be generated."
            fai_ch = SAMTOOLS_FAIDX(fasta_ch).fai
            ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
        }
    } else {
        log.info "FASTA index not specified. Will be generated."
        fai_ch = SAMTOOLS_FAIDX(fasta_ch).fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }
    
    // Sequence dictionary - generate if missing
    def dict_ch
    if (params.genomes[genome_name].containsKey('dict')) {
        def dict_file = file(params.genomes[genome_name].dict, checkIfExists: false)
        if (dict_file.exists()) {
            log.info "Using existing sequence dictionary: ${dict_file}"
            dict_ch = Channel.value(dict_file)
        } else {
            log.info "Sequence dictionary file not found: ${dict_file}. Will be generated."
            dict_ch = GATK4_CREATESEQUENCEDICTIONARY(fasta_ch).dict
            ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
        }
    } else {
        log.info "Sequence dictionary not specified. Will be generated."
        dict_ch = GATK4_CREATESEQUENCEDICTIONARY(fasta_ch).dict
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }
    
    // dbSNP VCF and index
    def dbsnp_ch = Channel.empty()
    def dbsnp_tbi_ch = Channel.empty()
    if (params.genomes[genome_name].containsKey('dbsnp')) {
        def dbsnp_file = file(params.genomes[genome_name].dbsnp, checkIfExists: false)
        if (dbsnp_file.exists()) {
            dbsnp_ch = Channel.value(dbsnp_file)
            
            // Check/create index for dbSNP
            if (params.genomes[genome_name].containsKey('dbsnp_tbi')) {
                def dbsnp_tbi_file = file(params.genomes[genome_name].dbsnp_tbi, checkIfExists: false)
                if (dbsnp_tbi_file.exists()) {
                    dbsnp_tbi_ch = Channel.value(dbsnp_tbi_file)
                } else {
                    log.info "dbSNP index not found: ${dbsnp_tbi_file}. Will be generated."
                    dbsnp_tbi_ch = TABIX_TABIX(dbsnp_ch).tbi
                    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
                }
            } else {
                log.info "dbSNP index not specified. Will be generated."
                dbsnp_tbi_ch = TABIX_TABIX(dbsnp_ch).tbi
                ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
            }
        } else {
            log.warn "dbSNP file not found: ${dbsnp_file}. Will proceed without dbSNP annotations."
        }
    } else {
        log.info "dbSNP not specified in genome config"
    }
    
    // Intervals files (probes, targets, etc.)
    def intervals_ch = Channel.empty()
    def bait_intervals_ch = Channel.empty()
    def target_intervals_ch = Channel.empty()
    def exome_bed_ch = Channel.empty()
    
    if (params.genomes[genome_name].containsKey('intervals')) {
        def intervals_file = file(params.genomes[genome_name].intervals, checkIfExists: false)
        if (intervals_file.exists()) {
            intervals_ch = Channel.value(intervals_file)
        } else {
            log.warn "Intervals file not found: ${intervals_file}. Some processes may fail."
        }
    }
    
    if (params.genomes[genome_name].containsKey('bait_intervals')) {
        def bait_file = file(params.genomes[genome_name].bait_intervals, checkIfExists: false)
        if (bait_file.exists()) {
            bait_intervals_ch = Channel.value(bait_file)
        } else {
            log.warn "Bait intervals file not found: ${bait_file}. Some processes may fail."
        }
    }
    
    if (params.genomes[genome_name].containsKey('target_intervals')) {
        def target_file = file(params.genomes[genome_name].target_intervals, checkIfExists: false)
        if (target_file.exists()) {
            target_intervals_ch = Channel.value(target_file)
        } else {
            log.warn "Target intervals file not found: ${target_file}. Some processes may fail."
        }
    }
    
    if (params.genomes[genome_name].containsKey('exome_bed')) {
        def bed_file = file(params.genomes[genome_name].exome_bed, checkIfExists: false)
        if (bed_file.exists()) {
            exome_bed_ch = Channel.value(bed_file)
        } else {
            log.warn "Exome BED file not found: ${bed_file}. Some processes may fail."
        }
    }
    
    // Known indels with index checking/creation
    def known_indels_ch = Channel.empty()
    def known_indels_tbi_ch = Channel.empty()
    
    if (params.genomes[genome_name].containsKey('known_indels')) {
        def known_indels_file = file(params.genomes[genome_name].known_indels, checkIfExists: false)
        if (known_indels_file.exists()) {
            known_indels_ch = Channel.value(known_indels_file)
            
            // Check/create index for known indels
            if (params.genomes[genome_name].containsKey('known_indels_tbi')) {
                def known_indels_tbi_file = file(params.genomes[genome_name].known_indels_tbi, checkIfExists: false)
                if (known_indels_tbi_file.exists()) {
                    known_indels_tbi_ch = Channel.value(known_indels_tbi_file)
                } else {
                    log.info "Known indels index not found: ${known_indels_tbi_file}. Will be generated."
                    known_indels_tbi_ch = TABIX_TABIX(known_indels_ch).tbi
                    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
                }
            } else {
                log.info "Known indels index not specified. Will be generated."
                known_indels_tbi_ch = TABIX_TABIX(known_indels_ch).tbi
                ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
            }
        } else {
            log.warn "Known indels file not found: ${known_indels_file}. BQSR may be less accurate."
        }
    }
    
    // Germline resource with index checking/creation
    def germline_resource_ch = Channel.empty()
    def germline_resource_tbi_ch = Channel.empty()
    
    if (params.genomes[genome_name].containsKey('germline_resource')) {
        def germline_file = file(params.genomes[genome_name].germline_resource, checkIfExists: false)
        if (germline_file.exists()) {
            germline_resource_ch = Channel.value(germline_file)
            
            // Check/create index for germline resource
            if (params.genomes[genome_name].containsKey('germline_resource_tbi')) {
                def germline_tbi_file = file(params.genomes[genome_name].germline_resource_tbi, checkIfExists: false)
                if (germline_tbi_file.exists()) {
                    germline_resource_tbi_ch = Channel.value(germline_tbi_file)
                } else {
                    log.info "Germline resource index not found: ${germline_tbi_file}. Will be generated."
                    germline_resource_tbi_ch = TABIX_TABIX(germline_resource_ch).tbi
                    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
                }
            } else {
                log.info "Germline resource index not specified. Will be generated."
                germline_resource_tbi_ch = TABIX_TABIX(germline_resource_ch).tbi
                ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
            }
        } else {
            log.warn "Germline resource file not found: ${germline_file}. Somatic variant calling may be affected."
        }
    }
    
    // Panel of normals with index checking/creation
    def pon_ch = Channel.empty()
    def pon_tbi_ch = Channel.empty()
    
    if (params.genomes[genome_name].containsKey('pon')) {
        def pon_file = file(params.genomes[genome_name].pon, checkIfExists: false)
        if (pon_file.exists()) {
            pon_ch = Channel.value(pon_file)
            
            // Check/create index for panel of normals
            if (params.genomes[genome_name].containsKey('pon_tbi')) {
                def pon_tbi_file = file(params.genomes[genome_name].pon_tbi, checkIfExists: false)
                if (pon_tbi_file.exists()) {
                    pon_tbi_ch = Channel.value(pon_tbi_file)
                } else {
                    log.info "Panel of normals index not found: ${pon_tbi_file}. Will be generated."
                    pon_tbi_ch = TABIX_TABIX(pon_ch).tbi
                    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
                }
            } else {
                log.info "Panel of normals index not specified. Will be generated."
                pon_tbi_ch = TABIX_TABIX(pon_ch).tbi
                ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
            }
        } else {
            log.warn "Panel of normals file not found: ${pon_file}. Somatic variant calling may have more false positives."
        }
    }
    
    // Variants for contamination
    def variants_for_contamination_ch = Channel.empty()
    def variants_for_contamination_tbi_ch = Channel.empty()
    
    if (params.genomes[genome_name].containsKey('variants_for_contamination')) {
        def contam_file = file(params.genomes[genome_name].variants_for_contamination, checkIfExists: false)
        if (contam_file.exists()) {
            variants_for_contamination_ch = Channel.value(contam_file)
            
            // Check/create index for contamination variants
            if (params.genomes[genome_name].containsKey('variants_for_contamination_tbi')) {
                def contam_tbi_file = file(params.genomes[genome_name].variants_for_contamination_tbi, checkIfExists: false)
                if (contam_tbi_file.exists()) {
                    variants_for_contamination_tbi_ch = Channel.value(contam_tbi_file)
                } else {
                    log.info "Contamination variants index not found: ${contam_tbi_file}. Will be generated."
                    variants_for_contamination_tbi_ch = TABIX_TABIX(variants_for_contamination_ch).tbi
                    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
                }
            } else {
                log.info "Contamination variants index not specified. Will be generated."
                variants_for_contamination_tbi_ch = TABIX_TABIX(variants_for_contamination_ch).tbi
                ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
            }
        } else {
            log.warn "Contamination variants file not found: ${contam_file}. Contamination estimation will be skipped."
        }
    }
    
    // Funcotator database
    def funcotator_db_ch = Channel.empty()
    if (params.genomes[genome_name].containsKey('funcotator')) {
        def funcotator_db = file(params.genomes[genome_name].funcotator, checkIfExists: false)
        if (funcotator_db.exists() && funcotator_db.isDirectory()) {
            funcotator_db_ch = Channel.value(funcotator_db)
        } else {
            log.warn "Funcotator database not found: ${funcotator_db}. Funcotator annotation will be skipped."
        }
    }
    
    // Annovar database
    def annovar_db_ch = Channel.empty()
    if (params.genomes[genome_name].containsKey('annovar_db')) {
        def annovar_db = file(params.genomes[genome_name].annovar_db, checkIfExists: false)
        if (annovar_db.exists() && annovar_db.isDirectory()) {
            annovar_db_ch = Channel.value(annovar_db)
        } else {
            log.warn "Annovar database not found: ${annovar_db}. Annovar annotation will be skipped."
        }
    }
    
    // Prepare Annovar parameters
    def annovar_buildver_ch = params.genomes[genome_name].containsKey('annovar_buildver') ? 
                            Channel.value(params.genomes[genome_name].annovar_buildver) : 
                            Channel.value('hg38')
                            
    def annovar_protocols_ch = params.genomes[genome_name].containsKey('annovar_protocols') ? 
                            Channel.value(params.genomes[genome_name].annovar_protocols.join(',')) : 
                            Channel.empty()
                            
    def annovar_operation_ch = params.genomes[genome_name].containsKey('annovar_operation') ? 
                            Channel.value(params.genomes[genome_name].annovar_operation.join(',')) : 
                            Channel.empty()

    // Summary log of available resources
    log.info "-" * 60
    log.info "Genome Reference Resources Summary"
    log.info "-" * 60
    log.info "Genome name:           ${genome_name}"
    log.info "FASTA:                 ${params.genomes[genome_name].fasta}"
    log.info "Intervals:             ${params.genomes[genome_name].containsKey('intervals') ? params.genomes[genome_name].intervals : 'Not provided'}"
    log.info "dbSNP:                 ${params.genomes[genome_name].containsKey('dbsnp') ? params.genomes[genome_name].dbsnp : 'Not provided'}"
    log.info "Known indels:          ${params.genomes[genome_name].containsKey('known_indels') ? params.genomes[genome_name].known_indels : 'Not provided'}"
    log.info "Germline resource:     ${params.genomes[genome_name].containsKey('germline_resource') ? params.genomes[genome_name].germline_resource : 'Not provided'}"
    log.info "Panel of normals:      ${params.genomes[genome_name].containsKey('pon') ? params.genomes[genome_name].pon : 'Not provided'}"
    log.info "Contamination vars:    ${params.genomes[genome_name].containsKey('variants_for_contamination') ? params.genomes[genome_name].variants_for_contamination : 'Not provided'}"
    log.info "-" * 60
    
    emit:
    fasta                          = fasta_ch
    fai                            = fai_ch
    dict                           = dict_ch
    dbsnp                          = dbsnp_ch
    dbsnp_tbi                      = dbsnp_tbi_ch
    intervals                      = intervals_ch
    bait_intervals                 = bait_intervals_ch
    target_intervals               = target_intervals_ch
    exome_bed                      = exome_bed_ch
    known_indels                   = known_indels_ch
    known_indels_tbi               = known_indels_tbi_ch
    germline_resource              = germline_resource_ch
    germline_resource_tbi          = germline_resource_tbi_ch
    pon                            = pon_ch
    pon_tbi                        = pon_tbi_ch
    variants_for_contamination     = variants_for_contamination_ch
    variants_for_contamination_tbi = variants_for_contamination_tbi_ch
    funcotator_db                  = funcotator_db_ch
    annovar_db                     = annovar_db_ch
    annovar_buildver               = annovar_buildver_ch
    annovar_protocols              = annovar_protocols_ch
    annovar_operation              = annovar_operation_ch
    versions                       = ch_versions
}
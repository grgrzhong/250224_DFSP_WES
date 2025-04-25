/*
 * Mutect2 somatic variant calling workflow
 * Handles both paired tumor-normal and tumor-only samples
 */
include { GATK4_GETPILEUPSUMMARIES                                      } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as  GATK4_GETPILEUPSUMMARIES_NORMAL  } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_CALCULATECONTAMINATION                                  } from '../../modules/variant_calling/gatk4/calculatecontamination'
include { GATK4_MUTECT2 as GATK4_MUTECT2_TUMOR_NORMAL                   } from '../../modules/variant_calling/gatk4/mutect2'
include { GATK4_MUTECT2 as GATK4_MUTECT2_TUMOR_ONLY                     } from '../../modules/variant_calling/gatk4/mutect2'
include { GATK4_LEARNREADORIENTATIONMODEL                               } from '../../modules/variant_calling/gatk4/learnreadorientationmodel'
include { GATK4_FILTERMUTECTCALLS                                       } from '../../modules/variant_calling/gatk4/filtermutectcalls'
include { BCFTOOLS_NORM                                                 } from '../../modules/variant_calling/bcftools/norm'
include { BCFTOOLS_VIEW                                                 } from '../../modules/variant_calling/bcftools/view'
include { GATK4_FUNCOTATOR                                              } from '../../modules/variant_calling/gatk4/funcotator'
include { ANNOVAR                                                       } from '../../modules/variant_calling/annovar'
include { TABIX_TABIX                                                   } from '../../modules/utilities/tabix/tabix'

// input
params.input            = "/home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test.csv"
params.outdir           = "${launchDir}/results"
params.publish_dir_mode = "copy"

workflow MUTECT2_CALL {
    take:
    input_samples        // Channel: [ val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai) ]
    fasta                // Channel: path(reference_fasta)
    fai                  // Channel: path(reference_fasta_index)
    dict                 // Channel: path(reference_dict)
    germline_resource    // Channel: path(gnomad_vcf)
    germline_resource_tbi// Channel: path(gnomad_vcf_index)
    panel_of_normals     // Channel: path(pon_vcf)
    panel_of_normals_tbi // Channel: path(pon_vcf_index)
    intervals            // Channel: path(intervals_bed)
    pileup_variants      // Channel: path(small_exac_common)
    pileup_variants_tbi  // Channel: path(small_exac_common_index)
    funcotator_resources // Channel: path(funcotator_data_sources)
    annovar_db           // Channel: path(annovar_humandb)
    gene_xref            // Channel: path(gene_fullxref.txt)

    main:
    ch_versions = Channel.empty()
    
    // Branch samples into paired and unpaired
    input_samples
        .branch {
            paired: it[3] != null  // Check if normal_bam exists
            unpaired: true
        }
        .set { sample_branches }
    
    // 1. Get Pileup Summaries for tumor samples
    GATK4_GETPILEUPSUMMARIES(
        input_samples.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            // Extract tumor sample for pileup
            [meta, tumor_bam, tumor_bai]
        },
        pileup_variants,
        pileup_variants_tbi,
        intervals
    )
    ch_versions = ch_versions.mix(GATK4_GETPILEUPSUMMARIES.out.versions)
    
    // 2a. Get Pileup Summaries for normal samples (paired only)
    GATK4_GETPILEUPSUMMARIES_NORMAL(
        sample_branches.paired.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            // Extract normal sample for pileup
            [meta, normal_bam, normal_bai]
        },
        pileup_variants,
        pileup_variants_tbi,
        intervals
    )
    ch_versions = ch_versions.mix(GATK4_GETPILEUPSUMMARIES_NORMAL.out.versions.first().ifEmpty([]))
    
    // Join tumor pileups with samples
    ch_tumor_pileups = input_samples.join(GATK4_GETPILEUPSUMMARIES.out.table, by: 0)
    
    // 2b. Calculate Contamination (different for paired vs unpaired)
    // For paired samples
    ch_paired_pileups = ch_tumor_pileups
        .join(GATK4_GETPILEUPSUMMARIES_NORMAL.out.table, by: 0)
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, tumor_table, normal_table ->
            [meta, tumor_table, normal_table]
        }
    
    GATK4_CALCULATECONTAMINATION(
        ch_paired_pileups.map { meta, tumor_table, normal_table ->
            [meta, tumor_table]
        },
        ch_paired_pileups.map { meta, tumor_table, normal_table ->
            [meta, normal_table]
        }
    )
    
    // For unpaired samples
    ch_unpaired_pileups = sample_branches.unpaired
        .join(GATK4_GETPILEUPSUMMARIES.out.table, by: 0)
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, tumor_table ->
            [meta, tumor_table]
        }
    
    GATK4_CALCULATECONTAMINATION.run(
        ch_unpaired_pileups,
        Channel.of() // Empty channel for matched pileup
    )
    
    // Combine all contamination results
    ch_contamination_table = GATK4_CALCULATECONTAMINATION.out.contamination
    ch_segments_table = GATK4_CALCULATECONTAMINATION.out.segments
    ch_versions = ch_versions.mix(GATK4_CALCULATECONTAMINATION.out.versions.first().ifEmpty([]))
    
    // 3a. Run Mutect2 for paired samples
    GATK4_MUTECT2_TUMOR_NORMAL(
        sample_branches.paired,
        fasta,
        fai,
        dict,
        panel_of_normals,
        panel_of_normals_tbi,
        germline_resource,
        germline_resource_tbi,
        intervals
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2_TUMOR_NORMAL.out.versions.first().ifEmpty([]))
    
    // 3b. Run Mutect2 for unpaired samples
    GATK4_MUTECT2_TUMOR_ONLY(
        sample_branches.unpaired,
        fasta,
        fai,
        dict,
        panel_of_normals,
        panel_of_normals_tbi,
        germline_resource,
        germline_resource_tbi,
        intervals
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2_TUMOR_ONLY.out.versions.first().ifEmpty([]))
    
    // Merge Mutect2 results from paired and unpaired samples
    ch_vcf = GATK4_MUTECT2_TUMOR_NORMAL.out.vcf.mix(GATK4_MUTECT2_TUMOR_ONLY.out.vcf)
    ch_tbi = GATK4_MUTECT2_TUMOR_NORMAL.out.tbi.mix(GATK4_MUTECT2_TUMOR_ONLY.out.tbi)
    ch_stats = GATK4_MUTECT2_TUMOR_NORMAL.out.stats.mix(GATK4_MUTECT2_TUMOR_ONLY.out.stats)
    ch_f1r2 = GATK4_MUTECT2_TUMOR_NORMAL.out.f1r2.mix(GATK4_MUTECT2_TUMOR_ONLY.out.f1r2)
    
    // 4. Learn Read Orientation Model
    GATK4_LEARNREADORIENTATIONMODEL(ch_f1r2)
    ch_versions = ch_versions.mix(GATK4_LEARNREADORIENTATIONMODEL.out.versions.first())
    
    // 5. Filter Mutect Calls
    // Join necessary channels for filtering
    ch_filter_input = ch_vcf
        .join(ch_tbi, by: 0)
        .join(ch_stats, by: 0)
        .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifact_prior, by: 0)
        .join(ch_contamination_table, by: 0)
        .join(ch_segments_table, by: 0)
    
    GATK4_FILTERMUTECTCALLS(
        ch_filter_input.map { meta, vcf, tbi, stats, artifact_prior, contamination, segments ->
            [meta, vcf, tbi, stats]
        },
        ch_filter_input.map { meta, vcf, tbi, stats, artifact_prior, contamination, segments ->
            artifact_prior
        },
        ch_filter_input.map { meta, vcf, tbi, stats, artifact_prior, contamination, segments ->
            contamination
        },
        ch_filter_input.map { meta, vcf, tbi, stats, artifact_prior, contamination, segments ->
            segments
        },
        fasta,
        fai,
        dict
    )
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions.first())
    
    // 6. Normalize variants with BCFtools
    BCFTOOLS_NORM(
        GATK4_FILTERMUTECTCALLS.out.vcf,
        fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())
    
    // 7. Filter for PASS variants
    BCFTOOLS_VIEW(
        BCFTOOLS_NORM.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())
    
    // Create indices for filtered VCFs
    TABIX_TABIX(BCFTOOLS_VIEW.out.vcf)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())
    
    // Combine VCF with its index for downstream processes
    ch_filtered_vcf_with_tbi = BCFTOOLS_VIEW.out.vcf.join(TABIX_TABIX.out.tbi, by: 0)
    
    // 8. Annotate with Funcotator
    GATK4_FUNCOTATOR(
        ch_filtered_vcf_with_tbi,
        fasta,
        fai,
        dict,
        funcotator_resources,
        intervals
    )
    ch_versions = ch_versions.mix(GATK4_FUNCOTATOR.out.versions.first())
    
    // 9. Annotate with Annovar
    ANNOVAR(
        BCFTOOLS_VIEW.out.vcf,
        annovar_db,
        gene_xref
    )
    ch_versions = ch_versions.mix(ANNOVAR.out.versions.first())
    
    emit:
    vcf_filtered          = BCFTOOLS_VIEW.out.vcf          // Channel: [ val(meta), path(vcf) ]
    vcf_tbi               = TABIX_TABIX.out.tbi            // Channel: [ val(meta), path(tbi) ]
    funcotator_maf        = GATK4_FUNCOTATOR.out.maf       // Channel: [ val(meta), path(maf) ]
    annovar_txt           = ANNOVAR.out.txt                // Channel: [ val(meta), path(txt) ]
    versions              = ch_versions                     // Channel: [ path(versions.yml) ]
}

workflow {
    // Check input
    if (params.input == null) {
        error "Please provide an input CSV file with --input"
    }
    
    // Setup reference channels
    fasta                = Channel.fromPath(params.genomes[params.genome].fasta, checkIfExists: true)
    fai                  = Channel.fromPath(params.genomes[params.genome].fai, checkIfExists: true)
    dict                 = Channel.fromPath(params.genomes[params.genome].dict, checkIfExists: true)
    germline_resource    = Channel.fromPath(params.genomes[params.genome].germline_resource, checkIfExists: true)
    germline_resource_tbi= Channel.fromPath(params.genomes[params.genome].germline_resource_tbi, checkIfExists: true)
    panel_of_normals     = params.genomes[params.genome].panel_of_normals ? 
                           Channel.fromPath(params.genomes[params.genome].panel_of_normals, checkIfExists: true) : 
                           Channel.empty()
    panel_of_normals_tbi = params.genomes[params.genome].panel_of_normals_tbi ? 
                           Channel.fromPath(params.genomes[params.genome].panel_of_normals_tbi, checkIfExists: true) : 
                           Channel.empty()
    pileup_variants      = Channel.fromPath(params.genomes[params.genome].pileup_variants, checkIfExists: true)
    pileup_variants_tbi  = Channel.fromPath(params.genomes[params.genome].pileup_variants_tbi, checkIfExists: true)
    intervals            = Channel.fromPath(params.genomes[params.genome].intervals, checkIfExists: true)
    funcotator_resources = Channel.fromPath(params.funcotator_resources, checkIfExists: true)
    annovar_db           = Channel.fromPath(params.annovar_db, checkIfExists: true)
    gene_xref            = Channel.fromPath(params.gene_xref, checkIfExists: true)

    // Parse input csv
    input_samples = Channel.fromPath(params.input)
        .ifEmpty { exit(1, "Samplesheet not found: ${params.input}") }
        .splitCsv(header: true)
        .map { row ->
            // Extract sample info
            def patient_id = row.patient ? row.patient.trim() : null
            def sample_id = row.sample ? row.sample.trim() : null
            def status = row.status ? row.status.trim() : null

            // Validate required fields
            if (!patient_id) {
                error("Missing or empty 'patient' field in row: ${row}")
            }
            if (!sample_id) {
                error("Missing or empty 'sample' field in row: ${row}")
            }
            if (!status) {
                error("Missing or empty 'status' field in row: ${row}")
            }
            if (!row.bam) {
                error("Missing 'bam' field in row: ${row}")
            }

            // Process BAM and BAI paths
            def bam = file(row.bam)
            def bai = row.bai ? file(row.bai) : file("${row.bam}.bai")

            // Check if files exist
            if (!bam.exists()) {
                error("BAM file not found: ${bam}")
            }
            if (!bai.exists()) {
                error("BAI file not found: ${bai}")
            }

            // Return a tuple with patient_id, sample_id, status, bam, and bai
            [
                patient_id: patient_id,
                sample_id: sample_id,
                status: status.toInteger(),
                bam: bam,
                bai: bai,
            ]
        }

    // Split samples into tumour and normal
    tumour_samples = input_samples.filter { it.status == 1 }
    normal_samples = input_samples.filter { it.status == 0 }

    // Create paired samples channel
    paired_samples = tumour_samples
        .map { tumour  -> [ tumour.patient_id, tumour] }
        .combine(
            normal_samples.map { normal -> [normal.patient_id, normal] },
            by: 0
        )
        .map {
            patient_id, tumour, normal -> 
            def meta = [
                id: "${tumour.sample_id}_vs_${normal.sample_id}",
                patient_id: patient_id,
                tumour_id: tumour.sample_id,
                normal_id: normal.sample_id,
                is_paired: true
            ]
            [meta, tumour.bam, tumour.bai, normal.bam, normal.bai]
        }
    
    // Create tumour-only samples channel
    normal_patient_ids = normal_samples
        .map { it.patient_id }
        .unique()
        .collect()
        .map { ids  -> ids.sort() }

    unpaired_samples = tumour_samples
        .branch {
            def normal_ids = normal_patient_ids.val
            paired: normal_ids.contains(it.patient_id)
            unpaired: true
        }
        .unpaired
        .map {
            tumour -> def meta = [
                id: tumour.sample_id,
                patient_id: tumour.patient_id,
                tumour_id: tumour.sample_id,
                normal_id: null,
                is_paired: false
            ]
            [meta, tumour.bam, tumour.bai, null, null]
        }
    
    // Combine paired and unpaired samples
    all_samples = paired_samples.mix(unpaired_samples)
    
    // Log sample info
    paired_samples.count().subscribe { count ->
        log.info("Found ${count} paired samples.")
    }
    
    unpaired_samples.count().subscribe { count ->
        log.info("Found ${count} unpaired samples.")
    }
    
    // Execute the mutation calling workflow
    MUTECT2_CALL(
        all_samples,
        fasta,
        fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi,
        intervals,
        pileup_variants,
        pileup_variants_tbi,
        funcotator_resources,
        annovar_db,
        gene_xref
    )
}
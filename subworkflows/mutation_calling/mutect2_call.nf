/*

    Mutect2 calling subworkflow
*/
include { PREPARE_GENOME                                                             } from '../../subworkflows/mutation_calling/prepare_genome.nf'
include { GATK4_GETPILEUPSUMMARIES as  GATK4_GETPILEUPSUMMARIES_TUMOUR_NORMAL        } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as  GATK4_GETPILEUPSUMMARIES_TUMOUR_ONLY          } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as  GATK4_GETPILEUPSUMMARIES_NORMAL               } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_CALCULATECONTAMINATION as GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL } from '../../modules/variant_calling/gatk4/calculatecontamination'
include { GATK4_CALCULATECONTAMINATION as GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY   } from '../../modules/variant_calling/gatk4/calculatecontamination'
include { GATK4_MUTECT2 as GATK4_MUTECT2_TUMOR_NORMAL                                } from '../../modules/variant_calling/gatk4/mutect2'
include { GATK4_MUTECT2 as GATK4_MUTECT2_TUMOR_ONLY                                  } from '../../modules/variant_calling/gatk4/mutect2'
include { GATK4_LEARNREADORIENTATIONMODEL                                            } from '../../modules/variant_calling/gatk4/learnreadorientationmodel'
include { GATK4_FILTERMUTECTCALLS                                                    } from '../../modules/variant_calling/gatk4/filtermutectcalls'
include { BCFTOOLS_NORM                                                              } from '../../modules/variant_calling/bcftools/norm'
include { BCFTOOLS_VIEW                                                              } from '../../modules/variant_calling/bcftools/view'
include { GATK4_FUNCOTATOR                                                           } from '../../modules/variant_calling/gatk4/funcotator'
include { ANNOVAR                                                                    } from '../../modules/variant_calling/annovar'


// input
params.input = "/home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test.csv"

workflow {
    // Check input
    if (params.input == null) {
        error "Please provide an input CSV file with --input"
    }
    
    // Prepare the genome and resource files
    PREPARE_GENOME(params.genome)
    
    // Extract reference channels from PREPARE_GENOME
    fasta                = PREPARE_GENOME.out.fasta
    fai                  = PREPARE_GENOME.out.fai
    dict                 = PREPARE_GENOME.out.dict
    germline_resource    = PREPARE_GENOME.out.germline_resource
    germline_resource_tbi= PREPARE_GENOME.out.germline_resource_tbi
    panel_of_normals     = PREPARE_GENOME.out.panel_of_normals
    panel_of_normals_tbi = PREPARE_GENOME.out.panel_of_normals_tbi
    pileup_variants      = PREPARE_GENOME.out.pileup_variants
    pileup_variants_tbi  = PREPARE_GENOME.out.pileup_variants_tbi
    intervals            = PREPARE_GENOME.out.intervals
    funcotator_resources = PREPARE_GENOME.out.funcotator_resources
    annovar_db           = PREPARE_GENOME.out.annovar_db

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
    
    all_samples.view()
    // Log sample info
    paired_samples.count().subscribe { count ->
        log.info("Found ${count} paired samples.")
    }
    
    unpaired_samples.count().subscribe { count ->
        log.info("Found ${count} unpaired samples.")
    }

    /*
        ===================== TUMOR-NORMAL PAIRED ANALYSIS ====================
    */
    if (paired_samples) {

        // Extract tumor samples from paired samples for pileup
        paired_tumor_samples = paired_samples
            .map { 
                meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                [meta, tumor_bam, tumor_bai]
        }

        // Get pileup summaries for tumor samples in paired mode
        GATK4_GETPILEUPSUMMARIES_TUMOUR_NORMAL(
            paired_tumor_samples,
            pileup_variants,
            pileup_variants_tbi,
            intervals
        )

        // Get pileup summaries for normal samples
        normal_pileups = paired_samples
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                def normal_meta = meta.clone()
                normal_meta.id = "${meta.normal_id}_pileup"
                [normal_meta, normal_bam, normal_bai]
        }

        GATK4_GETPILEUPSUMMARIES_NORMAL(
            normal_pileups,
            pileup_variants,
            pileup_variants_tbi,
            intervals
        )

        // Calculate contamination for paired samples
        paired_contamination_input = 
            GATK4_GETPILEUPSUMMARIES_TUMOUR_NORMAL.out.table
            .join(GATK4_GETPILEUPSUMMARIES_NORMAL.out.table, by: [0])
            .map { meta, tumor_table, normal_table ->
                [meta, tumor_table, normal_table]
            }

        GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL(paired_contamination_input)

        // Run Mutect2 for paired samples
        GATK4_MUTECT2_TUMOR_NORMAL(
            paired_samples,
            fasta,
            fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            intervals
        )

        // Learn read orientation models for paired samples
        GATK4_LEARNREADORIENTATIONMODEL(
            GATK4_MUTECT2_TUMOR_NORMAL.out.f1r2
        )

        // Filter Mutect2 calls for paired samples
        paired_filter_input = GATK4_MUTECT2_TUMOR_NORMAL.out.vcf
            .join(GATK4_MUTECT2_TUMOR_NORMAL.out.tbi, by: [0])
            .join(GATK4_MUTECT2_TUMOR_NORMAL.out.stats, by: [0])
            .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior, by: [0])
            .join(GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL.out.contamination, by: [0])
            .join(GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL.out.segmentation, by: [0])
            .map {
                meta, vcf, tbi, stats, artifactprior, contamination, segmentation ->
                [meta, vcf, tbi, stats, artifactprior, contamination, segmentation]
            }

        GATK4_FILTERMUTECTCALLS(
            paired_filter_input,
            fasta,
            fai,
            dict
        )

        // Store paired filtered VCFs for later processing
        paired_filtered_vcfs = GATK4_FILTERMUTECTCALLS.out.vcf
        paired_filtered_tbis = GATK4_FILTERMUTECTCALLS.out.tbi
    }

    /*
        =================== TUMOR-ONLY UNPAIRED ANALYSIS ====================
    */
    if (unpaired_samples) {
        // Extract tumor samples for unpaired analysis
        unpaired_tumor_samples = unpaired_samples
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            [meta, tumor_bam, tumor_bai]
        }

        // Get pileup summaries for tumor-only samples
        GATK4_GETPILEUPSUMMARIES_TUMOUR_ONLY(
            unpaired_tumor_samples,
            pileup_variants,
            pileup_variants_tbi,
            intervals
        )

        // Calculate contamination for tumor-only samples
        unpaired_contamination_input = GATK4_GETPILEUPSUMMARIES_TUMOUR_ONLY.out.table
            .map { meta, table ->
                [meta, table, []]  // Empty file for matched normal
            }

        GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY(unpaired_contamination_input)

        // Run Mutect2 for tumor-only samples
        GATK4_MUTECT2_TUMOR_ONLY(
            unpaired_tumor_samples,
            fasta,
            fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            intervals
        )

        // Learn read orientation models for tumor-only samples
        GATK4_LEARNREADORIENTATIONMODEL(
            GATK4_MUTECT2_TUMOR_ONLY.out.f1r2
        )

        // Filter Mutect2 calls for tumor-only samples
        unpaired_filter_input = GATK4_MUTECT2_TUMOR_ONLY.out.vcf
            .join(GATK4_MUTECT2_TUMOR_ONLY.out.tbi, by: [0])
            .join(GATK4_MUTECT2_TUMOR_ONLY.out.stats, by: [0])
            .join(GATK4_LEARNREADORIENTATIONMODEL.out.artifactprior, by: [0])
            .join(GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY.out.contamination, by: [0])
            .join(GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY.out.segmentation, by: [0])
            .map { meta, vcf, tbi, stats, artifactprior, contamination, segmentation ->
                [meta, vcf, tbi, stats, artifactprior, contamination, segmentation]
            }

        GATK4_FILTERMUTECTCALLS(
            unpaired_filter_input,
            fasta,
            fai,
            dict
        )

        // Store unpaired filtered VCFs for later processing
        unpaired_filtered_vcfs = GATK4_FILTERMUTECTCALLS.out.vcf
        unpaired_filtered_tbis = GATK4_FILTERMUTECTCALLS.out.tbi
    }

    /*
        ===================== COMBINED DOWNSTREAM PROCESSING ============================
    */
    // Combine all filtered VCFs
    filtered_vcfs = Channel.empty()
    filtered_tbis = Channel.empty()

    if (paired_samples) {
        filtered_vcfs = filtered_vcfs.mix(paired_filtered_vcfs)
        filtered_tbis = filtered_tbis.mix(paired_filtered_tbis)
    }

    if (unpaired_samples) {
        filtered_vcfs = filtered_vcfs.mix(unpaired_filtered_vcfs)
        filtered_tbis = filtered_tbis.mix(unpaired_filtered_tbis)
    }

    // Normalize the filtered VCFs
    vcfs_for_norm = filtered_vcfs.join(filtered_tbis, by: [0])
    
    BCFTOOLS_NORM(
        vcfs_for_norm,
        fasta
    )

    // Filter for PASS variants
    BCFTOOLS_VIEW(
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi, by: [0])
    )
    
    GATK4_FUNCOTATOR(
            BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi, by: [0]),
            fasta,
            fai,
            dict,
            funcotator_resources,
            intervals
        )

    // Annotate with Annovar if requested
    ANNOVAR(
            BCFTOOLS_VIEW.out.vcf,
            fasta,
            annovar_db
        )
    
}


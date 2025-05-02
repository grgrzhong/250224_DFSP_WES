/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                Mutect2 SNV/Indels calling workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load required modules
include { PREPARE_GENOME                                                                    } from '../../subworkflows/mutation_calling/prepare_genome.nf'
include { GATK4_GETPILEUPSUMMARIES as  GATK4_GETPILEUPSUMMARIES_PAIRED_TUMOUR               } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as  GATK4_GETPILEUPSUMMARIES_PAIRED_NORMAL               } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_GETPILEUPSUMMARIES as  GATK4_GETPILEUPSUMMARIES_TUMOUR_ONLY                 } from '../../modules/variant_calling/gatk4/getpileupsummaries'
include { GATK4_CALCULATECONTAMINATION as GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL        } from '../../modules/variant_calling/gatk4/calculatecontamination'
include { GATK4_CALCULATECONTAMINATION as GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY          } from '../../modules/variant_calling/gatk4/calculatecontamination'
include { GATK4_MUTECT2 as GATK4_MUTECT2_TUMOR_NORMAL                                       } from '../../modules/variant_calling/gatk4/mutect2'
include { GATK4_MUTECT2 as GATK4_MUTECT2_TUMOR_ONLY                                         } from '../../modules/variant_calling/gatk4/mutect2'
include { GATK4_LEARNREADORIENTATIONMODEL as GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_NORMAL  } from '../../modules/variant_calling/gatk4/learnreadorientationmodel'
include { GATK4_LEARNREADORIENTATIONMODEL as GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_ONLY    } from '../../modules/variant_calling/gatk4/learnreadorientationmodel'
include { GATK4_FILTERMUTECTCALLS                                                           } from '../../modules/variant_calling/gatk4/filtermutectcalls'
include { BCFTOOLS_NORM as MUTECT2_BCFTOOLS_NORM                                            } from '../../modules/variant_calling/bcftools/norm'
include { BCFTOOLS_VIEW as MUTECT2_BCFTOOLS_VIEW                                            } from '../../modules/variant_calling/bcftools/view'
include { TABIX as MUTECT2_TABIX                                                            } from "../../modules/variant_calling/tabix/tabix"

// input
params.input = "/home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test2.csv"
// params.input = "/home/zhonggr/projects/250224_DFSP_WES/data/test_data/csv/test.csv"

workflow MUTECT2_CALL{

    
    main:
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
        // funcotator_resources = PREPARE_GENOME.out.funcotator_resources
        // annovar_db           = PREPARE_GENOME.out.annovar_db

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
            paired_tumour_samples = paired_samples
                .map { 
                    meta, tumour_bam, tumour_bai, normal_bam, normal_bai ->
                    [meta, tumour_bam, tumour_bai]
            }
            paired_tumour_samples.view()

            // Get pileup summaries for tumor samples in paired mode
            GATK4_GETPILEUPSUMMARIES_PAIRED_TUMOUR(
                paired_tumour_samples,
                pileup_variants,
                pileup_variants_tbi
            )

            // Extract normal samples from paired samples for pileup
            paired_normal_samples = paired_samples
                .map { meta, tumour_bam, tumour_bai, normal_bam, normal_bai ->
                    def normal_meta = meta.clone()
                    normal_meta.tumour_id = normal_meta.normal_id
                    [normal_meta, normal_bam, normal_bai]
            }
            paired_normal_samples.view()

            GATK4_GETPILEUPSUMMARIES_PAIRED_NORMAL(
                paired_normal_samples,
                pileup_variants,
                pileup_variants_tbi
            )

            // Calculate contamination for paired samples
            GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL(
                GATK4_GETPILEUPSUMMARIES_PAIRED_TUMOUR.out.table
                .join(GATK4_GETPILEUPSUMMARIES_PAIRED_NORMAL.out.table, by: [0])
                .map { meta, tumor_table, normal_table ->
                    [meta, tumor_table, normal_table]
                }
            )

            // Run Mutect2 for paired samples
            GATK4_MUTECT2_TUMOR_NORMAL(
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

            GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_NORMAL(
                GATK4_MUTECT2_TUMOR_NORMAL.out.f1r2
            )

        }

        /*
            =================== TUMOR-ONLY UNPAIRED ANALYSIS ====================
        */
        if (unpaired_samples) {
        
            // Extract tumor samples for unpaired analysis
            unpaired_tumour_samples = unpaired_samples
                .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                [meta, tumor_bam, tumor_bai]
            }
            unpaired_tumour_samples.view()

            // Get pileup summaries for tumor-only samples
            GATK4_GETPILEUPSUMMARIES_TUMOUR_ONLY(
                unpaired_tumour_samples,
                pileup_variants,
                pileup_variants_tbi
            )

            // Calculate contamination for tumor-only samples
            GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY(
                GATK4_GETPILEUPSUMMARIES_TUMOUR_ONLY.out.table
                .map { meta, table ->
                    [meta, table, []]  // Empty file for matched normal
                }
            )

            // Run Mutect2 for tumor-only samples
            GATK4_MUTECT2_TUMOR_ONLY(
                unpaired_tumour_samples,
                unpaired_samples.map { 
                    meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
                    [meta, [], []]
                },
                fasta,
                fai,
                dict,
                germline_resource,
                germline_resource_tbi,
                panel_of_normals,
                panel_of_normals_tbi,
                intervals
            )

            GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_ONLY(
                GATK4_MUTECT2_TUMOR_ONLY.out.f1r2
            )
        }
        
        // Combine all Mutect2 outputs for further processing
        mutect2_read_orientation_models = GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_NORMAL.out.artifactprior
            .mix(GATK4_LEARNREADORIENTATIONMODEL_TUMOUR_ONLY.out.artifactprior)

        mutect2_vcf = GATK4_MUTECT2_TUMOR_NORMAL.out.vcf.mix(GATK4_MUTECT2_TUMOR_ONLY.out.vcf)

        mutect2_tbi = GATK4_MUTECT2_TUMOR_NORMAL.out.tbi.mix(GATK4_MUTECT2_TUMOR_ONLY.out.tbi)

        mutect2_stats = GATK4_MUTECT2_TUMOR_NORMAL.out.stats.mix(GATK4_MUTECT2_TUMOR_ONLY.out.stats)

        // Combine contamination tables
        mutect2_contamination_tables = GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL.out.contamination
            .mix(GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY.out.contamination)
        
        // Combine segmentation tables
        mutect2_segmentation_tables = GATK4_CALCULATECONTAMINATION_TUMOUR_NORMAL.out.segmentation
            .mix(GATK4_CALCULATECONTAMINATION_TUMOUR_ONLY.out.segmentation)
        
        // Prepare inputs for FilterMutectCalls
        filter_input = mutect2_vcf
            .join(mutect2_tbi, by: [0])
            .join(mutect2_stats, by: [0])
            .join(mutect2_read_orientation_models, by: [0])
            .join(mutect2_contamination_tables, by: [0])
            .join(mutect2_segmentation_tables, by: [0])
            .map { meta, vcf, tbi, stats, orientation_model, contamination, segmentation ->
                [meta, vcf, tbi, stats, orientation_model, contamination, segmentation]
            }
        
        filter_input.view()

        // Filter Mutect2 calls
        GATK4_FILTERMUTECTCALLS(
            filter_input,
            fasta,
            fai,
            dict
        )
        
        // Normalize variants with bcftools
        MUTECT2_BCFTOOLS_NORM(
            GATK4_FILTERMUTECTCALLS.out.vcf.join(GATK4_FILTERMUTECTCALLS.out.tbi, by: [0]),
            fasta,
            fai,
            dict
        )
        
        // Filter for PASS variants with bcftools
        MUTECT2_BCFTOOLS_VIEW(MUTECT2_BCFTOOLS_NORM.out.vcf)
        
        MUTECT2_TABIX(MUTECT2_BCFTOOLS_VIEW.out.vcf)

    // Define workflow outputs
    emit:
        vcf = MUTECT2_BCFTOOLS_VIEW.out.vcf
        tbi = MUTECT2_TABIX.out.tbi
}

workflow  {
    // Run the Mutect2 calling workflow
    MUTECT2_CALL()
    
}
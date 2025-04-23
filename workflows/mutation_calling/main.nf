#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

include { FASTP_TRIM                               } from '../../modules/variant_calling/fastp/main.nf'
include { FASTQC                                   } from '../../modules/variant_calling/fastqc/main.nf'
include { BWA_MEM                                  } from '../../modules/variant_calling/bwa/mem/main.nf'
include { TAG_UMI                                  } from '../../modules/variant_calling/tagumi/main.nf'
include { SAMTOOLS_SORT                            } from '../../modules/variant_calling/samtools/sort/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP } from '../../modules/variant_calling/samtools/index/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RECAL   } from '../../modules/variant_calling/samtools/index/main.nf'
include { GATK4_MARKDUPLICATES                     } from '../../modules/variant_calling/gatk4/markduplicates/main.nf'
include { SAMTOOLS_INDEX                           } from '../../modules/variant_calling/samtools/index/main.nf'
include { GATK4_BASERECALIBRATOR                   } from '../../modules/variant_calling/gatk4/baserecalibrator/main.nf'
include { GATK4_APPLYBQSR                          } from '../../modules/variant_calling/gatk4/applybqsr/main.nf'
include { GATK4_COLLECTHSMETRICS                   } from '../../modules/variant_calling/gatk4/collecthsmetrics/main.nf'
include { BAMTOOLS_STATS                           } from '../../modules/variant_calling/bamtools/stats/main.nf'

/*
========================================================================================
    PARAMETERS
========================================================================================
*/

// Set default input and output values
params.input            = params.input ?: "/home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test.csv"
params.outdir           = params.outdir ?: "${launchDir}/results"
params.publish_dir_mode = params.publish_dir_mode ?: "copy"


/*
========================================================================================
    WORKFLOW
========================================================================================
*/

// Define preprocessing workflow as a subworkflow
workflow {
    // Initialize version tracking
    // ch_versions = Channel.empty()

    // Setup reference channels
    fasta = file(params.genomes[params.genome].fasta)
    fai = file(params.genomes[params.genome].fai)
    dict = file(params.genomes[params.genome].dict)
    dbsnp = file(params.genomes[params.genome].dbsnp)
    dbsnp_tbi = file(params.genomes[params.genome].dbsnp_tbi)
    bait_intervals = file(params.genomes[params.genome].bait_intervals)
    target_intervals = file(params.genomes[params.genome].target_intervals)
    intervals = file(params.genomes[params.genome].intervals)

    // Verify essential reference files exist
    if (!fasta || !fasta.exists()) {
        exit(
            1,
            "ERROR: fasta reference file not found: " + "${params.genomes[params.genome]?.fasta ?: 'Not specified'}",
        )
    }
    if (!fai || !fai.exists()) {
        exit(
            1,
            "ERROR: fasta.fai file not found: " + "${params.genomes[params.genome]?.fai ?: 'Not specified'}",
        )
    }
    if (!dict || !dict.exists()) {
        exit(
            1,
            "ERROR: Dict file not found: " + "${params.genomes[params.genome]?.dict ?: 'Not specified'}",
        )
    }

    // Check for other required files
    if (!intervals || !intervals.exists()) {
        log.warn(
            "WARNING: Intervals file not found: " + "${params.genomes[params.genome]?.intervals ?: 'Not specified'}"
        )
    }
    if (!dbsnp || !dbsnp.exists()) {
        log.warn(
            "WARNING: dbSNP file not found: " + "${params.genomes[params.genome]?.dbsnp ?: 'Not specified'}"
        )
    }

    if (!bait_intervals || !bait_intervals.exists() || !target_intervals || !target_intervals.exists()) {
        log.warn(
            "WARNING: Bait or target intervals not found. " + "Hybrid metrics may fail."
        )
    }

    // Create BWA index channel with all index files
    fasta_index = Channel.fromPath("${fasta}.{amb,ann,bwt,pac,sa}", checkIfExists: true)
        .ifEmpty {
            exit(1, "ERROR: BWA indices not found for: ${fasta}")
        }
        .collect()


    // Create input sample channel
    Channel.fromPath(params.input)
        .ifEmpty { exit(1, "Sample sheet not found at ${params.input}") }
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            meta.patient_id = row.patient
            meta.status = row.status.toInteger()

            // Check that the fastq files exist and are different
            def fastq_1 = file(row.fastq_1)
            def fastq_2 = file(row.fastq_2)

            if (!fastq_1.exists()) {
                exit(1, "ERROR: Read 1 fastq file not found: ${row.fastq_1}")
            }
            if (!fastq_2.exists()) {
                exit(1, "ERROR: Read 2 fastq file not found: ${row.fastq_2}")
            }

            return [meta, fastq_1, fastq_2]
        }
        .set { ch_reads }

    // Trim reads and extract UMI information
    FASTP_TRIM(ch_reads)

    // Quality control on trimmed reads
    FASTQC(FASTP_TRIM.out.reads)

    // Align reads to reference genome
    BWA_MEM(
        FASTP_TRIM.out.reads,
        fasta,
        fasta_index,
    )

    // Tag UMIs in alignments
    TAG_UMI(BWA_MEM.out.bam)

    // Sort BAM by coordinates
    SAMTOOLS_SORT(TAG_UMI.out.bam)

    // Mark duplicates and creat index
    GATK4_MARKDUPLICATES(SAMTOOLS_SORT.out.bam)

    // Create bam index
    SAMTOOLS_INDEX_MARKDUP(GATK4_MARKDUPLICATES.out.bam)

    // Base recalibration
    GATK4_BASERECALIBRATOR(
        GATK4_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX_MARKDUP.out.bai),
        fasta,
        fai,
        dict,
        dbsnp,
        dbsnp_tbi,
        intervals,
    )

    // Apply BQSR to generate recalibrated BAM
    GATK4_APPLYBQSR(
        GATK4_MARKDUPLICATES.out.bam.join(SAMTOOLS_INDEX_MARKDUP.out.bai),
        GATK4_BASERECALIBRATOR.out.table.map { meta, table -> table },
        intervals,
    )

    // Create index for the recalibrated BAM
    SAMTOOLS_INDEX_RECAL(GATK4_APPLYBQSR.out.bam)

    // Calculate hybrid selection metrics
    GATK4_COLLECTHSMETRICS(
        GATK4_APPLYBQSR.out.bam.join(SAMTOOLS_INDEX_RECAL.out.bai),
        fasta,
        fai,
        dict,
        bait_intervals,
        target_intervals,
    )

    // Generate BAM statistics
    BAMTOOLS_STATS(
        GATK4_APPLYBQSR.out.bam.join(SAMTOOLS_INDEX_RECAL.out.bai)
    )
}

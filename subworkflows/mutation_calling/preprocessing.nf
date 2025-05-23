
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                Preprocesing raw fastq files for variant calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// IMPORT MODULES/SUBWORKFLOWS
include { FASTP_TRIM                               } from '../../modules/local/fastp/main.nf'
include { FASTQC                                   } from '../../modules/local/fastqc'
include { BWA_MEM                                  } from '../../modules/local/bwa/mem'
include { TAG_UMI                                  } from '../../modules/local/tagumi'
include { SAMTOOLS_SORT                            } from '../../modules/local/samtools/sort'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MARKDUP } from '../../modules/local/samtools/index'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_RECAL   } from '../../modules/local/samtools/index'
include { SAMTOOLS_INDEX                           } from '../../modules/local/samtools/index'
include { GATK4_MARKDUPLICATES                     } from '../../modules/local/gatk4/markduplicates'
include { GATK4_BASERECALIBRATOR                   } from '../../modules/local/gatk4/baserecalibrator'
include { GATK4_APPLYBQSR                          } from '../../modules/local/gatk4/applybqsr'
include { GATK4_COLLECTHSMETRICS                   } from '../../modules/local/gatk4/collecthsmetrics'
include { BAMTOOLS_STATS                           } from '../../modules/local/bamtools/stats'

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

// Define preprocessing workflow as a subworkflow
workflow PREPROCESSING {

    take:    
        ch_reads      // Input reads channel
        fasta
        fai           
        dict
        dbsnp
        dbsnp_tbi
        bait_intervals
        target_intervals
        intervals
        bwa_index     // BWA index files

    main:
    // Trim reads and extract UMI information
    FASTP_TRIM(ch_reads)

    // Quality control on trimmed reads
    FASTQC(FASTP_TRIM.out.reads)

    // Align reads to reference genome
    BWA_MEM(
        FASTP_TRIM.out.reads,
        fasta,
        bwa_index,
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
    
    emit:
        bam = GATK4_APPLYBQSR.out.bam
        bai = GATK4_APPLYBQSR.out.bai
    
}

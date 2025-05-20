
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                Copy number variation analysis using sequenza subworkflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Load required modules
include { SEQUENZAUTILS_BAM2SEQZ      } from "../../modules/local/sequenzautils/bam2seqz/main.nf"
include { SEQUENZAUTILS_SEQZBINNING   } from "../../modules/local/sequenzautils/seqzbinning/main.nf"

workflow CNV_SEQUENZA {

    take:
    bam_tumour_normal
    fasta
    wigfile
    window_size

    main:
    
    // Convert BAM to SEQZ format using sequenza-utils
    SEQUENZAUTILS_BAM2SEQZ(
        bam_tumour_normal,
        fasta,
        wigfile
    )
    // 
    SEQUENZAUTILS_SEQZBINNING(
        SEQUENZAUTILS_BAM2SEQZ.out.seqz
        .join(SEQUENZAUTILS_BAM2SEQZ.out.tbi),
        window_size
    )

    emit:
    seqz = SEQUENZAUTILS_BAM2SEQZ.out.seqz
}


#!/bin/bash
#SBATCH --job-name=RNA_STARFusion
#SBATCH --partition=amd
#SBATCH --time=48:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.out
#SBATCH --error=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Setup the working directory
# cd /home/zhonggr/projects/250224_DFSP_WES

# # Global nextflow configuration
# # export NXF_OFFLINE=true
# # export NXF_HOME=$HOME/.nextflow
# export NXF_DISABLE_CHECK_LATEST=true
# export NXF_OPTS="-Xms512m -Xmx8g"
# # export NXF_LOG_FILE="${PWD}/.nextflow.log"
# rm -f .nextflow.log*

# # Run the Nextflow pipeline in local mode
# nextflow run workflows/rna_fusion.nf \
#     -profile local

# Run the Nextflow pipeline in hpc mode
# nextflow run workflows/somatic_variant_calling.nf \
#     -profile hpc \
#     -resume \
#     --input /home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test1.csv

## Install STAR-Fusion
# conda create -n starfusion python=2.7
# conda activate starfusion
# conda install -c bioconda star-fusion
# conda install -c bioconda trinity
# conda install -c conda-forge -c bioconda samtools bzip2

# conda activate starfusion

## Reference directories in multiomics
# conda activate starfusion
# export TRINITY_HOME=/home/zhonggr/miniforge3/envs/starfusion/opt/trinity-2.8.5

singularity shell \
    --bind /mnt/m/Reference:/mnt/m/Reference \
    --bind /mnt/m/RNA-seq/STUMP/Input-trimmed:/mnt/m/RNA-seq/STUMP/Input-trimmed \
    --bind /mnt/f/projects/250224_DFSP_WES/outputs/stump/STAR-Fusion:/mnt/f/projects/250224_DFSP_WES/outputs/stump/STAR-Fusion \
    /home/zhonggr/projects/250224_DFSP_WES/containers/singularity/star-fusion.v1.15.0.simg

ref_dir=/mnt/m/Reference
STARINDEX=${ref_dir}/Gencode/STAR_index/
STARINDEX_HG19=${ref_dir}/Gencode/STAR_index_hg19/
REFERENCE=${ref_dir}/Gencode/gencode.hg38.v36.primary_assembly.fa
ANNOTATION=${ref_dir}/Gencode/gencode.v36.primary_assembly.annotation.gtf
# CTAT_RESOURCE_LIB=${ref_dir}/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/
CTAT_RESOURCE_LIB=${ref_dir}/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/
INPUT=/mnt/m/RNA-seq/STUMP/Input-trimmed

# Check if required directories and files exist
echo "Checking required directories and files..."

# Check reference directories
if [ ! -d "$ref_dir" ]; then
    echo "ERROR: Reference directory $ref_dir does not exist"
    exit 1
fi

if [ ! -d "$STARINDEX" ]; then
    echo "ERROR: STAR index directory $STARINDEX does not exist"
    exit 1
fi

# Check reference files
if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference genome file $REFERENCE does not exist"
    exit 1
fi

if [ ! -f "$ANNOTATION" ]; then
    echo "ERROR: Annotation file $ANNOTATION does not exist"
    exit 1
fi

if [ ! -d "$INPUT" ]; then
    echo "ERROR: Input directory $INPUT does not exist"
    exit 1
fi

if [ ! -d "$CTAT_RESOURCE_LIB" ]; then
    echo "ERROR: CTAT resource library $CTAT_RESOURCE_LIB does not exist" >&2
    exit 1
fi


# Check if input files exist
input_files=$(ls ${INPUT}/*.fastq.gz 2>/dev/null | grep "R1")
output_dir=/mnt/f/projects/250224_DFSP_WES/outputs/stump/STAR-Fusion

# singularity exec -e -B `pwd` -B /path/to/ctat_genome_lib_build_dir \
#         star-fusion-v$version.simg \
#         STAR-Fusion \
#         --left_fq reads_1.fq.gz \
#         --right_fq reads_2.fq.gz \
#         --genome_lib_dir /path/to/ctat_genome_lib_build_dir \
#         -O StarFusionOut \
#         --FusionInspector validate \
#         --examine_coding_effect \
#         --denovo_reconstruct

file=/mnt/m/RNA-seq/STUMP/Input-trimmed/S25_trimmed_R1.fastq.gz

for file in $(ls ${INPUT}/*.fastq.gz | grep "R1"); do 

    echo $file; 
    
    FILENAME=$(basename $file | cut -d "_" -f 1); 
    
    echo $FILENAME; 
    
    OUTPUT=${output_dir}/${FILENAME}/
    mkdir -p $OUTPUT
    echo $OUTPUT;

    STAR --genomeDir $STARINDEX \
        --readFilesIn $file ${file//R1/R2} \
        --outReadsUnmapped None \
        --runThreadN 8 \
        --twopassMode Basic \
        --readFilesCommand "gunzip -c" \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 8 \
        --chimOutJunctionFormat 1 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMattrRGline ID:GRPundef SM:$FILENAME \
        --chimMultimapScoreRange 3 \
        --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 \
        --chimNonchimScoreDropMin 10 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
        --outFileNamePrefix $OUTPUT \
        --outSAMtype BAM SortedByCoordinate \
        --outTmpDir /tmp/STAR_${FILENAME}/ \
        --quantMode GeneCounts \
        >& ${OUTPUT}/STAR.log

    # bam file needs SM tag for CTAT mutation

    # # --no_filter \
    STAR-Fusion --genome_lib_dir $CTAT_RESOURCE_LIB \
                -J $OUTPUT/Chimeric.out.junction \
                --left_fq $file\
                --right_fq ${file//R1/R2}\
                --output_dir $OUTPUT \
                --examine_coding_effect \
                --extract_fusion_reads \
                --FusionInspector inspect \
                >& ${OUTPUT}/STAR-Fusion.log

    samtools index ${OUTPUT}/Aligned.sortedByCoord.out.bam
            
    ## Set minimum FFPM to filter; Default 0.1; set to 0 to turn off
    ## --min_FFPM 100
    ##--max_sensitivity    includes options: --min_junction_reads 0 --min_sum_frags 1 --require_LDAS 0 --min_spanning_frags_only 1 --min_novel_junction_support 1 --skip_FFPM --no_single_fusion_per_breakpoint --skip_EM

    #mv $file ${file//R1/R2} $(dirname $file)/../Output/$FILENAME/

done

##--chimSegmentMin default: 0; int>=0: minimum length of chimeric segment length, if ==0, no chimeric output

##--chimJunctionOverhangMin default: 20; int>=0: minimum overhang for a chimeric junction

##--chimOutJunctionFormat; default: 0; int: formatting type for the Chimeric.out.junction file; 
##0 no comment lines/headers 1 comment lines at the end of the file: command line and Nreads: total, unique/multi-mapping

##--alignSJDBoverhangMin default: 3 int>0: minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments

##--alignSJoverhangMin default: 5 int>0: minimum overhang (i.e. block size) for spliced alignments

##--alignMatesGapMax default: 0 maximum gap between two mates, if 0, max intron gap will be determined by (2ˆwinBinNbits)*winAnchorDistNbins

##--alignIntronMax default: 0 maximum intron size, if 0, max intron size will be determined by (2ˆwinBinNbits)*winAnchorDistNbins

##--alignSJstitchMismatchNmax default: 0 -1 0 0; 4*int>=0: maximum number of mismatches for stitching of the splice junctions (-1: no limit). (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.

## --chimMultimapScoreRange default: 1 int>=0: the score range for multi-mapping chimeras below the best chimeric score. Only works with –chimMultimapNmax > 1

## --chimScoreJunctionNonGTAG default: -1 int: penalty for a non-GT/AG chimeric junction

##--chimMultimapNmax default: 0 int>=0: maximum number of chimeric multi-alignments 0 use the old scheme for chimeric detection which only considered unique alignments

## --chimNonchimScoreDropMin default: 20 52int>=0: to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value


### --outSAMstrandField intronMotif \ # include for potential use with StringTie for assembly
### --chimSegmentMin 12 \ # ** essential to invoke chimeric read detection & reporting **
### --chimOutJunctionFormat 1  \ # **essential** includes required metadata in Chimeric.junction.out file.
### --alignMatesGapMax 100000 \ # avoid readthru fusions within 100k
### --alignSJstitchMismatchNmax 5 -1 5 5 \  # settings improved certain chimera detections

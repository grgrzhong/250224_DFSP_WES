#!/bin/bash
#SBATCH --job-name=Test_NF
#SBATCH --partition=amd
#SBATCH --time=48:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/logs/slurm/%x_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Setup the working directory
cd /home/zhonggr/projects/250224_DFSP_WES

# Global nextflow configuration
# export NXF_OFFLINE=true
# export NXF_HOME=$HOME/.nextflow
export NXF_DISABLE_CHECK_LATEST=true
export NFX_OPTS="-Xms=512m -Xmx=8g"
export NXF_ANSI_LOG=true

# export NXF_LOG_FILE="${PWD}/logs/nextflow/.nextflow.log"
export NXF_LOG_FILE="${PWD}/.nextflow.log"
# rm -f ${NXF_LOG_FILE}

# Test the mutation_calling workflow
nextflow run workflows/mutation_calling/main.nf \
    -profile hpc \
    -resume \
    --input /home/zhonggr/projects/250224_DFSP_WES/data/WES/csv/test.csv \
    --outdir results
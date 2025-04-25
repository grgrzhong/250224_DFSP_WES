#!/bin/bash
#SBATCH --job-name=Test_NF
#SBATCH --partition=amd
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/$(date +%Y%m%d)_%j_%x.out
#SBATCH --error=/home/zhonggr/projects/250224_DFSP_WES/slurm/$(date +%Y%m%d)_%j_%x.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Setup the working directory
cd /home/zhonggr/projects/250224_DFSP_WES

# Global nextflow configuration
# export NXF_OFFLINE=true
# export NXF_HOME=$HOME/.nextflow
export NXF_DISABLE_CHECK_LATEST=true
export NXF_OPTS="-Xms512m -Xmx8g"
# export NXF_ANSI_LOG=true

# export NXF_LOG_FILE="${PWD}/logs/nextflow/.nextflow.log"
export NXF_LOG_FILE="${PWD}/.nextflow.log"
rm -f ${NXF_LOG_FILE}

# Test the mutation_calling workflow
nextflow run workflows/mutation_calling/main.nf \
    -profile hpc \
    --input /home/zhonggr/projects/250224_DFSP_WES/data/sarc/csv/samplesheet.csv \
    --outdir data/sarc

# run testing
rm -f ${NXF_LOG_FILE}
nextflow run subworkflows/mutation_calling/mutect2_call.nf \
    -profile hpc \
    -work-dir test_work \
    --outdir test_results

#!/bin/bash
#SBATCH --job-name=Test_NF
#SBATCH --partition=amd
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/slurm/%x_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zhonggr@hku.hk

# Setup the working directory
cd /home/zhonggr/projects/250224_DFSP_WES

# Global nextflow configuration
# export NXF_OFFLINE=true
# export NXF_HOME=$HOME/.nextflow
export NXF_DISABLE_CHECK_LATEST=true
export NXF_OPTS="-Xms512m -Xmx8g"
export NXF_LOG_FILE="${PWD}/.nextflow.log"
rm -f .nextflow.log*

nextflow run subworkflows/mutation_calling/mutect2_call.nf \
    -profile local \
    -resume \
    --input /home/zhonggr/projects/250224_DFSP_WES/data/wes/csv/test1.csv
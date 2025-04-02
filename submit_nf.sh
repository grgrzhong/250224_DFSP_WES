#!/bin/bash
#SBATCH --job-name=CNV_FACETS            
#SBATCH --partition=amd                  
#SBATCH --time=48:00:00                  
#SBATCH --qos=normal                     
#SBATCH --nodes=1                        
#SBATCH --ntasks=1                       
#SBATCH --cpus-per-task=32               
#SBATCH --mem=128G                       
#SBATCH --output=/home/zhonggr/projects/250224_DFSP_WES/logs/%x_%j.out 
#SBATCH --error=/home/zhonggr/projects/250224_DFSP_WES/logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL             
#SBATCH --mail-user=zhonggr@hku.hk

export NXF_OFFLINE=true
export NXF_HOME=$HOME/.nextflow
export NXF_DISABLE_CHECK_LATEST=true

# Set the working directory
cd /home/zhonggr/projects/250224_DFSP_WES

# Run the cnv_facets workflow
nextflow run workflows/mutation_calling/facets.nf \
    --outdir results \
    --samplesheet /home/zhonggr/projects/250224_DFSP_WES/data/WES/DFSP/samplesheet.csv
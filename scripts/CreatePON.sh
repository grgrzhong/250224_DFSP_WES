#!/bin/bash

# Use source to activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate varcall

work_dir="/home/zhonggr/projects/250224_DFSP_WES/data/mfs"
ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
REFERENCE=/home/zhonggr/projects/250224_DFSP_WES/reference/Gencode/gencode.hg38.v36.primary_assembly.fa
GERMLINE=/home/zhonggr/projects/250224_DFSP_WES/reference/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
ANNOTATION_FILE=/home/zhonggr//projects/250224_DFSP_WES/reference/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/
INTERVAL=/home/zhonggr/projects/250224_DFSP_WES/reference/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed
BAM=/home/zhonggr/projects/250224_DFSP_WES
PON_OUT=/media/maximus/Data/WES/PON-Mutect

export PATH=/media/maximus/Data/Software/gatk-4.4.0.0/:$PATH

CASELIST=$(ls $BAM | grep N)

#Step 1. Run Mutect2 in tumor-only mode for each normal sample.
for sample in $CASELIST; do 
echo $sample

gatk Mutect2 \
-R $REFERENCE \
-I $BAM/$sample/${sample}_recalibrated.bam \
-max-mnp-distance 0 \
-L $INTERVAL \
-O $PON_OUT/${sample}_pon.vcf.gz

done

#Step 2. Create a GenomicsDB from the normal Mutect2 calls.
cmd='gatk --java-options "-Xmx64g" GenomicsDBImport
-R $REFERENCE
-L $INTERVAL
--genomicsdb-workspace-path $PON_OUT/pon_db
--merge-input-intervals true'

for vcf_file in $(ls $PON_OUT/*.vcf.gz); do
 cmd+=" -V $vcf_file"
done

eval $cmd

#Step 3. Combine the normal calls using CreateSomaticPanelOfNormals.
#Take ~2 days
gatk --java-options -Xmx108g CreateSomaticPanelOfNormals \
-R $REFERENCE \
-L $INTERVAL \
-V gendb://$PON_OUT/pon_db \
-O $PON_OUT/pon.vcf.gz


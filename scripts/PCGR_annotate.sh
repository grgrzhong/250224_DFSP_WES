#!/bin/bash

SERVER=/run/user/1000/gvfs/smb-share:server=my-nas.local,share=multiomics
REFERENCE=/media/maximus/Data/Reference/Gencode/gencode.hg38.v36.primary_assembly.fa
GERMLINE=/media/maximus/Data/Reference/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
VEP_DIR=$SERVER/Reference/VEP_cache
REFDATA_DIR=$SERVER/Reference/PCGR_reference/20250314
INTERVAL=/media/maximus/Data/Reference/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed
BAM=/media/maximus/Data/WES/BAM
MUTECT=/media/maximus/Data/WES/Mutect-Call-APYL
PON=$SERVER/WES/DFSP/PON-Mutect/pon.vcf.gz


cd ~

# To activate this environment, use
conda activate /home/maximus/pcgr_conda/pcgr

SAMPLELIST=$(ls $MUTECT)

for SAMPLE in $SAMPLELIST; do 
echo $SAMPLE

python "$SERVER/Scripts/DNA analysis/Reformat_vcf.py" \
-I $MUTECT/${SAMPLE}/${SAMPLE}_normalized_filtered.vcf.gz \
-O $MUTECT/${SAMPLE}/${SAMPLE}_normalized_filtered_reformatted.vcf.gz 

tabix $MUTECT/${SAMPLE}/${SAMPLE}_normalized_filtered_reformatted.vcf.gz 

INPUT_VCF=$MUTECT/${SAMPLE}/${SAMPLE}_normalized_filtered_reformatted.vcf.gz 
mkdir -p $MUTECT/../PCGR/${SAMPLE}
OUTPUT=$MUTECT/../PCGR/${SAMPLE}

pcgr \
--input_vcf $INPUT_VCF \
--vep_dir $VEP_DIR \
--refdata_dir $REFDATA_DIR \
--output_dir $OUTPUT \
--genome_assembly grch38 \
--sample_id ${SAMPLE} \
--assay WES \
--effective_target_size_mb 34 \
--tumor_dp_tag TDP \
--tumor_af_tag TAF \
--control_dp_tag NDP \
--control_af_tag NAF \
--tumor_dp_min 20 \
--tumor_af_min 0.05 \
--control_dp_min 10 \
--control_af_max 0.01 \
--estimate_tmb \
--tmb_dp_min 20 \
--tmb_af_min 0.05 \
--estimate_msi \
--estimate_signatures \
--vcf2maf \
--ignore_noncoding \
--force_overwrite


##Tumour only calling
pcgr \
--input_vcf $INPUT_VCF \
--vep_dir $VEP_DIR \
--refdata_dir $REFDATA_DIR \
--output_dir $OUTPUT \
--genome_assembly grch38 \
--sample_id ${SAMPLE} \
--assay WES \
--effective_target_size_mb 34 \
--tumor_only \
--tumor_dp_tag TDP \
--tumor_af_tag TAF \
--tumor_dp_min 20 \
--tumor_af_min 0.05 \
--estimate_tmb \
--tmb_dp_min 20 \
--tmb_af_min 0.05 \
--estimate_msi \
--estimate_signatures \
--vcf2maf \
--ignore_noncoding \
--force_overwrite


done

#--input_cna ${INPUT_PATH}/${SAMPLE}_PCGR.seg.txt \
#--n_copy_gain 3 \

# To deactivate an active environment, use
conda deactivate



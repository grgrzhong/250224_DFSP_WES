#!/bin/bash

# Use source to activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate varcall

ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"
bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
vcf_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2"
work_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_filter"
mkdir -p ${work_dir}

# Sort the repeatmasker bed file and create index
# gunzip ${ref_dir}/RepeatMasker.bed.gz
# bgzip ${ref_dir}/RepeatMasker.bed
# sort -k1,1 -k2,2n ${ref_dir}/RepeatMasker.bed -o ${ref_dir}/RepeatMasker.sorted.bed
# bgzip -c ${ref_dir}/RepeatMasker.sorted.bed > ${ref_dir}/RepeatMasker.sorted.bed.gz
# tabix -p bed ${ref_dir}/RepeatMasker.sorted.bed.gz

# # Sort the blacklist bed file and create index
# gunzip ${ref_dir}/blacklist.bed.gz
# bgzip ${ref_dir}/blacklist.bed
# sort -k1,1 -k2,2n ${ref_dir}/blacklist.bed -o ${ref_dir}/blacklist.sorted.bed
# bgzip -c ${ref_dir}/blacklist.sorted.bed > ${ref_dir}/blacklist.sorted.bed.gz
# tabix -p bed ${ref_dir}/blacklist.sorted.bed.gz

REFERENCE=${ref_dir}/Gencode/gencode.hg38.v36.primary_assembly.fa
GERMLINE=${ref_dir}/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
ANNOTATION_FILE=${ref_dir}/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/
INTERVAL=${ref_dir}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed
PON=${ref_dir}/pon_dfsp/pon.vcf.gz

echo "Reference directoyr:  ${ref_dir}"
echo "Bam directory:        ${bam_dir}"
echo "Vcf directory:        ${vcf_dir}"
echo "Work directory:       ${work_dir}"
echo "Reference:            ${REFERENCE}"
echo "Germline:             ${GERMLINE}"
echo "Annotation file:      ${ANNOTATION_FILE}"
echo "Interval:             ${INTERVAL}"
echo "PON:                  ${PON}"

# List all unfiltered.vcf.gz files in the directory
# tumour_ids=$(find "$vcf_dir" -name "*_unfiltered.vcf.gz" | sort | sed 's|.*/||' | sed 's/_unfiltered.vcf.gz$//')

# find "$vcf_dir" -name "*_unfiltered.vcf.gz" | sed -E 's|.*/([A-Z]+-[0-9]+)_.*|\1|' | sort -u
# CASELIST=$(ls $vcf_dir | cut -d '-' -f1,2 | uniq)

echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > ${work_dir}/vcf.rm.header
echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\"EncodeDacMapability\">" > ${work_dir}/vcf.map.header

tumour_ids=DFSP-028-T

for tumour_id in $tumour_ids; do 
    
    # Extract case_id from tumour_id (everything before the second hyphen)
    case_id=$(echo $tumour_id | cut -d'-' -f1,2)
    
    echo ${case_id}
    
    echo $tumour_id
    normal_id=${case_id}-N

    echo ${normal_id}

    # Check if the normal sample directory exists
    mkdir -p ${work_dir}/${tumour_id}
    VAROUT=${work_dir}/${tumour_id}

    # Get Pileup Summaries
    echo $(date +"%F") $(date +"%T") "Getting Pileup Summaries of tumour sample ..."
    gatk GetPileupSummaries \
        -I $bam_dir/$tumour_id/${tumour_id}_recalibrated.bam \
        -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
        -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
        -O ${VAROUT}/${tumour_id}.getpileupsummaries.table

    #Check if presence of paired normal samples
    if [ -d ${bam_dir}/${normal_id} ]; then

        echo $(date +"%F") $(date +"%T") "Getting Pileup Summaries of normal sample ..."
        gatk GetPileupSummaries \
            -I $bam_dir/$normal_id/${normal_id}_recalibrated.bam \
            -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
            -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
            -O ${VAROUT}/${normal_id}.getpileupsummaries.table

        # Calculate contamination on paired samples
        echo $(date +"%F") $(date +"%T") "Calculating contamination of paired samples ..."
        gatk CalculateContamination \
            -I ${VAROUT}/${tumour_id}.getpileupsummaries.table \
            -matched ${VAROUT}/${normal_id}.getpileupsummaries.table \
            -O $VAROUT/${tumour_id}.contamination.table \
            -segments $VAROUT/${tumour_id}.segments.table

        
        # Call variant on paired tumour samples
        echo $(date +"%F") $(date +"%T") "Calling somatic variants on paired samples ...";
        gatk --java-options -Xmx8g Mutect2 \
            -I $bam_dir/${tumour_id}/${tumour_id}_recalibrated.bam \
            -I $bam_dir/${normal_id}/${normal_id}_recalibrated.bam \
            -normal $normal_id \
            -R $REFERENCE \
            -L $INTERVAL \
            --germline-resource $GERMLINE  \
            --panel-of-normals $PON \
            --f1r2-tar-gz $VAROUT/${tumour_id}.f1r2.tar.gz \
            --native-pair-hmm-threads 8 \
            --callable-depth 20 \
            -O $VAROUT/${tumour_id}_unfiltered.vcf.gz \
            -bamout $VAROUT/${tumour_id}_realigned.bam \
            >& ${VAROUT}/${tumour_id}.Mutect2.log

    else
        # Calculate contamination based on tumour samples only
        echo $(date +"%F") $(date +"%T") "Calculating contamination of tumour-only samples ..."
        gatk CalculateContamination \
            -I ${VAROUT}/${tumour_id}.getpileupsummaries.table \
            -O $VAROUT/${tumour_id}.contamination.table \
            -segments $VAROUT/${tumour_id}.segments.table
        
        echo $(date +"%F") $(date +"%T") "########Calling somatic variants on unpaired samples ...";

        gatk --java-options -Xmx8g Mutect2 \
            -I $bam_dir/${tumour_id}/${tumour_id}_recalibrated.bam \
            -R $REFERENCE \
            -L $INTERVAL \
            --germline-resource $GERMLINE  \
            --panel-of-normals $PON \
            --f1r2-tar-gz $VAROUT/${tumour_id}.f1r2.tar.gz \
            --callable-depth 20 \
            -O $VAROUT/${tumour_id}_unfiltered.vcf.gz \
            -bamout $VAROUT/${tumour_id}_realigned.bam \
            >& ${VAROUT}/${tumour_id}.Mutect2.log

    fi

    ##--af-of-alleles-not-in-resource 0.0000025

    # echo $(date +"%F") $(date +"%T") "Learning Read Orientation Model......"
    gatk LearnReadOrientationModel \
        -I ${vcf_dir}/${tumour_id}/${tumour_id}.f1r2.tar.gz \
        -O $VAROUT/${tumour_id}.read-orientation-model.tar.gz

    # Filter Mutect calls
    echo $(date +"%F") $(date +"%T") "Filtering Mutect calls ..."
    gatk --java-options -Xmx4g FilterMutectCalls \
        --variant $vcf_dir/${tumour_id}/${tumour_id}_unfiltered.vcf.gz \
        --stats $vcf_dir/${tumour_id}/${tumour_id}_unfiltered.vcf.gz.stats \
        --reference $REFERENCE \
        --ob-priors $VAROUT/${tumour_id}.read-orientation-model.tar.gz \
        --contamination-table $VAROUT/${tumour_id}.contamination.table \
        --tumor-segmentation $VAROUT/${tumour_id}.segments.table \
        --min-allele-fraction 0.01 \
        --unique-alt-read-count 1 \
        --output $VAROUT/${tumour_id}_filtered.vcf.gz
        # >& ${VAROUT}/${tumour_id}.FilterMutectCalls.log

    # Normalize reads
    echo $(date +"%F") $(date +"%T") "Normalizing Mutect calls ..."
    bcftools norm \
        $VAROUT/${tumour_id}_filtered.vcf.gz \
        -m-both -f $REFERENCE \
        -Oz -o $VAROUT/${tumour_id}_normalized.vcf.gz 

    echo $(date +"%F") $(date +"%T") "Normalizing filtered Mutect calls ..."
    
    bcftools view \
    -f PASS $VAROUT/${tumour_id}_normalized.vcf.gz \
    -o $VAROUT/${tumour_id}_normalized_filtered.vcf.gz

    # Filter repeatmasker and blacklist regions
    echo $(date +"%F") $(date +"%T") "Annotating repeatmasker regions ..."

    bcftools annotate \
        $VAROUT/${tumour_id}_normalized_filtered.vcf.gz \
        --header-lines ${work_dir}/vcf.rm.header \
        --annotations ${ref_dir}/RepeatMasker.bed.gz \
        --columns CHROM,FROM,TO,RepeatMasker \
        --output $VAROUT/${tumour_id}_repeatmasker.vcf.gz

    echo $(date +"%F") $(date +"%T") "Annotating blacklist regions ..."
    bcftools annotate \
        $VAROUT/${tumour_id}_repeatmasker.vcf.gz \
        --header-lines ${work_dir}/vcf.map.header \
        --annotations ${ref_dir}/blacklist.bed.gz \
        --columns CHROM,FROM,TO,EncodeDacMapability \
        --output-type z \
        --output $VAROUT/${tumour_id}_repeatmasker_blacklist.vcf.gz
    
    # Filter out variants in RepeatMasker or Mapability
    echo $(date +"%F") $(date +"%T") "Filtering RepeatMasker and blacklist regions ..."
    bcftools filter \
        $VAROUT/${tumour_id}_repeatmasker_blacklist.vcf.gz \
        -e 'INFO/RepeatMasker != "." || INFO/EncodeDacMapability != "."' \
        -Oz \
        -o $VAROUT/${tumour_id}.final.filtered.vcf.gz

    tabix $VAROUT/${tumour_id}.final.filtered.vcf.gz

    # Annotate variants by Funcotator
    echo $(date +"%F") $(date +"%T") "Annotating variants with Funcotator ..."
    gatk Funcotator \
        -R $REFERENCE \
        -V $VAROUT/${tumour_id}.final.filtered.vcf.gz \
        -O $VAROUT/${tumour_id}_annotated.maf.gz \
        -L $INTERVAL \
        --output-file-format MAF \
        --data-sources-path $ANNOTATION_FILE \
        --ref-version hg38 \
        --remove-filtered-variants true
        # >& ${VAROUT}/${tumour_id}.Funcotator.log &

    less -S $VAROUT/${tumour_id}_annotated.maf.gz | \
        grep -v "#" > $VAROUT/${tumour_id}_annotated.tsv

    #Annotate by Annovar
    echo $(date +"%F") $(date +"%T") "Annotating variants with Annovar ..."
    
    perl ${ref_dir}/annovar/table_annovar.pl $VAROUT/${tumour_id}.final.filtered.vcf.gz \
        ${ref_dir}/annovar/humandb/ \
        -buildver hg38 \
        -out $VAROUT/${tumour_id} \
        -remove \
        -protocol refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70 \
        -operation gx,r,f,f,f,f,f \
        -nastring . -polish -xreffile ${ref_dir}/annovar/example/gene_fullxref.txt \
        --otherinfo --vcfinput \
        >& ${VAROUT}/${tumour_id}.annovar.log

    less -S $VAROUT/${tumour_id}.hg38_multianno.txt | \
        awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "AD", "AF", "DP"}; NR >1 {split($NF, a, ":"); $(NF+1)=a[2]; $(NF+1)=a[3]; $(NF+1)=a[4]; print}' > \
        $VAROUT/${tumour_id}_annovar.txt

done
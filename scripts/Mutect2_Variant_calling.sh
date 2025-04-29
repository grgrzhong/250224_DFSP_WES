#!/bin/bash

# Use source to activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate varcall

work_dir="/home/zhonggr/projects/250224_DFSP_WES/data/mfs"
ref_dir="/home/zhonggr/projects/250224_DFSP_WES/data/reference"

REFERENCE=${ref_dir}/Gencode/gencode.hg38.v36.primary_assembly.fa
GERMLINE=${ref_dir}/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
ANNOTATION_FILE=${ref_dir}/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/
INTERVAL=${ref_dir}/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed
PON=${ref_dir}/pon_dfsp/pon.vcf.gz

BAM=${work_dir}/bam

CASELIST=$(ls $BAM | cut -d '-' -f1,2 | uniq)

# for case in $CASELIST; do 
#     echo $case

#     TUMOURS=$(ls $BAM | grep ${case}-T)
#     NORMAL=${case}-N

#     for TUMOUR in $TUMOURS; do

#         echo $TUMOUR
#         echo $NORMAL

#         mkdir -p ${work_dir}/mutect2/${TUMOUR}
#         VAROUT=${work_dir}/mutect2/${TUMOUR}

#         # Get Pileup Summaries
#         echo $(date +"%F") $(date +"%T") "#####Getting Pileup Summaries of tumour sample......"
#         gatk GetPileupSummaries \
#             -I $BAM/$TUMOUR/${TUMOUR}_recalibrated.bam \
#             -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
#             -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
#             -O ${VAROUT}/${TUMOUR}.getpileupsummaries.table

#         #Check if presence of paired normal samples
#         if [ -d ${BAM}/${NORMAL} ]; then

#             echo $(date +"%F") $(date +"%T") "#####Getting Pileup Summaries of normal sample......"
#             gatk GetPileupSummaries \
#                 -I $BAM/$NORMAL/${NORMAL}_recalibrated.bam \
#                 -V ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
#                 -L ${ref_dir}/GetPileupSummary/small_exac_common_3.hg38.vcf.gz \
#                 -O ${VAROUT}/${NORMAL}.getpileupsummaries.table

#             # Calculate contamination on paired samples
#             echo $(date +"%F") $(date +"%T") "#######Calculating contamination of paired samples......."
#             gatk CalculateContamination \
#                 -I ${VAROUT}/${TUMOUR}.getpileupsummaries.table \
#                 -matched ${VAROUT}/${NORMAL}.getpileupsummaries.table \
#                 -O $VAROUT/${TUMOUR}.contamination.table \
#                 -segments $VAROUT/${TUMOUR}.segments.table

#             # Call variant on paired tumour samples
#             echo $(date +"%F") $(date +"%T") "########Calling somatic variants on paired samples........."
#             gatk --java-options -Xmx128g Mutect2 \
#                 -I $BAM/${TUMOUR}/${TUMOUR}_recalibrated.bam \
#                 -I $BAM/${NORMAL}/${NORMAL}_recalibrated.bam \
#                 -normal $NORMAL \
#                 -R $REFERENCE \
#                 -L $INTERVAL \
#                 --germline-resource $GERMLINE  \
#                 --panel-of-normals $PON \
#                 --f1r2-tar-gz $VAROUT/${TUMOUR}.f1r2.tar.gz \
#                 --native-pair-hmm-threads 8 \
#                 --callable-depth 20 \
#                 -O $VAROUT/${TUMOUR}_unfiltered.vcf.gz \
#                 -bamout $VAROUT/${TUMOUR}_realigned.bam \
#                 >& ${VAROUT}/${TUMOUR}.Mutect2Call.log

#         else
#             # Calculate contamination based on tumour samples only
#             echo $(date +"%F") $(date +"%T") "#######Calculating contamination of tumour-only samples...."
#             gatk CalculateContamination \
#                 -I ${VAROUT}/${TUMOUR}.getpileupsummaries.table \
#                 -O $VAROUT/${TUMOUR}.contamination.table \
#                 -segments $VAROUT/${TUMOUR}.segments.table

#             # Call variant on unpaired tumour samples
#             echo $(date +"%F") $(date +"%T") "########Calling somatic variants on unpaired samples......."
#             gatk --java-options -Xmx128g Mutect2 \
#                 -I $BAM/${TUMOUR}/${TUMOUR}_recalibrated.bam \
#                 -R $REFERENCE \
#                 -L $INTERVAL \
#                 --germline-resource $GERMLINE  \
#                 --panel-of-normals $PON \
#                 --f1r2-tar-gz $VAROUT/${TUMOUR}.f1r2.tar.gz \
#                 --callable-depth 20 \
#                 -O $VAROUT/${TUMOUR}_unfiltered.vcf.gz \
#                 -bamout $VAROUT/${TUMOUR}_realigned.bam \
#                 >& ${VAROUT}/${TUMOUR}.Mutect2Call.log &
#         fi
#         ##--af-of-alleles-not-in-resource 0.0000025

#         echo $(date +"%F") $(date +"%T") "##########Learning Read Orientation Model......"
#         gatk LearnReadOrientationModel \
#             -I $VAROUT/${TUMOUR}.f1r2.tar.gz \
#             -O $VAROUT/${TUMOUR}.read-orientation-model.tar.gz

#         #Filter Mutect calls
#         echo $(date +"%F") $(date +"%T") "##########Filtering Mutect calls.............."
#         gatk --java-options -Xmx64g FilterMutectCalls \
#             -V $VAROUT/${TUMOUR}_unfiltered.vcf.gz \
#             -R $REFERENCE \
#             --ob-priors $VAROUT/${TUMOUR}.read-orientation-model.tar.gz \
#             --contamination-table $VAROUT/${TUMOUR}.contamination.table \
#             --tumor-segmentation $VAROUT/${TUMOUR}.segments.table \
#             --min-allele-fraction 0.01 \
#             --unique-alt-read-count 1 \
#             --stats $VAROUT/${TUMOUR}_unfiltered.vcf.gz.stats \
#             -O $VAROUT/${TUMOUR}_filtered.vcf.gz \
#             >& ${VAROUT}/${TUMOUR}.Mutect2Filter.log &

#         #Normalize reads
#         bcftools norm -m-both -f $REFERENCE -Oz -o $VAROUT/${TUMOUR}_normalized.vcf.gz $VAROUT/${TUMOUR}_filtered.vcf.gz

#         bcftools view -f PASS $VAROUT/${TUMOUR}_normalized.vcf.gz -o $VAROUT/${TUMOUR}_normalized_filtered.vcf.gz

#         tabix $VAROUT/${TUMOUR}_normalized_filtered.vcf.gz
#         rm $VAROUT/${TUMOUR}_normalized.vcf.gz

#         #Annotate variants by Funcotator
#         echo $(date +"%F") $(date +"%T") "###########Annotating variants.................."
#         gatk Funcotator \
#             -R $REFERENCE \
#             -V $VAROUT/${TUMOUR}_normalized_filtered.vcf.gz \
#             -O $VAROUT/${TUMOUR}_annotated.maf.gz \
#             -L $INTERVAL \
#             --output-file-format MAF \
#             --data-sources-path $ANNOTATION_FILE \
#             --ref-version hg38 \
#             --remove-filtered-variants true \
#             >& ${VAROUT}/${TUMOUR}.Funcotator.log &

#         less -S $VAROUT/${TUMOUR}_annotated.maf.gz | grep -v "#" > $VAROUT/${TUMOUR}_annotated.tsv

#         #Annotate by Annovar
#         perl ${ref_dir}/annovar/table_annovar.pl $VAROUT/${TUMOUR}_normalized_filtered.vcf.gz \
#             ${ref_dir}/annovar/humandb/ \
#             -buildver hg38 -out $VAROUT/${TUMOUR} -remove \
#             -protocol refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70 \
#             -operation gx,r,f,f,f,f,f \
#             -nastring . -polish -xreffile /media/maximus/Data/Software/annovar/example/gene_fullxref.txt \
#             --otherinfo --vcfinput \
#             >& ${VAROUT}/${TUMOUR}.Annovar.log &

#         less -S $VAROUT/${TUMOUR}.hg38_multianno.txt | \
#             awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "AD", "AF", "DP"}; NR >1 {split($NF, a, ":"); $(NF+1)=a[2]; $(NF+1)=a[3]; $(NF+1)=a[4]; print}' > \
#             $VAROUT/${TUMOUR}_annovar.txt

#     done
# done

for case in $CASELIST; do 
    echo $case

    TUMOURS=$(ls $BAM | grep ${case}-T)
    NORMAL=${case}-N

    for TUMOUR in $TUMOURS; do

        echo $TUMOUR
        echo $NORMAL

        mkdir -p ${work_dir}/mutect2/${TUMOUR}
        VAROUT=${work_dir}/mutect2/${TUMOUR}

        # #Filter Mutect calls
        # echo $(date +"%F") $(date +"%T") "##########Filtering Mutect calls.............."
        # gatk --java-options -Xmx8g FilterMutectCalls \
        #     -V $VAROUT/${TUMOUR}_unfiltered.vcf.gz \
        #     -R $REFERENCE \
        #     --ob-priors $VAROUT/${TUMOUR}.read-orientation-model.tar.gz \
        #     --contamination-table $VAROUT/${TUMOUR}.contamination.table \
        #     --tumor-segmentation $VAROUT/${TUMOUR}.segments.table \
        #     --min-allele-fraction 0.01 \
        #     --unique-alt-read-count 1 \
        #     --stats $VAROUT/${TUMOUR}_unfiltered.vcf.gz.stats \
        #     -O $VAROUT/${TUMOUR}_filtered.vcf.gz \
        #     >& ${VAROUT}/${TUMOUR}.Mutect2Filter.log &

        # # Normalize reads
        # bcftools norm -m-both -f $REFERENCE -Oz -o $VAROUT/${TUMOUR}_normalized.vcf.gz $VAROUT/${TUMOUR}_filtered.vcf.gz

        # bcftools view -f PASS $VAROUT/${TUMOUR}_normalized.vcf.gz -o $VAROUT/${TUMOUR}_normalized_filtered.vcf.gz

        # tabix $VAROUT/${TUMOUR}_normalized_filtered.vcf.gz
        # rm $VAROUT/${TUMOUR}_normalized.vcf.gz

        # #Annotate variants by Funcotator
        # echo $(date +"%F") $(date +"%T") "###########Annotating variants.................."
        # gatk Funcotator \
        #     -R $REFERENCE \
        #     -V $VAROUT/${TUMOUR}_normalized_filtered.vcf.gz \
        #     -O $VAROUT/${TUMOUR}_annotated.maf.gz \
        #     -L $INTERVAL \
        #     --output-file-format MAF \
        #     --data-sources-path $ANNOTATION_FILE \
        #     --ref-version hg38 \
        #     --remove-filtered-variants true \
        #     >& ${VAROUT}/${TUMOUR}.Funcotator.log &

        # less -S $VAROUT/${TUMOUR}_annotated.maf.gz | grep -v "#" > $VAROUT/${TUMOUR}_annotated.tsv

        #Annotate by Annovar
        perl ${ref_dir}/annovar/table_annovar.pl $VAROUT/${TUMOUR}_normalized_filtered.vcf.gz \
            ${ref_dir}/annovar/humandb/ \
            -buildver hg38 -out $VAROUT/${TUMOUR} -remove \
            -protocol refGene,cytoBand,dbnsfp33a,gnomad_exome,avsnp150,clinvar_20221231,cosmic70 \
            -operation gx,r,f,f,f,f,f \
            -nastring . -polish -xreffile ${ref_dir}/annovar/example/gene_fullxref.txt \
            --otherinfo --vcfinput \
            >& ${VAROUT}/${TUMOUR}.Annovar.log &

        less -S $VAROUT/${TUMOUR}.hg38_multianno.txt | awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "AD", "AF", "DP"}; NR >1 {split($NF, a, ":"); $(NF+1)=a[2]; $(NF+1)=a[3]; $(NF+1)=a[4]; print}' > $VAROUT/${TUMOUR}_annovar.txt

    done
done
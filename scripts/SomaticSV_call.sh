#!/bin/bash

SERVER=/run/user/1000/gvfs/smb-share:server=maximus-nas.local,share=genomics
REFERENCE=/media/Data/Reference/Gencode/gencode.hg38.v36.primary_assembly.fa
GERMLINE=/media/Data/Reference/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
ANNOTATION_FILE=/media/Data/Reference/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/
REGIONS=/media/Data/Reference/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed.gz
BAM=/media/Data/WES/BAM
EXCLUDE=/media/Data/Software/delly/excludeTemplates/human.hg38.excl.tsv

CASELIST=$(ls $BAM | cut -d '-' -f1,2 | uniq)

export PATH=/media/Data/delly/src:$PATH
export ANNOTSV=~/Software/AnnotSV

for case in $CASELIST; do 
echo $case
TUMOURS=$(ls $BAM | grep ${case}-T)
NORMAL=${case}-N;

for TUMOUR in $TUMOURS; do

echo $TUMOUR
echo $NORMAL

mkdir -p ${BAM}/../SV-call/${TUMOUR}
OUTPUT=${BAM}/../SV-call/${TUMOUR}

#Call somaticSV using Delly
mkdir -p ${BAM}/../SV-call/${TUMOUR}/Delly
OUTPUTDELLY=${BAM}/../SV-call/${TUMOUR}/Delly

echo $(date +"%F") $(date +"%T") "##########Running Delly......";
/media/Data/Software/delly_v1.2.6_linux_x86_64bit call \
-x $EXCLUDE \
-o ${OUTPUTDELLY}/${TUMOUR}.bcf \
-g $REFERENCE \
${BAM}/${TUMOUR}/${TUMOUR}_recalibrated.bam \
${BAM}/${NORMAL}/${NORMAL}_recalibrated.bam

echo -e "${TUMOUR}\ttumor\n${NORMAL}\tcontrol" > ${OUTPUTDELLY}/samples.tsv

/media/Data/Software/delly_v1.2.6_linux_x86_64bit filter \
-f somatic \
-o ${OUTPUTDELLY}/${TUMOUR}.somatic.filter.bcf \
-s ${OUTPUTDELLY}/samples.tsv \
${OUTPUTDELLY}/${TUMOUR}.bcf

bcftools view ${OUTPUTDELLY}/${TUMOUR}.somatic.filter.bcf -Oz > ${OUTPUTDELLY}/${TUMOUR}.somaticSV.delly.vcf.gz;

tabix ${OUTPUTDELLY}/${TUMOUR}.somaticSV.delly.vcf.gz;

#Call somaticSV using Manta
mkdir -p ${BAM}/../SV-call/${TUMOUR}/Manta
OUTPUTMANTA=${BAM}/../SV-call/${TUMOUR}/Manta

echo $(date +"%F") $(date +"%T") "##########Running Manta......";
python2 /media/Data/Software/manta-1.6.0.centos6_x86_64/bin/configManta.py \
--normalBam=$BAM/${NORMAL}/${NORMAL}_recalibrated.bam \
--tumourBam=$BAM/${TUMOUR}/${TUMOUR}_recalibrated.bam \
--exome \
--referenceFasta=$REFERENCE \
--runDir=$OUTPUTMANTA \
--callRegions=$REGIONS

echo $(date +"%F") $(date +"%T") "##########Running Manta2......";
python2 $OUTPUTMANTA/runWorkflow.py -j 8

echo $(date +"%F") $(date +"%T") "##########Running Manta Conversion......";
python2 /media/Data/Software/manta-1.6.0.centos6_x86_64/libexec/convertInversion.py /usr/bin/samtools $REFERENCE $OUTPUTMANTA/results/variants/somaticSV.vcf.gz
tabix -f --preset vcf $OUTPUTMANTA/results/variants/somaticSV.vcf.gz

mv $OUTPUTMANTA/results/variants/candidateSmallIndels.vcf.gz \
  $OUTPUTMANTA/${TUMOUR}.candidateSmallIndels.vcf.gz
mv $OUTPUTMANTA/results/variants/candidateSmallIndels.vcf.gz.tbi \
  $OUTPUTMANTA/${TUMOUR}.candidateSmallIndels.vcf.gz.tbi
mv $OUTPUTMANTA/results/variants/candidateSV.vcf.gz \
  $OUTPUTMANTA/${TUMOUR}.candidateSV.vcf.gz
mv $OUTPUTMANTA/results/variants/candidateSV.vcf.gz.tbi \
  $OUTPUTMANTA/${TUMOUR}.candidateSV.vcf.gz.tbi
mv $OUTPUTMANTA/results/variants/diploidSV.vcf.gz \
  $OUTPUTMANTA/${TUMOUR}.diploidSV.vcf.gz
mv $OUTPUTMANTA/results/variants/diploidSV.vcf.gz.tbi \
  $OUTPUTMANTA/${TUMOUR}.diploidSV.vcf.gz.tbi
mv $OUTPUTMANTA/results/variants/somaticSV.vcf.gz \
  $OUTPUTMANTA/${TUMOUR}.somaticSV.manta.vcf.gz
mv $OUTPUTMANTA/results/variants/somaticSV.vcf.gz.tbi \
  $OUTPUTMANTA/${TUMOUR}.somaticSV.manta.vcf.gz.tbi
  
#Merge outputs from Delly and Manta
echo $(date +"%F") $(date +"%T") "##########Merge output ......";
bcftools view \
--samples ${TUMOUR},${NORMAL} \
--output-type z \
--output-file ${OUTPUTMANTA}/${TUMOUR}.manta.swap.vcf.gz \
${OUTPUTMANTA}/${TUMOUR}.somaticSV.manta.vcf.gz 
    
tabix --preset vcf ${OUTPUTMANTA}/${TUMOUR}.manta.swap.vcf.gz

bcftools concat \
--allow-overlaps \
--output-type z \
--output ${OUTPUT}/${TUMOUR}.delly.manta.unfiltered.vcf.gz \
${OUTPUTDELLY}/${TUMOUR}.somaticSV.delly.vcf.gz \
${OUTPUTMANTA}/${TUMOUR}.manta.swap.vcf.gz

tabix --preset vcf ${OUTPUT}/${TUMOUR}.delly.manta.unfiltered.vcf.gz

echo $(date +"%F") $(date +"%T") "########Filtering somatic SV......";
bcftools filter \
--include 'FILTER="PASS"' \
${OUTPUT}/${TUMOUR}.delly.manta.unfiltered.vcf.gz | \
bcftools sort \
--output-type z \
--output-file ${OUTPUT}/${TUMOUR}.delly.manta.vcf.gz 

tabix --preset vcf ${OUTPUT}/${TUMOUR}.delly.manta.vcf.gz

echo $(date +"%F") $(date +"%T") "########Annotating somatic SV......"
$ANNOTSV/bin/AnnotSV \
-SVinputFile ${OUTPUT}/${TUMOUR}.delly.manta.vcf.gz \
-outputFile ${OUTPUT}/${TUMOUR}.somaticSV.annotated.tsv \
-genomeBuild GRCh38 \
-annotationMode both \
-SVminSize 50 \
>& ${OUTPUT}/${TUMOUR}.AnnotSV.log &

done
done

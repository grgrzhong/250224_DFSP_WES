#!/bin/bash

SERVER=/run/user/1000/gvfs/smb-share:server=maximus-nas.local,share=genomics
REFERENCE=/media/maximus/Data/Reference/Gencode/gencode.hg38.v36.primary_assembly.fa
BAITINTERVAL=/media/maximus/Data/Reference/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-probes-hg38.interval_list
TARGETINTERVAL=/media/maximus/Data/Reference/Exome/xgen-exome-hyb-panel-v2/hg38/xgen-exome-hyb-panel-v2-targets-hg38.interval_list
INTERVAL=/media/maximus/Data/Reference/Exome/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged/xgen-exome-hyb-panel-v2-hg38_200bp_sorted_merged.bed
GERMLINE=/media/maximus/Data/Reference/Population_database/somatic-hg38_af-only-gnomad.hg38.vcf.gz
ANNOTATION_FILE=/media/maximus/Data/Reference/Funocator_Datasource/funcotator_dataSources.v1.7.20200521s/
DBSNP=/media/maximus/Data/Reference/Population_database/dbSNP.vcf.gz
INPUT=/media/maximus/Data/WES/1.Fastq-trimmed

export PATH=/media/maximus/Data/Software/gatk-4.4.0.0/:$PATH
#alias fgbio="java -Xmx64g -jar /media/Data/Software/fgbio-2.2.1.jar"

for file in $(ls ${INPUT}/*.fastq.gz | cut -d "_" -f 1 | uniq); do 
FILENAME=$(basename $file); 
echo $FILENAME;
mkdir -p ${INPUT}/../BAM/${FILENAME};
OUTPUT=${INPUT}/../BAM/${FILENAME};
LOGFILE=${OUTPUT}/${FILENAME}_log;

{
#Alignment by bwa
echo $(date +"%F") $(date +"%T") "Aligning '$FILENAME' to Reference genome ........";
bwa mem -M -t 16 \
-R "@RG\tID:$FILENAME\tLB:XGenV2\tPL:ILLUMINA\tPM:NOVASEQ\tSM:$FILENAME\tPU:NA" \
$REFERENCE ${file}_trimmed_1.fastq.gz ${file}_trimmed_2.fastq.gz | samtools view -Sb - > ${OUTPUT}/${FILENAME}.bam

#Extract and tag umi information
echo $(date +"%F") $(date +"%T") "Extract and tag UMI......................";
python UMI.py -I ${OUTPUT}/${FILENAME}.bam -O ${OUTPUT}/${FILENAME}_umi.bam

##Sort by coordinate
echo $(date +"%F") $(date +"%T") "Sort by coordinate......................";
samtools sort ${OUTPUT}/${FILENAME}_umi.bam -o ${OUTPUT}/${FILENAME}_sorted.bam;

echo $(date +"%F") $(date +"%T") "Mark duplicates......................";
gatk --java-options -Xmx96g MarkDuplicates \
-I ${OUTPUT}/${FILENAME}_sorted.bam \
-M ${OUTPUT}/${FILENAME}_metrics.txt \
-O ${OUTPUT}/${FILENAME}_marked.bam \
--BARCODE_TAG "RX";

echo $(date +"%F") $(date +"%T") "Making index......................";
samtools index ${OUTPUT}/${FILENAME}_marked.bam;

echo $(date +"%F") $(date +"%T") "Base recalibration......................";
gatk BaseRecalibrator \
-I ${OUTPUT}/${FILENAME}_marked.bam \
-R $REFERENCE \
-L $INTERVAL \
-O ${OUTPUT}/${FILENAME}_recal_table.table \
--known-sites $DBSNP ;

gatk ApplyBQSR \
-I ${OUTPUT}/${FILENAME}_marked.bam \
-O ${OUTPUT}/${FILENAME}_recalibrated.bam \
-L $INTERVAL \
-bqsr ${OUTPUT}/${FILENAME}_recal_table.table \
--create-output-bam-md5;

echo $(date +"%F") $(date +"%T") "Collect metrics....................";
gatk CollectHsMetrics \
-I ${OUTPUT}/${FILENAME}_recalibrated.bam \
-O ${OUTPUT}/${FILENAME}_hs_metrics.txt \
-R $REFERENCE \
-BI $BAITINTERVAL \
-TI $TARGETINTERVAL

#gatk CollectWgsMetrics \
#-I ${OUTPUT}/${FILENAME}_recalibrated.bam \
#-O ${OUTPUT}/${FILENAME}_wgs_metrics.txt \
#-R $REFERENCE;

# Generate alignment stats
bamtools stats -in ${OUTPUT}/${FILENAME}_recalibrated.bam > ${OUTPUT}/${FILENAME}_aln_stat.txt;

rm ${OUTPUT}/${FILENAME}.bam ${OUTPUT}/${FILENAME}_umi.bam ${OUTPUT}/${FILENAME}_sorted.bam ${OUTPUT}/${FILENAME}_marked.bam ${OUTPUT}/${FILENAME}_marked.bam.bai;

} &> $LOGFILE;

done

#!/bin/bash

INPUT=/media/Data/WES/Input
OUTPUT=/media/Data/WES/1.Fastq-trimmed

for file in $(ls ${INPUT}/*.fastq.gz | cut -d "_" -f 1 | uniq); do 
FILENAME=`echo $(basename $file)`; 
echo $FILENAME;

/media/Data/Software/fastp \
-i ${file}_1.fastq.gz \
-I ${file}_2.fastq.gz \
-o ${OUTPUT}/${FILENAME}_trimmed_1.fastq.gz \
-O ${OUTPUT}/${FILENAME}_trimmed_2.fastq.gz \
--detect_adapter_for_pe \
--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--trim_poly_g \
--trim_poly_x \
--umi \
--umi_loc per_read \
--umi_len 8 \
-j ${OUTPUT}/${FILENAME}.json \
-h ${OUTPUT}/${FILENAME}.html \
-w 8

#--dont_eval_duplication

done

mkdir -p ${INPUT}/../FastQC-trimmed

/media/Data/Software/FastQC/fastqc \
${OUTPUT}/*.fastq.gz \
-o ${INPUT}/../FastQC-trimmed \
--memory 8192 -t 8 



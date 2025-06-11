#!/bin/bash

#############################################################################
conda activate varcall

# remove the final vcf header lines that are not needed for PCGR annotation
mutect2_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/mutect2"

sample_list=$(ls "${mutect2_dir}")
echo "Sample list: ${sample_list}"

for sample in ${sample_list}; do
    
    echo "Processing sample: ${sample}"
    
    # Define input and output VCF paths
    input_vcf="${mutect2_dir}/${sample}/${sample}.final.vcf.gz"
    output_vcf="${mutect2_dir}/${sample}/${sample}.final.cleaned.vcf.gz"
    
    # Check if input VCF exists
    if [ ! -f "${input_vcf}" ]; then
        echo "Input VCF not found: ${input_vcf}"
        continue
    fi
    
    # Remove specific INFO fields from the VCF header
    # bcftools annotate \
    #     --remove INFO/RepeatMasker,INFO/EncodeDacMapability \
    #     --output-type z \
    #     --output "${output_vcf}" \
    #     "${input_vcf}"
    
    # ## remove the final bcf and rename the cleaned final vcf
    # rm -f "${mutect2_dir}/${sample}/${sample}.final.bcf.gz"
    # rm -f "${mutect2_dir}/${sample}/${sample}.final.bcf.gz.tbi"

    # mv "${output_vcf}" "${mutect2_dir}/${sample}/${sample}.final.vcf.gz"

    # Index the cleaned VCF file
    rm -rf "${mutect2_dir}/${sample}/${sample}.final.vcf.gz.tbi"
    tabix -p vcf "${mutect2_dir}/${sample}/${sample}.final.vcf.gz"
    
    # echo "Cleaned VCF created: ${output_vcf}"
done

vcf_old=/home/zhonggr/projects/250224_DFSP_WES/data/wes/mutect2_old/DFSP-001-T/DFSP-001-T.final.vcf.gz
vcf_new=/home/zhonggr/projects/250224_DFSP_WES/data/wes/mutect2/DFSP-001-T/DFSP-001-T.final.vcf.gz
vcf_cleaned=/home/zhonggr/projects/250224_DFSP_WES/data/wes/mutect2/DFSP-001-T/DFSP-001-T.final.cleaned.vcf.gz


bcftools view -h "${vcf_old}" | grep "##INFO"
bcftools view -h "${vcf_new}" | grep "##INFO"

# bcftools annotate \
#     --remove INFO/RepeatMasker,INFO/EncodeDacMapability \
#     --output-type z \
#     --output "${vcf_cleaned}" \
#     "${vcf_new}"

# bcftools view -h "${vcf_cleaned}" | grep "##INFO"
# bcftools view --no-header "${vcf_cleaned}" | wc -l
# bcftools view --no-header "${vcf_new}" | wc -l
# bcftools view --no-header "${vcf_old}" | wc -l
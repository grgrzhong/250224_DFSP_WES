#!/bin/bash
# filepath: /home/zhonggr/projects/250224_DFSP_WES/scripts/normal_vaf.sh

cur_dir=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T
tumour_id=DFSP-001-T
bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
normal_id="DFSP-001-N"

echo "=== DEBUGGING VARIANT KEYS ==="
echo "First VCF variant:"
bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${cur_dir}/${tumour_id}.final.vcf.gz | head -1
echo "First Annovar variant:"
awk 'NR==2 {print $1":"$2":"$4":"$5}' ${cur_dir}/${tumour_id}.hg38_multianno.txt

echo "=== VCF SAMPLE ORDER ==="
bcftools query -l ${cur_dir}/${tumour_id}.final.vcf.gz

echo "=== DEBUGGING FORMAT FIELD ==="
echo "Sample data from first variant:"
bcftools query -f '%CHROM\t%POS\t[%SAMPLE=%GT:%AD:%AF:%DP]\n' ${cur_dir}/${tumour_id}.final.vcf.gz | head -1

# Check if this is a paired sample
if [ -d ${bam_dir}/${normal_id} ]; then
    echo ""
    echo "Processing paired sample: ${tumour_id} with normal ${normal_id}"
    
    # For paired samples: the last field in ANNOVAR contains the VCF format data
    # Structure: GT:AD:AF:DP:F1R2:F2R1:FAD:SB\t[normal_data]\t[tumor_data]\t[additional_fields]
    awk 'BEGIN {FS=OFS="\t"}
    NR==1 {
        # Add columns for extracted data
        print $0, "AD", "AF", "DP", "TUMOR_AD", "TUMOR_AF", "TUMOR_DP", "NORMAL_AD", "NORMAL_AF", "NORMAL_DP", "SAMPLE_TYPE"
    }
    NR > 1 {
        # Split the last field by tabs first to separate format from sample data
        n = split($NF, vcf_parts, "\t")
        
        if (n >= 3) {
            # vcf_parts[1] = FORMAT string (GT:AD:AF:DP:F1R2:F2R1:FAD:SB)
            # vcf_parts[2] = Normal sample data (0/0:38,0:0.045:38:1,0:6,0:20,0:19,19,0,0)
            # vcf_parts[3] = Tumor sample data (0/1:52,2:0.065:54:6,0:4,1:28,1:26,26,1,1)
            
            # Parse normal sample data (first sample)
            split(vcf_parts[2], normal_fields, ":")
            normal_gt = normal_fields[1]
            normal_ad = normal_fields[2]
            normal_af = normal_fields[3]
            normal_dp = normal_fields[4]
            
            # Parse tumor sample data (second sample)
            split(vcf_parts[3], tumor_fields, ":")
            tumor_gt = tumor_fields[1]
            tumor_ad = tumor_fields[2]
            tumor_af = tumor_fields[3]
            tumor_dp = tumor_fields[4]
            
        } else {
            # Fallback: if the format is different, try to parse as concatenated string
            # This happens when the VCF data is stored differently in ANNOVAR
            full_string = $NF
            
            # Look for pattern like: GT:AD:AF:DP:...SAMPLE1_DATA...SAMPLE2_DATA
            # Try to identify where normal and tumor data start
            if (match(full_string, /0\/0:[0-9,]+:[0-9.]+:[0-9]+/)) {
                # Found normal sample pattern (usually 0/0)
                normal_match = substr(full_string, RSTART, RLENGTH)
                split(normal_match, normal_fields, ":")
                normal_ad = normal_fields[2]
                normal_af = normal_fields[3]
                normal_dp = normal_fields[4]
                
                # Look for tumor pattern after normal
                remaining = substr(full_string, RSTART + RLENGTH)
                if (match(remaining, /0\/1:[0-9,]+:[0-9.]+:[0-9]+/)) {
                    tumor_match = substr(remaining, RSTART, RLENGTH)
                    split(tumor_match, tumor_fields, ":")
                    tumor_ad = tumor_fields[2]
                    tumor_af = tumor_fields[3]
                    tumor_dp = tumor_fields[4]
                } else {
                    tumor_ad = tumor_af = tumor_dp = "NA"
                }
            } else {
                normal_ad = normal_af = normal_dp = "NA"
                tumor_ad = tumor_af = tumor_dp = "NA"
            }
        }
        
        # Print original data + extracted values
        print $0, tumor_ad, tumor_af, tumor_dp, tumor_ad, tumor_af, tumor_dp, normal_ad, normal_af, normal_dp, "PAIRED"
    }' ${cur_dir}/${tumour_id}.hg38_multianno.txt > ${cur_dir}/${tumour_id}.annovar_with_normal.txt

else
    echo "Processing tumor-only sample: ${tumour_id}"
    
    # For tumor-only samples
    awk 'BEGIN {FS=OFS="\t"}
    NR==1 {
        print $0, "AD", "AF", "DP", "TUMOR_AD", "TUMOR_AF", "TUMOR_DP", "NORMAL_AD", "NORMAL_AF", "NORMAL_DP", "SAMPLE_TYPE"
    }
    NR > 1 {
        # Extract tumor data only
        split($NF, format_fields, ":")
        tumor_ad = format_fields[2]
        tumor_af = format_fields[3]
        tumor_dp = format_fields[4]
        
        print $0, tumor_ad, tumor_af, tumor_dp, tumor_ad, tumor_af, tumor_dp, "NA", "NA", "NA", "TUMOR_ONLY"
    }' ${cur_dir}/${tumour_id}.hg38_multianno.txt > ${cur_dir}/${tumour_id}.annovar_with_normal.txt

fi

# Enhanced verification using actual VCF data
echo ""
echo "=== VERIFICATION WITH VCF COMPARISON ==="
echo "Comparing first 3 variants between VCF and processed ANNOVAR:"

# Get VCF data for comparison
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE]\t[%AD]\t[%AF]\t[%DP]\n' ${cur_dir}/${tumour_id}.final.vcf.gz | head -6 > /tmp/vcf_data.tmp

echo "VCF sample data (first 3 variants):"
cat /tmp/vcf_data.tmp

echo ""
echo "ANNOVAR processed data (first 3 variants):"
awk 'NR>1 && NR<=4 {
    n = NF
    tumor_ad = $(n-7)
    tumor_af = $(n-6)
    tumor_dp = $(n-5)
    normal_ad = $(n-4)
    normal_af = $(n-3)
    normal_dp = $(n-2)
    
    print $1 "\t" $2 "\t" $4 "\t" $5 "\tNORMAL\t" normal_ad "\t" normal_af "\t" normal_dp
    print $1 "\t" $2 "\t" $4 "\t" $5 "\tTUMOR\t" tumor_ad "\t" tumor_af "\t" tumor_dp
    print ""
}' ${cur_dir}/${tumour_id}.annovar_with_normal.txt

# Summary statistics
echo "=== SUMMARY STATISTICS ==="
awk 'NR>1 {
    n = NF
    tumor_af = $(n-6)
    normal_af = $(n-3)
    sample_type = $n
    
    total++
    if (sample_type == "PAIRED") {
        paired++
        if (normal_af != "NA" && normal_af + 0 <= 0.01) low_normal_vaf++
        if (normal_af != "NA" && normal_af + 0 == 0) zero_normal_vaf++
        if (normal_af != "NA" && normal_af + 0 > 0.05) high_normal_vaf++
    }
    
    if (tumor_af != "NA" && tumor_af + 0 > 0.05) high_tumor_vaf++
    
} END {
    print "Total variants: " total
    print "Paired samples: " paired
    print "High tumor VAF (>5%): " high_tumor_vaf
    if (paired > 0) {
        print "Normal VAF = 0: " zero_normal_vaf
        print "Normal VAF <= 1%: " low_normal_vaf
        print "Normal VAF > 5%: " high_normal_vaf
        print "Percentage with low normal VAF (<=1%): " (low_normal_vaf/paired)*100 "%"
        print "Potential germline variants (normal VAF > 5%): " high_normal_vaf
    }
}' ${cur_dir}/${tumour_id}.annovar_with_normal.txt

# Clean up
rm -f /tmp/vcf_data.tmp

echo ""
echo "Output file created: ${cur_dir}/${tumour_id}.annovar_with_normal.txt"
echo "Use this file for downstream analysis with proper tumor/normal VAF separation."
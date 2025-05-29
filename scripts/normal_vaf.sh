
vcf=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-037-T/DFSP-037-T.final.vcf.gz 

vcf=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T/DFSP-001-T.final.vcf.gz

bcftools view -h ${vcf} | grep "^#CHROM"

bcftools view -h ${vcf} | grep "^##FORMAT=<ID=AF"

# Extract sample data to see what's actually in the normal vs tumor columns
echo "Sample order:"
bcftools query -l ${vcf}

echo -e "\nFirst 5 variants with detailed sample information:"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE]\t[%GT]\t[%AD]\t[%AF]\t[%DP]\n' ${vcf} | head -10

echo -e "\nFormatted view of first 3 variants:"
bcftools query -f '%CHROM:%POS %REF>%ALT\n[%SAMPLE: GT=%GT AD=%AD AF=%AF DP=%DP]\n' ${vcf} | head -12

# Extract normal sample AF values specifically
echo -e "\nNormal sample (DFSP-001-N) AF values for first 10 variants:"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AF]\n' -s DFSP-001-N ${vcf} | head -10

# Check how many variants have normal AF <= 0.01 (1%)
echo -e "\nAnalyzing normal sample VAF distribution:"
bcftools query -f '[%AF]\n' -s DFSP-001-N ${vcf} | \
awk '{
    total++
    if ($1 + 0 <= 0.01) low_vaf++
    if ($1 + 0 == 0) zero_vaf++
} 
END {
    print "Total variants: " total
    print "Normal VAF = 0: " zero_vaf
    print "Normal VAF <= 0.01 (1%): " low_vaf
    print "Percentage with VAF <= 1%: " (low_vaf/total)*100 "%"
}'

###############################################################################
## Test the extraction of normal sample VAFs from a VCF file
###############################################################################

# Extract relevant columns from the Annovar output with consistent structure
echo $(date +"%F") $(date +"%T") "Processing Annovar output with normal sample analysis..."

# multi_annno=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T/DFSP-001-T.hg38_multianno.txt
cur_dir=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T
tumour_id=DFSP-001-T


# Extract relevant columns from the Annovar output with consistent structure
echo $(date +"%F") $(date +"%T") "Processing Annovar output with normal sample analysis..."

cur_dir=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T
tumour_id=DFSP-001-T
bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
normal_id="DFSP-001-N"

# Debug first to see what's happening
echo "=== DEBUGGING VARIANT KEYS ==="
echo "First VCF variant:"
bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${cur_dir}/${tumour_id}.final.vcf.gz | head -1
echo "First Annovar variant:"
awk 'NR==2 {print $1":"$2":"$4":"$5}' ${cur_dir}/${tumour_id}.hg38_multianno.txt

# Extract relevant columns from the Annovar output with consistent structure
echo $(date +"%F") $(date +"%T") "Processing Annovar output with normal sample analysis..."

cur_dir=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T
tumour_id=DFSP-001-T
bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
normal_id="DFSP-001-N"

# Debug first to see what's happening
echo "=== DEBUGGING VARIANT KEYS ==="
echo "First VCF variant:"
bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${cur_dir}/${tumour_id}.final.vcf.gz | head -1
echo "First Annovar variant:"
awk 'NR==2 {print $1":"$2":"$4":"$5}' ${cur_dir}/${tumour_id}.hg38_multianno.txt

# Extract relevant columns from the Annovar output with consistent structure
echo $(date +"%F") $(date +"%T") "Processing Annovar output with normal sample analysis..."

cur_dir=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T
tumour_id=DFSP-001-T
bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
normal_id="DFSP-001-N"

# Debug first to see what's happening
echo "=== DEBUGGING VARIANT KEYS ==="
echo "First VCF variant:"
bcftools query -f '%CHROM:%POS:%REF:%ALT\n' ${cur_dir}/${tumour_id}.final.vcf.gz | head -1
echo "First Annovar variant:"
awk 'NR==2 {print $1":"$2":"$4":"$5}' ${cur_dir}/${tumour_id}.hg38_multianno.txt

echo "First variant FORMAT field:"
awk 'NR==2 {print "Full last field: " $NF; split($NF, a, ":"); print "a[1]=" a[1] " a[2]=" a[2] " a[3]=" a[3] " a[4]=" a[4]}' ${cur_dir}/${tumour_id}.hg38_multianno.txt

# Extract normal sample data from Annovar output
echo $(date +"%F") $(date +"%T") "Processing Annovar output to extract normal sample data..."

cur_dir=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T
tumour_id=DFSP-001-T
normal_id="DFSP-001-N"

# Extract normal sample data from Annovar output
echo $(date +"%F") $(date +"%T") "Processing Annovar output to extract normal sample data..."

cur_dir=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/mutect2_with_black_repeat_filter_new/DFSP-001-T
tumour_id=DFSP-001-T
bam_dir="/home/zhonggr/projects/250224_DFSP_WES/data/wes/preprocessing/recalibrated"
normal_id="DFSP-001-N"

# Check if this is a paired sample
if [ -d ${bam_dir}/${normal_id} ]; then
    echo "Processing paired sample: ${tumour_id} with normal ${normal_id}"
    
    # For paired samples: extract normal sample data from the Annovar output
    # The last field contains format like: GT:AD:AF:DP for both normal and tumor
    less -S ${cur_dir}/${tumour_id}.hg38_multianno.txt | \
    awk 'BEGIN {FS=OFS="\t"}
    NR==1 {
        # Add original AD, AF, DP columns plus tumor and normal-specific columns
        print $0, "AD", "AF", "DP", "TUMOR_AD", "TUMOR_AF", "TUMOR_DP", "NORMAL_AD", "NORMAL_AF", "NORMAL_DP", "SAMPLE_TYPE"
    }
    NR > 1 {
        # The last field contains sample data for tumor (and normal if paired)
        # Format: tumor_GT:tumor_AD:tumor_AF:tumor_DP:normal_GT:normal_AD:normal_AF:normal_DP:...
        
        # Split the last field by colons
        split($NF, format_fields, ":")
        
        # Extract tumor data (same as original script)
        tumor_ad = format_fields[2]   # AD is position 2
        tumor_af = format_fields[3]   # AF is position 3  
        tumor_dp = format_fields[4]   # DP is position 4
        
        # Extract normal data (assuming normal data follows tumor data)
        # This depends on how mutect2 outputs the format field for paired samples
        # Common patterns: GT:AD:AF:DP (tumor) then GT:AD:AF:DP (normal)
        # Or: GT:AD:AF:DP:F1R2:F2R1:... (tumor) then similar for normal
        
        # Try to find normal data - usually starts after position 4
        normal_ad = "NA"
        normal_af = "NA" 
        normal_dp = "NA"
        
        # Look for normal sample data pattern
        # If there are more than 8 fields, normal data likely starts around position 6-8
        if (length(format_fields) > 8) {
            # Common pattern: positions 6, 7, 8 for normal AD, AF, DP
            normal_ad = format_fields[6]
            normal_af = format_fields[7]
            normal_dp = format_fields[8]
        }
        
        # Print original data + original AD,AF,DP + tumor columns + normal columns
        print $0, tumor_ad, tumor_af, tumor_dp, tumor_ad, tumor_af, tumor_dp, normal_ad, normal_af, normal_dp, "PAIRED"
    }' > ${cur_dir}/${tumour_id}.annovar_test.txt

else
    echo "Processing tumor-only sample: ${tumour_id}"
    
    # For tumor-only samples: no normal data available
    less -S ${cur_dir}/${tumour_id}.hg38_multianno.txt | \
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
        
        # Print original data + original AD,AF,DP + tumor columns + normal columns (NA)
        print $0, tumor_ad, tumor_af, tumor_dp, tumor_ad, tumor_af, tumor_dp, "NA", "NA", "NA", "TUMOR_ONLY"
    }' > ${cur_dir}/${tumour_id}.annovar_test.txt

fi

# Debug: Check the format field structure
echo ""
echo "=== DEBUGGING FORMAT FIELD ==="
echo "First variant format field:"
awk 'NR==2 {
    split($NF, a, ":")
    print "Full format field: " $NF
    print "Number of colon-separated fields: " length(a)
    for(i=1; i<=length(a); i++) {
        print "Field " i ": " a[i]
    }
}' ${cur_dir}/${tumour_id}.hg38_multianno.txt

# Verification
echo ""
echo "=== VERIFICATION ==="
echo "Check data extraction (first 3 variants):"
awk 'NR>1 && NR<=4 {
    n = NF
    orig_ad = $(n-10)
    orig_af = $(n-9) 
    orig_dp = $(n-8)
    tumor_ad = $(n-7)
    tumor_af = $(n-6)
    tumor_dp = $(n-5)
    normal_ad = $(n-4)
    normal_af = $(n-3)
    normal_dp = $(n-2)
    
    print "Variant " NR-1 ":"
    print "  Original: AD=" orig_ad " AF=" orig_af " DP=" orig_dp
    print "  Tumor:    AD=" tumor_ad " AF=" tumor_af " DP=" tumor_dp  
    print "  Normal:   AD=" normal_ad " AF=" normal_af " DP=" normal_dp
    print "  Match:    " (orig_ad == tumor_ad && orig_af == tumor_af && orig_dp == tumor_dp ? "YES" : "NO")
    print ""
}' ${cur_dir}/${tumour_id}.annovar_test.txt
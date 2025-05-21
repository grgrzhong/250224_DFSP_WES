#!/bin/bash

# Activate conda environment
conda activate renv

# Export facets segments
input="/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/cnv/facets"
Rscript ./modules/R/export_facets.R --input ${input}

# Annotate segments
conda activate varcall

seg_file=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/cnv/facets/DFSP-362-T-P2/DFSP-362-T-P2.filtered.seg
annotation=/home/zhonggr/projects/250224_DFSP_WES/data/reference/Gencode/annotation_protein_coding.bed
output_file=/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_calling/cnv/facets/DFSP-362-T-P2/DFSP-362-T-P2.filtered.annotated.tsv

cat -t ${seg_file}

bedtools intersect -wa -wb \
    -b ${annotation} \
    -a ${seg_file} \
    > ${output_file}


# Create a temporary BED file from SEG file
temp_bed=$(mktemp)

# Convert SEG to BED format with special handling for chromosome names
awk -F'\t' '
  # Skip header
  NR==1 {next}
  
  # Process each data row
  {
    # Extract fields
    sample=$1
    chrom=$2
    start=$3-1  # Convert to 0-based
    end=$4
    score=$5
    segcn=$6
    svtype=$7
    
    # Handle special chromosome cases
    if (chrom == "NA") {
      # Skip rows with NA chromosomes - they cause problems with bedtools
      next
    }
    else if (chrom == "X") {
      # Keep X as-is
      chrom="X"
    }
    else if (chrom == "Y") {
      # Keep Y as-is
      chrom="Y"
    }
    
    # Print in BED format
    print chrom"\t"start"\t"end"\t"sample"\t"score"\t"segcn"\t"svtype
  }
' "${seg_file}" > "${temp_bed}"

# Check if conversion was successful
if [ ! -s "${temp_bed}" ]; then
    echo "Error: Failed to convert SEG file to BED format"
    exit 1
fi

echo "Successfully created temporary BED file for bedtools"
echo "First few lines of the BED file:"
head -n 3 "${temp_bed}"
echo "Total lines in BED file: $(wc -l < ${temp_bed})"

# Check if annotation file exists
if [ ! -f "$annotation" ]; then
    echo "ERROR: Annotation file not found: $annotation"
    exit 1
fi

# Show first few lines of annotation file
echo "First few lines of annotation file:"
head -n 3 "$annotation"

# Run bedtools intersect with proper format
bedtools intersect -wa -wb \
    -a "${temp_bed}" \
    -b "${annotation}" \
    > "${output_file}"

# Check if the intersect was successful
if [ $? -eq 0 ] && [ -s "${output_file}" ]; then
    echo "Successfully created annotated file: ${output_file}"
    
    # Count lines in the output file
    count=$(wc -l < "${output_file}")
    echo "Found ${count} gene-CNV intersections"
else
    echo "Error: bedtools intersect failed or produced empty output"
    # Show any error output
    bedtools intersect -wa -wb -a "${temp_bed}" -b "${annotation}" 2>&1 | head -n 5
    exit 1
fi

# Clean up temporary file
rm "${temp_bed}"

echo "Annotation complete!"
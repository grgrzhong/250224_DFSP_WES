#!/bin/bash

# Set directories
PROJECT_DIR="/home/zhonggr/projects/250224_DFSP_WES"
# WORK_DIR="data/WES/DFSP"
WORK_DIR="data/WES/SARC"
BAM_DIR=$PROJECT_DIR/$WORK_DIR/BAM

SAMPLE_SHEET="${PROJECT_DIR}/${WORK_DIR}/samplesheet_sarc.csv"

# Ensure output directory exists
mkdir -p "$(dirname "$SAMPLE_SHEET")"

# Create temporary file
temp_file="${SAMPLE_SHEET}.tmp"

# Create header
echo "patient,sample,status,bam,bai" > "$temp_file"

# Process each sample directory
for sample_dir in $(ls "$BAM_DIR"); do
    # Extract patient ID by removing everything after and including the last dash
    patient=$(echo "$sample_dir" | sed -E 's/^(.*?)-(T|N).*$/\1/')
    
    # Determine status (0 for normal, 1 for tumor) based on sample name
    if [[ "$sample_dir" == *"-N"* ]]; then
        status="0"
    elif [[ "$sample_dir" == *"-T"* ]]; then
        status="1"
    else
        echo "Warning: Cannot determine status for $sample_dir, skipping"
        continue
    fi
    
    # Find the BAM file
    bam_file="${BAM_DIR}/${sample_dir}/${sample_dir}_recalibrated.bam"
    bai_file="${bam_file}_recalibrated.bai"
    
    # Check if files exist
    if [[ ! -f "$bam_file" ]]; then
        echo "Warning: BAM file not found: $bam_file"
        # Try to find BAM file in the sample directory
        bam_file=$(find "${BAM_DIR}/${sample_dir}" -name "*.bam" | head -n 1)
        if [[ -z "$bam_file" ]]; then
            echo "Error: No BAM file found for $sample_dir, skipping"
            continue
        fi
    fi
    
    # Check for BAI file
    if [[ ! -f "$bai_file" ]]; then
        # Try with .bai extension
        bai_file="${bam_file%.bam}.bai"
        if [[ ! -f "$bai_file" ]]; then
            echo "Warning: BAI file not found for $bam_file"
            # Try to find BAI file in the sample directory
            bai_file=$(find "${BAM_DIR}/${sample_dir}" -name "*.bai" | head -n 1)
            if [[ -z "$bai_file" ]]; then
                echo "Error: No BAI file found for $sample_dir, skipping"
                continue
            fi
        fi
    fi
    
    # Add to sample sheet
    echo "${patient},${sample_dir},${status},${bam_file},${bai_file}" >> "$temp_file"
done

# Sort file (keeping header)
(head -n 1 "$temp_file" && tail -n +2 "$temp_file" | sort -t',' -k1,1 -k2,2) > "$SAMPLE_SHEET"
rm -f "$temp_file"

echo "Created samplesheet at: $SAMPLE_SHEET"

# Print statistics
echo -e "\nSamplesheet statistics:"
echo "Total patients: $(cut -d',' -f1 "$SAMPLE_SHEET" | tail -n +2 | sort -u | wc -l)"
echo "Total samples: $(tail -n +2 "$SAMPLE_SHEET" | wc -l)"
echo "Total tumor samples: $(grep -c ",1," "$SAMPLE_SHEET")"
echo "Total normal samples: $(grep -c ",0," "$SAMPLE_SHEET")"

# Print patients without normal samples
echo -e "\nPatients without normal samples:"
cut -d',' -f1,3 "$SAMPLE_SHEET" | tail -n +2 | awk -F',' '
    {
        patient[$1]++
        if ($2 == "0") normal[$1]++
    }
    END {
        for (p in patient) {
            if (!(p in normal)) print p
        }
    }'
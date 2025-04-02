#!/bin/bash

# filepath: /home/zhonggr/projects/250224_DFSP_WES/scripts/create_samplesheet.sh

# Set directories
data_dir="${1:-/home/zhonggr/projects/250224_DFSP_WES/data/WES/DFSP}"
bam_dir="${data_dir}/bam"
fastq_dir="${data_dir}/fastq"  # This contains all FASTQ files, not organized by sample

SAMPLE_SHEET="${data_dir}/samplesheet.csv"

# Ensure output directory exists
mkdir -p "$(dirname "$SAMPLE_SHEET")"

# Create temporary file
temp_file="${SAMPLE_SHEET}.tmp"

# Create header
echo "patient,sample,status,fastq_1,fastq_2,bam,bai" > "$temp_file"

# Process each sample directory in BAM_DIR
for sample_dir in $(ls "$bam_dir"); do
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
    
    # Initialize variables
    fastq_1=""
    fastq_2=""
    bam_file=""
    bai_file=""
    
    # Find the BAM file
    potential_bam="${bam_dir}/${sample_dir}/${sample_dir}_recalibrated.bam"
    if [[ -f "$potential_bam" ]]; then
        bam_file="$potential_bam"
    else
        # Try to find any BAM file in the sample directory
        found_bam=$(find "${bam_dir}/${sample_dir}" -name "*.bam" | head -n 1)
        if [[ -n "$found_bam" ]]; then
            bam_file="$found_bam"
        else
            echo "Warning: No BAM file found for $sample_dir"
        fi
    fi
    
    # Find the BAI file
    if [[ -n "$bam_file" ]]; then
        # Check for BAM.BAI file
        potential_bai="${bam_file}.bai"
        if [[ -f "$potential_bai" ]]; then
            bai_file="$potential_bai"
        else
            # Check for BAM file with .bai extension
            potential_bai="${bam_file%.bam}.bai"
            if [[ -f "$potential_bai" ]]; then
                bai_file="$potential_bai"
            else
                # Try to find any BAI file in the sample directory
                found_bai=$(find "${bam_dir}/${sample_dir}" -name "*.bai" | head -n 1)
                if [[ -n "$found_bai" ]]; then
                    bai_file="$found_bai"
                else
                    echo "Warning: No BAI file found for $sample_dir"
                fi
            fi
        fi
    fi
    
    # Find FASTQ files - now looking in the main fastq directory
    # We'll search for files that contain the sample name in their filename
    echo "Looking for FASTQ files for sample: $sample_dir"
    
    # Search patterns to try (from most specific to least specific)
    patterns=(
        "${sample_dir}_R1"       # Sample_R1 pattern
        "${sample_dir}_1"        # Sample_1 pattern
        "${sample_dir}.*_R1"     # Sample followed by anything then _R1
        "${sample_dir}.*_1"      # Sample followed by anything then _1
        "${sample_dir}"          # Just the sample name
    )
    
    # Find R1 (read 1) FASTQ file
    for pattern in "${patterns[@]}"; do
        if [[ -z "$fastq_1" ]]; then
            # Try to find files with this pattern
            potential_files=$(find "$fastq_dir" -name "*${pattern}*.fastq.gz" 2>/dev/null)
            
            # If no .fastq.gz files, try .fq.gz
            if [[ -z "$potential_files" ]]; then
                potential_files=$(find "$fastq_dir" -name "*${pattern}*.fq.gz" 2>/dev/null)
            fi
            
            # Choose the first matching file
            for file in $potential_files; do
                # Check if this file looks like an R1/read 1 file
                if [[ "$file" =~ _R1_ || "$file" =~ _1\. ]]; then
                    fastq_1="$file"
                    echo "  Found R1: $fastq_1"
                    break
                fi
            done
            
            # If we found a match, break the outer loop too
            if [[ -n "$fastq_1" ]]; then
                break
            fi
        fi
    done
    
    # Find R2 (read 2) FASTQ file - similar approach to R1
    # If we found an R1 file, try to find the matching R2 by replacing R1/1 with R2/2
    if [[ -n "$fastq_1" ]]; then
        # First try a direct replacement
        potential_r2="${fastq_1/_R1_/_R2_}"
        if [[ -f "$potential_r2" ]]; then
            fastq_2="$potential_r2"
            echo "  Found R2 (direct replacement): $fastq_2"
        else
            potential_r2="${fastq_1/_1\./_2\.}"
            if [[ -f "$potential_r2" ]]; then
                fastq_2="$potential_r2"
                echo "  Found R2 (direct replacement): $fastq_2"
            fi
        fi
    fi
    
    # If we still don't have an R2 file, search for it as we did for R1
    if [[ -z "$fastq_2" ]]; then
        patterns=(
            "${sample_dir}_R2"
            "${sample_dir}_2"
            "${sample_dir}.*_R2"
            "${sample_dir}.*_2"
            "${sample_dir}"
        )
        
        for pattern in "${patterns[@]}"; do
            if [[ -z "$fastq_2" ]]; then
                potential_files=$(find "$fastq_dir" -name "*${pattern}*.fastq.gz" 2>/dev/null)
                
                # If no .fastq.gz files, try .fq.gz
                if [[ -z "$potential_files" ]]; then
                    potential_files=$(find "$fastq_dir" -name "*${pattern}*.fq.gz" 2>/dev/null)
                fi
                
                for file in $potential_files; do
                    if [[ "$file" =~ _R2_ || "$file" =~ _2\. ]]; then
                        fastq_2="$file"
                        echo "  Found R2: $fastq_2"
                        break
                    fi
                done
                
                if [[ -n "$fastq_2" ]]; then
                    break
                fi
            fi
        done
    fi
    
    # Last resort: try to find any remaining FASTQ files that might match by listing all related files
    if [[ -z "$fastq_1" || -z "$fastq_2" ]]; then
        echo "  Still searching for FASTQ files..."
        # Get all files that might relate to this sample
        related_files=$(find "$fastq_dir" -name "*${sample_dir}*" | grep -E '\.fastq\.gz$|\.fq\.gz$' | sort)
        
        # If we have exactly two files, and we need both R1 and R2, just use them
        if [[ -z "$fastq_1" && -z "$fastq_2" ]]; then
            file_count=$(echo "$related_files" | wc -l)
            if [[ $file_count -eq 2 ]]; then
                fastq_1=$(echo "$related_files" | head -n 1)
                fastq_2=$(echo "$related_files" | tail -n 1)
                echo "  Using the only two related files found:"
                echo "    R1: $fastq_1"
                echo "    R2: $fastq_2" 
            else
                echo "  Found $file_count related files, cannot determine which to use."
            fi
        # If we just need R1, use the first file
        elif [[ -z "$fastq_1" ]]; then
            fastq_1=$(echo "$related_files" | grep -E '_R1_|_1\.' | head -n 1)
            echo "  Found R1 from related files: $fastq_1"
        # If we just need R2, use the first file that looks like R2
        elif [[ -z "$fastq_2" ]]; then
            fastq_2=$(echo "$related_files" | grep -E '_R2_|_2\.' | head -n 1)
            echo "  Found R2 from related files: $fastq_2"
        fi
    fi
    
    # Add to sample sheet (use empty strings if files aren't found)
    echo "${patient},${sample_dir},${status},${fastq_1},${fastq_2},${bam_file},${bai_file}" >> "$temp_file"
    
    # Report what we found for this sample
    echo "Sample ${sample_dir} added with:"
    if [[ -n "$fastq_1" ]]; then echo "  FASTQ_1: $fastq_1"; else echo "  FASTQ_1: Missing"; fi
    if [[ -n "$fastq_2" ]]; then echo "  FASTQ_2: $fastq_2"; else echo "  FASTQ_2: Missing"; fi
    if [[ -n "$bam_file" ]]; then echo "  BAM: $bam_file"; else echo "  BAM: Missing"; fi
    if [[ -n "$bai_file" ]]; then echo "  BAI: $bai_file"; else echo "  BAI: Missing"; fi
    echo "-----------------------------------------"
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

# Count samples with missing files
echo "Samples missing FASTQ_1: $(awk -F, '{if(NR>1 && $4=="") print}' "$SAMPLE_SHEET" | wc -l)"
echo "Samples missing FASTQ_2: $(awk -F, '{if(NR>1 && $5=="") print}' "$SAMPLE_SHEET" | wc -l)"
echo "Samples missing BAM: $(awk -F, '{if(NR>1 && $6=="") print}' "$SAMPLE_SHEET" | wc -l)"
echo "Samples missing BAI: $(awk -F, '{if(NR>1 && $7=="") print}' "$SAMPLE_SHEET" | wc -l)"

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

# Print sample counts by patient
echo -e "\nSample counts by patient:"
cut -d',' -f1,3 "$SAMPLE_SHEET" | tail -n +2 | awk -F',' '
    {
        patient[$1]++
        if ($2 == "0") normal[$1]++
        else if ($2 == "1") tumor[$1]++
    }
    END {
        printf "%-20s %-15s %-15s %-15s\n", "Patient", "Total", "Normal", "Tumor"
        for (p in patient) {
            printf "%-20s %-15s %-15s %-15s\n", p, patient[p], normal[p] ? normal[p] : 0, tumor[p] ? tumor[p] : 0
        }
    }'
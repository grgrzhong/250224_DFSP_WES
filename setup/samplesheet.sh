#!/bin/bash

# filepath: /home/zhonggr/projects/250224_DFSP_WES/samplesheet.sh

# Set directories
data_dir="${1:-/home/zhonggr/projects/250224_DFSP_WES/data/sarc}"
# data_dir="${1:-/home/zhonggr/projects/250224_DFSP_WES/data/wes}"
fastq_dir="${data_dir}/preprocessing/fastq"  # This contains all FASTQ files, not organized by sample
bam_dir="${data_dir}/preprocessing/recalibrated"

SAMPLE_SHEET="${data_dir}/csv/samplesheet.csv"

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
    echo "Looking for FASTQ files for sample: $sample_dir"
    
    # Define common patterns for R1/R2 files
    R1_patterns=("_R1_" "_R1." "_1.")
    R2_patterns=("_R2_" "_R2." "_2.")
    
    # Find R1 FASTQ file
    for r1_pattern in "${R1_patterns[@]}"; do
        if [[ -z "$fastq_1" ]]; then
            # Look for files matching both the sample name and the R1 pattern
            found_files=$(find "$fastq_dir" -name "*${sample_dir}*${r1_pattern}*" \( -name "*.fastq.gz" -o -name "*.fq.gz" \) 2>/dev/null | sort)
            
            if [[ -n "$found_files" ]]; then
                # Take the first matching file
                fastq_1=$(echo "$found_files" | head -n 1)
                echo "  Found R1: $fastq_1 (pattern: $r1_pattern)"
                break
            fi
        fi
    done
    
    # Find R2 FASTQ file - if we found R1, try to find its pair
    if [[ -n "$fastq_1" ]]; then
        # Try to find the matching R2 based on the R1 filename
        for i in "${!R1_patterns[@]}"; do
            r1_pattern="${R1_patterns[$i]}"
            r2_pattern="${R2_patterns[$i]}"
            
            if [[ "$fastq_1" == *"$r1_pattern"* ]]; then
                potential_r2="${fastq_1/$r1_pattern/$r2_pattern}"
                if [[ -f "$potential_r2" ]]; then
                    fastq_2="$potential_r2"
                    echo "  Found R2 (paired with R1): $fastq_2"
                    break
                fi
            fi
        done
    fi
    
    # If we still don't have an R2, search independently
    if [[ -z "$fastq_2" ]]; then
        for r2_pattern in "${R2_patterns[@]}"; do
            if [[ -z "$fastq_2" ]]; then
                found_files=$(find "$fastq_dir" -name "*${sample_dir}*${r2_pattern}*" \( -name "*.fastq.gz" -o -name "*.fq.gz" \) 2>/dev/null | sort)
                
                if [[ -n "$found_files" ]]; then
                    fastq_2=$(echo "$found_files" | head -n 1)
                    echo "  Found R2: $fastq_2 (pattern: $r2_pattern)"
                    break
                fi
            fi
        done
    fi
    
    # Additional validation - check if R1 and R2 appear to be properly paired
    if [[ -n "$fastq_1" && -n "$fastq_2" ]]; then
        # Extract base names for comparison (removing R1/R2 and extension parts)
        f1_base=$(basename "$fastq_1" | sed -E 's/_R[12]_.*|_[12]\..*//')
        f2_base=$(basename "$fastq_2" | sed -E 's/_R[12]_.*|_[12]\..*//')
        
        # Ensure the base parts match
        if [[ "$f1_base" != "$f2_base" ]]; then
            echo "  WARNING: R1 and R2 files appear to be from different samples!"
            echo "    R1 base: $f1_base"
            echo "    R2 base: $f2_base"
        fi
    fi
    
    # Last resort: if we still don't have both FASTQ files, try to find any files containing the sample name
    if [[ -z "$fastq_1" || -z "$fastq_2" ]]; then
        echo "  Still searching for FASTQ files..."
        # Get all files that might relate to this sample
        related_files=$(find "$fastq_dir" -name "*${sample_dir}*" | grep -E '\.fastq\.gz$|\.fq\.gz$' | sort)
        
        # If we have exactly two files, and we need both R1 and R2, just use them
        if [[ -z "$fastq_1" && -z "$fastq_2" ]]; then
            file_count=$(echo "$related_files" | wc -l)
            if [[ $file_count -eq 2 ]]; then
                # Try to determine which is R1 and which is R2
                r1_file=""
                r2_file=""
                while read -r file; do
                    for r1p in "${R1_patterns[@]}"; do
                        if [[ "$file" == *"$r1p"* ]]; then
                            r1_file="$file"
                            break
                        fi
                    done
                    for r2p in "${R2_patterns[@]}"; do
                        if [[ "$file" == *"$r2p"* ]]; then
                            r2_file="$file"
                            break
                        fi
                    done
                done <<< "$related_files"
                
                # If we identified both, use them
                if [[ -n "$r1_file" && -n "$r2_file" ]]; then
                    fastq_1="$r1_file"
                    fastq_2="$r2_file"
                # Otherwise just use first and second
                else
                    fastq_1=$(echo "$related_files" | head -n 1)
                    fastq_2=$(echo "$related_files" | tail -n 1)
                fi
                echo "  Using two related files found:"
                echo "    R1: $fastq_1"
                echo "    R2: $fastq_2" 
            else
                echo "  Found $file_count related files, cannot determine which to use."
            fi
        # If we just need R1, use the first file that matches R1 patterns
        elif [[ -z "$fastq_1" ]]; then
            for r1p in "${R1_patterns[@]}"; do
                r1_file=$(echo "$related_files" | grep "$r1p" | head -n 1)
                if [[ -n "$r1_file" ]]; then
                    fastq_1="$r1_file"
                    echo "  Found R1 from related files: $fastq_1"
                    break
                fi
            done
            # If still not found, use first file
            if [[ -z "$fastq_1" && -n "$related_files" ]]; then
                fastq_1=$(echo "$related_files" | head -n 1)
                echo "  Using first related file as R1: $fastq_1"
            fi
        # If we just need R2, use the first file that matches R2 patterns
        elif [[ -z "$fastq_2" ]]; then
            for r2p in "${R2_patterns[@]}"; do
                r2_file=$(echo "$related_files" | grep "$r2p" | head -n 1)
                if [[ -n "$r2_file" ]]; then
                    fastq_2="$r2_file"
                    echo "  Found R2 from related files: $fastq_2"
                    break
                fi
            done
            # If still not found, use first file not matching fastq_1
            if [[ -z "$fastq_2" && -n "$related_files" ]]; then
                fastq_2=$(echo "$related_files" | grep -v "$fastq_1" | head -n 1)
                echo "  Using non-R1 related file as R2: $fastq_2"
            fi
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
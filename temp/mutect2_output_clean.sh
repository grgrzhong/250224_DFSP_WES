#!/bin/bash

#############################################################################
# Clean the funcotator output files
#############################################################################
# Create destination directory if it doesn't exist
mkdir -p ./data/wes/annotation/funcotator/

# Loop through each sample directory
for sample_dir in ./data/wes/variant_calling/mutect2/DFSP-*; do
    sample=$(basename "$sample_dir")
    echo "Processing sample: $sample"
    
    # Copy annotated.maf.gz files
    find "$sample_dir" -name "*annotated.maf.gz" -exec cp {} ./data/wes/annotation/funcotator/${sample} \;
    
    # Copy annotated.tsv files
    find "$sample_dir" -name "*annotated.tsv" -exec cp {} ./data/wes/annotation/funcotator/${sample} \;
done

echo "Files copied successfully."

#############################################################################
# Clean the annovar output files
#############################################################################
# Create destination directory if it doesn't exist
mkdir -p ./data/wes/annotation/annovar/

# Loop through each sample directory
for sample_dir in ./data/wes/variant_calling/mutect2/DFSP-*; do
    sample=$(basename "$sample_dir")
    echo "Processing sample: $sample"
    
    # Copy annotated.maf.gz files
    find "$sample_dir" -name "*Annovar.log" -exec cp {} ./data/wes/annotation/annovar/${sample} \;
    
    find "$sample_dir" -name "*.avinput" -exec cp {} ./data/wes/annotation/annovar/${sample} \;
    
    find "$sample_dir" -name "*hg38_multianno.txt" -exec cp {} ./data/wes/annotation/annovar/${sample} \;
    
    find "$sample_dir" -name "*hg38_multianno.vcf" -exec cp {} ./data/wes/annotation/annovar/${sample} \;
    
    find "$sample_dir" -name "*annovar.txt" -exec cp {} ./data/wes/annotation/annovar/${sample} \;

done

echo "Files copied successfully."


#############################################################################
# Clean the funcotator output files
#############################################################################
BASE_DIR="/lustre1/g/path_my/250224_DFSP_WES/data/wes/annotation/annovar"

# Find all sample directories
for sample_dir in "$BASE_DIR"/DFSP-*; do
    if [ -d "$sample_dir" ]; then
        echo "Processing directory: $(basename "$sample_dir")"
        
        # Find files with _annovar.txt pattern and rename them
        find "$sample_dir" -name "*_annovar.txt" | while read file; do
            dir=$(dirname "$file")
            filename=$(basename "$file")
            new_filename="${filename/_annovar/}"
            new_path="$dir/$new_filename"
            
            echo "Renaming: $filename → $new_filename"
            mv "$file" "$new_path"
        done
    fi
done

echo "Renaming complete."

#############################################################################
# Remove specific .txt files
#############################################################################
BASE_DIR="/lustre1/g/path_my/250224_DFSP_WES/data/wes/annotation/annovar"

# Find all sample directories
for sample_dir in "$BASE_DIR"/DFSP-*; do
    if [ -d "$sample_dir" ]; then
        sample=$(basename "$sample_dir")
        echo "Processing directory: $sample"
        
        # Find and remove only the plain .txt files (not those with _annovar or multianno)
        target_file="$sample_dir/$sample.txt"
        if [ -f "$target_file" ]; then
            echo "Removing: $target_file"
            rm -f "$target_file"
        else
            echo "File not found: $target_file"
        fi
    fi
done

echo "Removal complete."
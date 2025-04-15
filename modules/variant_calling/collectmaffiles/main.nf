#!/usr/bin/env nextflow
// filepath: modules/maf_merger.nf

process MERGE_MAF_FILES {
    tag "merge_maf"
    publishDir "${params.outdir}/merged_mafs", mode: 'copy'
    
    container params.containers.python // If using containers
    
    input:
    path(maf_files)
    val(output_name)
    
    output:
    path "${output_name}", emit: merged_maf
    path "merge_stats.txt", emit: stats
    
    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import glob
    import os
    
    maf_files = ["${maf_files.join('", "')}"]
    output_name = "${output_name}"
    
    print(f"Processing {len(maf_files)} MAF files")
    
    # Store header lines from the first file
    header_lines = []
    with open(maf_files[0], 'r') as first_file:
        for line in first_file:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                break
    
    # Initialize merged dataframe
    merged_df = None
    file_stats = {}
    
    # Process each MAF file
    for maf_file in maf_files:
        print(f"Processing file: {maf_file}")
        try:
            # Count header lines to skip
            lines_to_skip = 0
            with open(maf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        lines_to_skip += 1
                    else:
                        break
            
            # Read the MAF file
            df = pd.read_csv(maf_file, sep='\t', comment='#', skiprows=lines_to_skip)
            file_stats[os.path.basename(maf_file)] = len(df)
            
            # Add to merged dataframe
            if merged_df is None:
                merged_df = df.copy()
            else:
                merged_df = pd.concat([merged_df, df], ignore_index=True)
                
        except Exception as e:
            print(f"Error processing {maf_file}: {e}")
            file_stats[os.path.basename(maf_file)] = "ERROR: " + str(e)
    
    # Write the merged MAF file with headers
    with open(output_name, 'w') as f:
        # Write header comments
        for header in header_lines:
            f.write(header)
    
    # Append the data
    merged_df.to_csv(output_name, sep='\t', index=False, mode='a')
    
    # Write statistics
    with open("merge_stats.txt", 'w') as stats_file:
        stats_file.write(f"Total MAF files processed: {len(maf_files)}\\n")
        stats_file.write(f"Total variants in merged file: {len(merged_df)}\\n\\n")
        stats_file.write("Individual file statistics:\\n")
        for file, count in file_stats.items():
            stats_file.write(f"{file}: {count}\\n")
    
    print(f"Merged MAF file created with {len(merged_df)} variants")
    """
}

#!/usr/bin/env nextflow
// filepath: main.nf

nextflow.enable.dsl=2

// Import modules
include { MERGE_MAF_FILES } from './modules/maf_merger'

// Default parameters
params.input_dir = null
params.outdir = 'results'
params.containers = [python: 'python:3.9']

// Validate inputs
if (!params.input_dir) {
    error "Input directory not specified. Please provide --input_dir parameter."
}

workflow {
    // Find all MAF files according to expected pattern
    maf_files_ch = Channel.fromPath("${params.input_dir}/**/variants/*.maf")
                         .collect()
                         .ifEmpty { error "No MAF files found in ${params.input_dir}" }
    
    // Run the merge process
    MERGE_MAF_FILES(maf_files_ch, "merged_cohort.maf")
}
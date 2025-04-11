process TAG_UMI {
    tag "${meta.id}"
    label 'process_medium'
    
    // container "${params.singularity_container_dir}/pysam-0.23.0.sif"
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("*.umi.bam"), emit: bam
    path "versions.yml", emit: versions
    
    script:

    """
    #!/usr/bin/env python3

    import pysam
    import argparse

    input_bam = "${bam}"
    output_bam = "${meta.id}.umi.bam"

    with pysam.AlignmentFile(input_bam, "rb") as infile, \\
        pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
        for read in infile:
            read_name = read.query_name
            
            # Extract UMI from read name
            umi = read_name.split(":")[-1]
            
            # Replace the UMI delimiter
            umi = umi.replace("_", "-")
            
            # Remove UMI from read names
            read.query_name = ":".join(read.query_name.split(":")[:-1])
            
            # Set UMI as a custom tag "RX"
            read.set_tag("RX", umi)
            outfile.write(read)
            
    # Create versions.yml file
    with open("versions.yml", "w") as f:
        f.write("TAG_UMI:\\n" + 
                "  python: \$(python --version | sed 's/Python //g')\\n" + 
                "  pysam: \$(python -c 'import pysam; print(pysam.__version__)')\\n")
    """
}
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Define required parameters
params.singularity_container_dir = "/home/zhonggr/projects/250224_DFSP_WES/containers/singularity"
params.genome = "GRCh38"
params.genomes = [
    "GRCh38": [
        "fasta": "/home/zhonggr/projects/250224_DFSP_WES/data/Reference/Gencode/gencode.hg38.v36.primary_assembly.fa"
    ]
]
params.outdir = "results"

process BWA_MEM {

    // Use Docker container as fallback
    container "/home/zhonggr/projects/250224_DFSP_WES/containers/singularity/bwa-0.7.19.sif"
    
    publishDir([
        [
            path: { "${params.outdir}/preprocessing/bam/test}" },
            mode: "copy",
            pattern: "*.bam"
        ]
    ])
    
    input:
    tuple val(meta), path(reads1), path(reads2)
    path fasta
    path index
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions
    path "debug.log", emit: log
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    echo "Starting BWA for ${prefix}" > debug.log
    echo "Reads: ${reads1} ${reads2}" >> debug.log
    echo "Reference: ${fasta}" >> debug.log
    echo "Index files: \$(ls -la)" >> debug.log
    
    bwa mem \\
        -M \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${prefix}\\tLB:XGenV2\\tPL:ILLUMINA\\tPM:NOVASEQ\\tSM:${prefix}\\tPU:NA" \\
        ${fasta} \\
        ${reads1} ${reads2} | \\
        samtools view -Sb - > ${prefix}.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -e '^Version' | sed 's/Version: //g')
        samtools: \$(samtools --version | grep -e '^samtools' | sed 's/samtools //g')
    END_VERSIONS
    """
}

workflow {
    // Create test sample with absolute paths
    ch_reads = Channel
    .of([
        [id: 'test_sample'], 
        file('/lustre1/g/path_my/250224_DFSP_WES/work/55/41ec08660aea6d828cabbaa794db5b/DFSP-261-N_trimmed_1.fastq.gz'),
        file('/lustre1/g/path_my/250224_DFSP_WES/work/55/41ec08660aea6d828cabbaa794db5b/DFSP-261-N_trimmed_2.fastq.gz')
    ])
    
    ch_reads.view { println "Input reads: ${it}" }

    // Reference - properly check if file exists
    ref_file = file("/home/zhonggr/projects/250224_DFSP_WES/data/Reference/Gencode/gencode.hg38.v36.primary_assembly.fa")
    
    // Check file existence correctly - FIX HERE
    if(!ref_file.exists()) {  // Changed from file(ref_path).exists() to ref_file.exists()
        println "ERROR: Reference file doesn't exist: ${ref_file}"
        exit 1
    }
    
    println "Reference file exists: ${ref_file}"
    
    // Create a proper index channel
    ch_fasta_index = Channel
        .fromPath("${ref_file}.{amb,ann,bwt,fai,pac,sa}")
        .ifEmpty { 
            println "WARNING: No index files found for: ${ref_file}" 
            println "Continuing with test file"
            return file(workflow.projectDir + "/nextflow.config")
        }
        .collect()
    
    ch_fasta_index.view { println "Index files: ${it}" }

    BWA_MEM(
        ch_reads,
        ref_file,
        ch_fasta_index
    )
    
    BWA_MEM.out.log.view { println "Debug log: ${it.text}" }
}
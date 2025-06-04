nextflow.enable.dsl = 2

// Parameters
// params.input = "/mnt/f/projects/250224_DFSP_WES/data/rna/sample_info/samplesheet.csv"
params.input = "/mnt/f/projects/250224_DFSP_WES/data/rna/sample_info/test.csv"
params.ref_dir = "/mnt/m/Reference"
params.star_index = "${params.ref_dir}/Gencode/STAR_index/"
params.ctat_resource_lib = "${params.ref_dir}/GRCh38_gencode_v33_CTAT_lib_Apr062020.plug-n-play/ctat_genome_lib_build_dir/"
// params.outdir = "/mnt/f/projects/250224_DFSP_WES/outputs/stump/STAR-Fusion"

// Load required processes
include { STAR_ALIGNMENT } from '../modules/local/starfusion/alignment/main.nf'
include { STAR_FUSION } from '../modules/local/starfusion/fusion/main.nf'

workflow {

    star_index = file(params.star_index)
    ctat_resource_lib = file(params.ctat_resource_lib)

    // Read input CSV
    sample_list = Channel.fromPath(params.input)
        .ifEmpty { exit(1, "Samplesheet not found: ${params.input}") }
        .splitCsv(header: true, quote: '"', strip: true)
    
    input_samples = sample_list
        .map{
            row -> 
            // Helper function to clean csv
            def cleanValue = { value ->
                if (value == null || value == 'NA' || value == 'NULL' || value == '') {
                    return null
                }
                // Remove quotes and trim whitespace
                return value.toString().replaceAll('^"(.*)"$', '$1').trim()
            }

            // Extract sample info
            def sample_id = cleanValue(row.sample)

            // Clean and validate file paths
            def fastq_1_clean = cleanValue(row.fastq_1)
            def fastq_2_clean = cleanValue(row.fastq_2)
            def fastq_trimmed_1_clean = cleanValue(row.fastq_trimmed_1)
            def fastq_trimmed_2_clean = cleanValue(row.fastq_trimmed_2)
            
           // Determine available input types
            def has_fastq = fastq_1_clean && fastq_2_clean
            def has_fastq_trimmed = fastq_trimmed_1_clean && fastq_trimmed_2_clean
            
            // Exit if no input is provided
            if (!has_fastq && !has_fastq_trimmed) {
                error("Either FASTQ or Trimmed FASTQ must be provided for sample ${sample_id}")
            }

            def meta = [
                id: sample_id,
                sample_id: sample_id
            ]

            // Initialize all file variables with null
            def fastq_1 = null
            def fastq_2 = null
            def fastq_trimmed_1 = null
            def fastq_trimmed_2 = null
            

            if (has_fastq) {
                // Validate FASTQ files exist
                try {
                    fastq_1 = file(fastq_1_clean, checkIfExists: true)
                    fastq_2 = file(fastq_2_clean, checkIfExists: true)
                } catch (Exception e) {
                    error("FASTQ file not found for sample ${sample_id}: ${e.message}")
                }
            }
            
            if (has_fastq_trimmed) {
                // Validate FASTQ files exist
                try {
                    fastq_trimmed_1 = file(fastq_trimmed_1_clean, checkIfExists: true)
                    fastq_trimmed_2 = file(fastq_trimmed_2_clean, checkIfExists: true)
                } catch (Exception e) {
                    error("FASTQ file not found for sample ${sample_id}: ${e.message}")
                }
            }                

            // Return unified format with all possible inputs
            return [meta, fastq_1, fastq_2, fastq_trimmed_1, fastq_trimmed_2]

        }

    fastq = input_samples
        .filter {
            _meta, fastq_1, fastq_2, _fastq_trimmed_1, _fastq_trimmed_2 -> 
            fastq_1 != null && fastq_2 != null
        }
        .map {
            meta, fastq_1, fastq_2, _fastq_trimmed_1, _fastq_trimmed_2 -> 
            [meta, fastq_1, fastq_2]
        }
    
    fastq_trimmed = input_samples
        .filter {
            _meta, _fastq_1, _fastq_2, fastq_trimmed_1, fastq_trimmed_2 -> 
            fastq_trimmed_1 != null && fastq_trimmed_2 != null
        }
        .map {
            meta, _fastq_1, _fastq_2, fastq_trimmed_1, fastq_trimmed_2 -> 
            [meta, fastq_trimmed_1, fastq_trimmed_2]
        }     

    fastq_trimmed.view()

    // Run STAR alignment
    STAR_ALIGNMENT(
        fastq_trimmed,
        star_index
    )

    // Run STAR-Fusion
    star_fusion_input = STAR_ALIGNMENT.out.junction
        .join(fastq_trimmed)
        .map { meta, junction, fq1, fq2 -> 
            tuple(meta, fq1, fq2, junction) 
        }

    STAR_FUSION(
        star_fusion_input,
        ctat_resource_lib
    )
}

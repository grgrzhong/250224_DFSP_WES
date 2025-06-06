workflow PREPARE_SAMPLE {
    take:
    input_csv  // Path to the input CSV file

    main:
    // Create unified channel from the CSV with support for both input types
    sample_list = Channel.fromPath(input_csv)
        .ifEmpty { exit(1, "Samplesheet not found: ${input_csv}") }
        .splitCsv(header: true, quote: '"', strip: true)
        // .view()
    
    // sample_list.view()

    input_samples = sample_list
        .map { 
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
            def patient_id = cleanValue(row.patient)
            def sample_id = cleanValue(row.sample)
            def status = cleanValue(row.status)

            // Validate required fields
            if (!patient_id) error("Missing or empty 'patient' field in row: ${row}")
            if (!sample_id) error("Missing or empty 'sample' field in row: ${row}")
            if (!status) error("Missing or empty 'status' field in row: ${row}")

            // Clean and validate file paths
            def fastq_1_clean = cleanValue(row.fastq_1)
            def fastq_2_clean = cleanValue(row.fastq_2)
            def bam_clean = cleanValue(row.bam)
            def bai_clean = cleanValue(row.bai)

            // Determine available input types
            def has_fastq = fastq_1_clean && fastq_2_clean
            def has_bam = bam_clean

            // Exit if no input is provided
            if (!has_fastq && !has_bam) {
                error("Either FASTQ or BAM files must be provided for sample ${sample_id}")
            }
            
            // Convert status to integer with validation
            def status_int
            try {
                status_int = status.toInteger()
                if (status_int != 0 && status_int != 1) {
                    error("Status must be 0 (normal) or 1 (tumour) for sample ${sample_id}, got: ${status}")
                }
            } catch (NumberFormatException _e) {
                error("Invalid status value '${status}' for sample ${sample_id}. Must be 0 or 1")
            }
            
            def meta = [
                patient_id: patient_id,
                sample_id: sample_id,
                status: status_int
            ]
            
            // Initialize all file variables with null
            def fastq_1 = null
            def fastq_2 = null
            def bam = null
            def bai = null

            if (has_fastq) {
                // Validate FASTQ files exist
                try {
                    fastq_1 = file(fastq_1_clean, checkIfExists: true)
                    fastq_2 = file(fastq_2_clean, checkIfExists: true)
                } catch (Exception e) {
                    error("FASTQ file not found for sample ${sample_id}: ${e.message}")
                }
            }

            if (has_bam) {
                // Validate BAM file exists
                try {
                    bam = file(bam_clean, checkIfExists: true)
                    // BAI is optional, but if provided, check it exists
                    if (bai_clean) {
                        bai = file(bai_clean, checkIfExists: true)
                    } else {
                        // Try to find BAI file automatically
                        def auto_bai = file("${bam_clean}.bai")
                        if (auto_bai.exists()) {
                            bai = auto_bai
                            log.info("Found BAI file automatically for ${sample_id}: ${auto_bai}")
                        } else {
                            log.warn("No BAI file found for ${sample_id}. BAM indexing may be required.")
                        }
                    }
                } catch (Exception e) {
                    error("BAM file not found for sample ${sample_id}: ${e.message}")
                }
            }                

            // Return unified format with all possible inputs
            return [meta, fastq_1, fastq_2, bam, bai]
    
        }

    // Prepare the fastq channel for preprocessing and mapping
    fastq = input_samples
        .filter { 
            _meta, fastq_1, fastq_2, _bam, _bai ->
            fastq_1 != null && fastq_2 != null
        }
        .map { 
            meta, fastq_1, fastq_2, _bam, _bai ->
            [meta, fastq_1, fastq_2]
        }
        // .view()
    
    
    tumour_samples = input_samples
        .filter { 
            meta, _fastq_1, _fastq_2, _bam, _bai ->
            meta.status == 1 
        }
        // .view()

    normal_samples = input_samples
        .filter { 
            meta, _fastq_1, _fastq_2, _bam, _bai ->
            meta.status == 0 
        }
        // .view()


    // Prepare bam input for tumour_normal paired samples
    bam_tumour_normal = tumour_samples
        .map {
            meta, fastq_1, fastq_2, bam, bai -> 
            [meta.patient_id, [meta, fastq_1, fastq_2, bam, bai]] 
        }
        .combine(
            normal_samples
                .map { 
                    meta, fastq_1, fastq_2, bam, bai -> 
                    [meta.patient_id, [meta, fastq_1, fastq_2, bam, bai]] 
                },
            by: 0
        )
        .map { 
            patient_id, tumour, normal ->

            def tumour_meta = tumour[0]
            def normal_meta = normal[0]

            def meta = [
                id: tumour_meta.sample_id,
                patient_id: patient_id,
                tumour_id: tumour_meta.sample_id,
                normal_id: normal_meta.sample_id,
                is_paired: true
            ]
            
            return [meta, tumour[3], tumour[4], normal[3], normal[4]]
        }
        .ifEmpty {
            Channel.empty()
        }
        // .view()

    // Create a channel of patient IDs that have normal samples
    normal_patient_ids = normal_samples
        .map { meta, _fastq_1, _fastq_2, _bam, _bai -> meta.patient_id }
        .unique()
        .collect()
        .ifEmpty([]) // Empty list if no normal samples

    // Create bam tumour-only samples channel
    bam_tumour_only = tumour_samples
        .map { meta, fastq_1, fastq_2, bam, bai -> [meta.patient_id, [meta, fastq_1, fastq_2, bam, bai]] }
        .combine(normal_patient_ids)
        .filter { patient_id, _tumour_data, normal_ids -> !(patient_id in normal_ids) }
        .map { patient_id, tumour_data, _normal_ids ->
            def tumour_meta = tumour_data[0]
            
            def meta = [
                id: tumour_meta.sample_id,
                patient_id: patient_id,
                tumour_id: tumour_meta.sample_id,
                normal_id: null,
                is_paired: false
            ]
            
            return [meta, tumour_data[3], tumour_data[4]] // [meta, bam, bai]
        }
        .ifEmpty {
            Channel.empty()
        }
        .view()

    // Emit all channels
    emit:
        input_samples       = input_samples
        tumour_samples      = tumour_samples
        normal_samples      = normal_samples
        fastq               = fastq
        bam_tumour_normal   = bam_tumour_normal
        bam_tumour_only     = bam_tumour_only
}
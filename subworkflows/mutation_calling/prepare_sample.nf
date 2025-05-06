workflow PREPARE_SAMPLE {
    take:
        input_csv  // Path to the input CSV file

    main:
        
        // Create a channel for fastq files
        reads = Channel.fromPath(input_csv)
            .ifEmpty { exit(1, "Sample sheet not found at ${input_csv}") }
            .splitCsv(header: true)
            .filter { row -> row.containsKey('fastq_1') && row.containsKey('fastq_2') }
            .map { row ->
                def meta = [:]
                meta.id = row.sample
                meta.patient_id = row.patient
                meta.status = row.status.toInteger()

                // Check that the fastq files exist
                def fastq_1 = file(row.fastq_1)
                def fastq_2 = file(row.fastq_2)

                if (!fastq_1.exists()) {
                    error("Read 1 fastq file not found: ${row.fastq_1}")
                }
                if (!fastq_2.exists()) {
                    error("Read 2 fastq file not found: ${row.fastq_2}")
                }

                return [meta, fastq_1, fastq_2]
            }

        // Create a channel for BAM files
        input_samples = Channel.fromPath(input_csv)
            .ifEmpty { exit(1, "Samplesheet not found: ${input_csv}") }
            .splitCsv(header: true)
            .map { row ->
                // Extract sample info
                def patient_id = row.patient ? row.patient.trim() : null
                def sample_id = row.sample ? row.sample.trim() : null
                def status = row.status ? row.status.trim() : null

                // Validate required fields
                if (!patient_id) {
                    error("Missing or empty 'patient' field in row: ${row}")
                }
                if (!sample_id) {
                    error("Missing or empty 'sample' field in row: ${row}")
                }
                if (!status) {
                    error("Missing or empty 'status' field in row: ${row}")
                }
                if (!row.bam) {
                    error("Missing 'bam' field in row: ${row}")
                }

                // Process BAM and BAI paths
                def bam = file(row.bam)
                def bai = row.bai ? file(row.bai) : file("${row.bam}.bai")

                // Check if files exist
                if (!bam.exists()) {
                    error("BAM file not found: ${bam}")
                }
                if (!bai.exists()) {
                    error("BAI file not found: ${bai}")
                }

                // Return a tuple with patient_id, sample_id, status, bam, and bai
                [
                    patient_id: patient_id,
                    sample_id: sample_id,
                    status: status.toInteger(),
                    bam: bam,
                    bai: bai,
                ]
            }

        // Split samples into tumour and normal
        tumour_samples = input_samples.filter { it.status == 1 }
        normal_samples = input_samples.filter { it.status == 0 }

        // Create paired samples channel
        paired_samples = tumour_samples
            .map { tumour -> [ tumour.patient_id, tumour] }
            .combine(
                normal_samples.map { normal -> [normal.patient_id, normal] },
                by: 0
            )
            .map {
                patient_id, tumour, normal -> 
                def meta = [
                    id: "${tumour.sample_id}_vs_${normal.sample_id}",
                    patient_id: patient_id,
                    tumour_id: tumour.sample_id,
                    normal_id: normal.sample_id,
                    is_paired: true
                ]
                [meta, tumour.bam, tumour.bai, normal.bam, normal.bai]
            }
        
        // Create tumour-only samples channel
        normal_patient_ids = normal_samples
            .map { it.patient_id }
            .unique()
            .collect()
            .map { ids -> ids.sort() }

        unpaired_samples = tumour_samples
            .branch {
                def normal_ids = normal_patient_ids.val
                paired: normal_ids.contains(it.patient_id)
                unpaired: true
            }
            .unpaired
            .map {
                tumour -> def meta = [
                    id: tumour.sample_id,
                    patient_id: tumour.patient_id,
                    tumour_id: tumour.sample_id,
                    normal_id: null,
                    is_paired: false
                ]
                [meta, tumour.bam, tumour.bai, null, null]
            }
        
        // Combine paired and unpaired samples
        all_samples = paired_samples.mix(unpaired_samples)
        
        // Extract tumor samples from paired samples for various analyses
        paired_tumour_samples = paired_samples
            .map { meta, tumour_bam, tumour_bai, _normal_bam, _normal_bai ->
                [meta, tumour_bam, tumour_bai]
            }

        // Extract normal samples from paired samples
        paired_normal_samples = paired_samples
            .map { meta, _tumour_bam, _tumour_bai, normal_bam, normal_bai ->
                def normal_meta = meta.clone()
                normal_meta.tumour_id = normal_meta.normal_id
                [normal_meta, normal_bam, normal_bai]
            }

        // Extract tumor samples for unpaired analysis
        unpaired_tumour_samples = unpaired_samples
            .map { meta, tumor_bam, tumor_bai, _normal_bam, _normal_bai ->
                [meta, tumor_bam, tumor_bai]
            }

        // Log sample counts
        paired_samples.count().subscribe { count ->
            log.info("Found ${count} paired samples.")
        }
        
        unpaired_samples.count().subscribe { count ->
            log.info("Found ${count} unpaired samples.")
        }

    emit:
        reads                   = reads
        all_samples             = all_samples
        tumour_samples          = tumour_samples
        normal_samples          = normal_samples
        paired_samples          = paired_samples
        unpaired_samples        = unpaired_samples
        paired_tumour_samples   = paired_tumour_samples
        paired_normal_samples   = paired_normal_samples
        unpaired_tumour_samples = unpaired_tumour_samples
}
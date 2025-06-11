nextflow.enable.dsl = 2

// Parameters
params.input_vcf = "/home/zhonggr/projects/250224_DFSP_WES/data/wes/mutect2/DFSP-001-T/DFSP-001-T.final.vcf.gz"
params.sample_name = "DFSP-001-T"
params.ref_dir = "/mnt/m/Reference"
params.vep_dir = "${params.ref_dir}/VEP_cache"
params.pcgr_ref_dir = "${params.ref_dir}/PCGR_reference/20250314"
params.panel_of_normal = "/mnt/m/WES/DFSP/PON-Mutect/pon.vcf.gz"
params.outdir = "results"

include { FORMAT_VCF_TUMOUR_NORMAL } from '../modules/local/pcgr/format_vcf/tumour_normal/main.nf'
include { FORMAT_VCF_TUMOUR_ONLY   }  from '../modules/local/pcgr/format_vcf/tumour_only/main.nf'

workflow  {

    vcf= file(params.input_vcf)
    vcf_index = vcf + ".tbi"
    sample_name = params.sample_name

    input_ch = Channel.of(
        [
            [id: sample_name, normal_id: "DFSP-001-N"],
            vcf,
            vcf_index
        ]
    )

    FORMAT_VCF_TUMOUR_NORMAL(input_ch)
}
#!/bin/bash

python /home/zhonggr/projects/250224_DFSP_WES/data/reference/InterVar-2.2.1/Intervar.py \
    -b hg19 \
    -i /home/zhonggr/projects/250224_DFSP_WES/test/DFSP-001-T/DFSP-001-T.final.vcf.gz \
    --input_type=VCF \
    -o /home/zhonggr/projects/250224_DFSP_WES/test/DFSP-010-T-P1/DFSP-010-T-P1 \
    --table_annovar /home/zhonggr/projects/250224_DFSP_WES/data/annovar/table_annovar.pl \
    --convert2annovar /home/zhonggr/projects/250224_DFSP_WES/data/annovar/convert2annovar.pl \
    --annotate_variation /home/zhonggr/projects/250224_DFSP_WES/data/annovar/annotate_variation.pl \
    -d /home/zhonggr/projects/250224_DFSP_WES/data/reference/annovar/humandb/ \
    -t /home/zhonggr/projects/250224_DFSP_WES/data/reference/InterVar-2.2.1/intervardb
    # --skip_annovar

python /home/zhonggr/projects/250224_DFSP_WES/data/reference/InterVar-2.2.1/Intervar.py \
    -b hg38 \
    -i /home/zhonggr/projects/250224_DFSP_WES/test/DFSP-031-T/DFSP-031-T.hg38_multianno.txt \
    --input_type=AVinput  \
    -o /home/zhonggr/projects/250224_DFSP_WES/test/DFSP-031-T/DFSP-031-T \
    --table_annovar /home/zhonggr/projects/250224_DFSP_WES/data/annovar/table_annovar.pl \
    --convert2annovar /home/zhonggr/projects/250224_DFSP_WES/data/annovar/convert2annovar.pl \
    --annotate_variation /home/zhonggr/projects/250224_DFSP_WES/data/annovar/annotate_variation.pl \
    -d /home/zhonggr/projects/250224_DFSP_WES/data/reference/annovar/humandb/ \
    -t /home/zhonggr/projects/250224_DFSP_WES/data/reference/InterVar-2.2.1/intervardb \
    --skip_annovar
    
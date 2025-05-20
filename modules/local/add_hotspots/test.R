library(tidyverse)
library(readxl)

# Read ANNOVAR output (adjust separator if needed)
annovar_file <- "/home/zhonggr/projects/250224_DFSP_WES/data/wes/variant_annotation/annovar/DFSP-001-T/DFSP-001-T_annovar.txt"

annovar_data <- read_tsv(
    annovar_file, 
    col_names = TRUE, 
    comment = "#"
)



# Read OncoKB (TSV)
oncokb <- read_tsv("allAnnotatedVariants.tsv")

# Read Cancer Hotspots (Excel)
hotspots <- read_excel("hotspots_v2.xls", sheet = 1)
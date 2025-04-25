library(maftools)
library(qs)
library(fs)
library(here)
library(tidyverse)
library(readxl)

## Read the annovar output
input_dir <- here("data/wes/variant_calling/mutect2")
input_files <- dir_ls(input_dir, recurse = TRUE, glob = "*annovar.txt")

read_tsv(input_files[1])

maf_data <- annovarToMaf(
  annovar = input_files,
  refBuild = "hg38",
  tsbCol = "sample_id",
  table = "ensGene",
  MAFobj = FALSE
)

maf_data |> as_tibble()

## Read the funcotator output
input_dir <- here("data/wes/variant_calling/mutect2")
input_files <- dir_ls(input_dir, recurse = TRUE, glob = "*annotated.tsv")

input_data <- map2(
  input_files, names(input_files),
  ~ {
    maf_data <- read_tsv(.x, show_col_types = FALSE)
    maf_data$Tumor_Sample_Barcode <- basename(.y) |> str_remove("\\..*$")

    return(maf_data)
  }
)

maf_data <- merge_mafs(mafs = input_data, verbose = TRUE)
qsave(maf_data, here("data/wes/annotation/summary/merged_maf.qs"))


maf_data <- qread(here("data/wes/annotation/summary/merged_maf.qs"))
str(maf_data)

sample_info <- read_excel(
  here("data/clinical/DFSP-Multiomics-Sample list (updated 2024.09).xlsx")
) |> 
  select(Sample.ID, Specimen.Class, Histology.Nature)

unique(sample_info$Histology.Nature)



group_info <- read_excel(here("data/clinical/groups.xlsx"))

sample_info |> 
  filter(Histology.Nature %in% group_info$group1)

# getSampleSummary(maf_data)
# getGeneSummary(maf_data)
# getClinicalData(maf_data)
# getFields(maf_data)
# write.mafSummary(maf = maf_data, basename = "maf_data")

plotmafSummary(
  maf = maf_data,
  rmOutlier = TRUE,
  addStat = "median",
  dashboard = TRUE,
  titvRaw = FALSE
)

oncoplot(maf = maf_data, top = 20)

maf_data.titv <- titv(maf = maf_data, plot = FALSE, useSyn = TRUE)
# plot titv summary
plotTiTv(res = maf_data.titv)

maf_data@data |>
  as_tibble() |>
  arrange(Hugo_Symbol) |>
  filter(
    Hugo_Symbol %in% c("MYOD1", "SUZ12", "NF1")
  )

lollipopPlot(
  maf = maf_data,
  gene = "MYOD1",
  AACol = "Protein_Change",
  showMutationRate = TRUE,
  labelPos = "all"
)

lollipopPlot(
  maf = maf_data,
  gene = "SUZ12",
  AACol = "Protein_Change",
  showMutationRate = TRUE,
  labelPos = "all"
)

maf_data@data |>
  as_tibble() |>
  arrange(Hugo_Symbol) |>
  filter(
    # Variant_Classification %in% c("Splice_Site"),
    Hugo_Symbol %in% c("MYOD1", "SUZ12", "NF1")
  ) |>
  # view()
  select(
    Hugo_Symbol, Tumor_Sample_Barcode,
    Genome_Change, cDNA_Change, Protein_Change, Codon_Change
  )

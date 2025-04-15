library(maftools)
library(qs)
library(fs)
library(here)
library(tidyverse)
library(readxl)

## Read the annovar output
input_dir <- here("data/SARC/annotation/annovar")
input_files <- dir_ls(input_dir, recurse = TRUE, glob = "*.txt")

read_tsv(input_files[1])

maf_data <- annovarToMaf(
  annovar = input_files, 
  refBuild = 'hg38',
  tsbCol = 'sample_id', 
  table = 'ensGene',
  MAFobj = FALSE
)

maf_data |> as_tibble()

## Read the funcotator output
input_dir <- here("data/SARC/variant_calling/funcotator")
input_files <- dir_ls(input_dir, recurse = TRUE, glob = "*.tsv")

input_data <- map2(
  input_files, names(input_files),
  ~ {
    maf_data <- read_tsv(.x, show_col_types = FALSE)
    maf_data$Tumor_Sample_Barcode <- basename(.y) |> str_remove("\\..*$")
  
    return(maf_data)
  }
)

maf_data <- merge_mafs(mafs = input_data, verbose = TRUE)
str(maf_data)

getSampleSummary(maf_data)
getGeneSummary(maf_data)
getClinicalData(maf_data)
getFields(maf_data)
write.mafSummary(maf = maf_data, basename = 'maf_data')

plotmafSummary(
  maf = maf_data, 
  rmOutlier = TRUE, 
  addStat = 'median', 
  dashboard = TRUE, 
  titvRaw = FALSE
)

oncoplot(maf = maf_data, top = 10)

maf_data.titv = titv(maf = maf_data, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = maf_data.titv)

maf_data@data |> 
  as_tibble() |>
  arrange(Hugo_Symbol) |> 
  filter(
    Hugo_Symbol %in% c("MYOD1", "SUZ12", "NF1")
  )

lollipopPlot(
  maf = maf_data,
  gene = 'MYOD1',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  labelPos = "all"
)

lollipopPlot(
  maf = maf_data,
  gene = 'SUZ12',
  AACol = 'Protein_Change',
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
    Hugo_Symbol,Tumor_Sample_Barcode, 
    Genome_Change, cDNA_Change, Protein_Change, Codon_Change
  )
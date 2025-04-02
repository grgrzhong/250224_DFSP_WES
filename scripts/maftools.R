library(maftools)
library(qs)
library(here)
library(tidyverse)
library(readxl)

# ann_dir <- here("data/WES/DFSP/variant_annotation/annovar/DFSP-001-T")
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

# read_tsv(laml.clin)
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
getSampleSummary(laml)
getGeneSummary(laml)
getClinicalData(laml)
getFields(laml)
write.mafSummary(maf = laml, basename = 'laml')

plotmafSummary(
  maf = laml, 
  rmOutlier = TRUE, 
  addStat = 'median', 
  dashboard = TRUE, 
  titvRaw = FALSE
)
oncoplot(maf = laml, top = 10)

laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)

lollipopPlot(
  maf = laml,
  gene = 'DNMT3A',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  labelPos = 882
)

my_data = data.frame(pos = sample.int(912, 15, replace = TRUE), count = sample.int(30, 15, replace = TRUE))
head(my_data)

lollipopPlot(data = my_data, gene = "DNMT3A")
plotProtein(gene = "TP53", refSeqID = "NM_000546")

brca <- system.file("extdata", "brca.maf.gz", package = "maftools")
brca = read.maf(maf = brca, verbose = FALSE)

rainfallPlot(maf = brca, detectChangePoints = TRUE, pointSize = 0.4)

laml.mutload = tcgaCompare(maf = laml, cohortName = 'Example-LAML', logscale = TRUE, capture_size = 50)

plotVaf(maf = laml, vafCol = 'i_TumorVAF_WU')

gistic_res_folder <- system.file("extdata", package = "maftools")
laml.gistic = readGistic(gisticDir = gistic_res_folder, isTCGA = TRUE)

gisticChromPlot(gistic = laml.gistic, markBands = "all")

coGisticChromPlot(gistic1 = laml.gistic, gistic2 = laml.gistic, g1Name = "AML-1", g2Name = "AML-2", type = 'Amp')

gisticBubblePlot(gistic = laml.gistic)

gisticOncoPlot(gistic = laml.gistic, clinicalData = getClinicalData(x = laml), clinicalFeatures = 'FAB_classification', sortByAnnotation = TRUE, top = 10)

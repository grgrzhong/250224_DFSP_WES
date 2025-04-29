library(maftools)
library(tidyverse)

dfsp.maf <- annovarToMaf(list.files("Y:/WES/DFSP - Pilot/Annovar", full.names = TRUE))
dfsp.clin.path <- "G:/Sarcoma/GLI1/GLI1-Metadata.tsv"
dfsp.clin<- read.delim("G:/Sarcoma/GLI1/GLI1-Metadata.tsv", sep = "\t")
sample.classic <- dfsp.clin$Tumor_Sample_Barcode[dfsp.clin$Histology.Subtype=="Classic"]
sample.pdgfb <- dfsp.clin$Tumor_Sample_Barcode[dfsp.clin$Histology.Subtype=="Classic" & dfsp.clin$Molecular.Subtype=="PDGFB"]
sample.pdgfd <- dfsp.clin$Tumor_Sample_Barcode[dfsp.clin$Histology.Subtype=="Classic" & dfsp.clin$Molecular.Subtype=="PDGFD"]
sample.fst <- dfsp.clin$Tumor_Sample_Barcode[dfsp.clin$Histology.Subtype=="FST"]
sample.met <- dfsp.clin$Tumor_Sample_Barcode[dfsp.clin$Origin=="Metastasis"]


dfsp.maf.filtered <- dfsp.maf %>% separate(AD, into = c("RAD", "VAD"), sep = ",", remove = TRUE)

dfsp.maf.filtered <- dfsp.maf.filtered %>% filter(ExonicFunc.refGene!= "synonymous SNV" & 
                                                    Func.refGene %in% c("exonic", "splicing") &
                                                    AF>=0.05 & VAD >=5  & 
                                                    DP >=10 & gnomAD_exome_ALL <0.001)

cancerGeneList <- read.delim("G:/Sarcoma/cancerGeneList.tsv", sep = "\t")
ActionableGene <- read.delim("G:/Sarcoma/oncokb_biomarker_drug_associations.tsv", sep = "\t")
ActionableGene.CNV <- ActionableGene[ActionableGene$Alterations %in% c("Amplification", "Deletion"),]

dfsp.maf.filtered.hotspot <- dfsp.maf.filtered[dfsp.maf.filtered$Hugo_Symbol %in% cancerGeneList$Hugo.Symbol,]

CNV <- read.delim("G:/Sarcoma/GLI1/CNV-GATK/GLI-004-T.called.annotated.bed", sep = "\t")
CNV_filtered <- CNV[CNV$X. != "0",]
CNV_filtered <- CNV_filtered[CNV_filtered$KLHL17 %in% ActionableGene.CNV$Gene,]

dfsp.maf.classic <- annovarToMaf(
  list.files("G:/Sarcoma/DFSP/Annovar", full.names = TRUE) [unlist(strsplit(list.files("G:/Sarcoma/DFSP/Annovar"), split = ".txt")) %in% sample.classic])
dfsp.maf.fst <- annovarToMaf(
  list.files("G:/Sarcoma/DFSP/Annovar", full.names = TRUE) [unlist(strsplit(list.files("G:/Sarcoma/DFSP/Annovar"), split = ".txt")) %in% sample.fst])
dfsp.maf.pdgfb <- annovarToMaf(
  list.files("G:/Sarcoma/DFSP/Annovar", full.names = TRUE) [unlist(strsplit(list.files("G:/Sarcoma/DFSP/Annovar"), split = ".txt")) %in% sample.pdgfb])
dfsp.maf.pdgfd <- annovarToMaf(
  list.files("G:/Sarcoma/DFSP/Annovar", full.names = TRUE) [unlist(strsplit(list.files("G:/Sarcoma/DFSP/Annovar"), split = ".txt")) %in% sample.pdgfd])
dfsp.maf.met <- annovarToMaf(
  list.files("G:/Sarcoma/DFSP/Annovar", full.names = TRUE) [unlist(strsplit(list.files("G:/Sarcoma/DFSP/Annovar"), split = ".txt")) %in% sample.met])



all.lesions <- "D:/Sarcoma/Methylation_profiling/DFSP/CNV/All/all_lesions.conf_99.txt"
amp.genes <-"D:/Sarcoma/Methylation_profiling/DFSP/CNV/amp_genes.conf_99.txt"
del.genes <- "D:/Sarcoma/Methylation_profiling/DFSP/CNV/del_genes.conf_99.txt"
scores.gis <- "D:/Sarcoma/Methylation_profiling/DFSP/CNV/scores.gistic"


#**GisticOncoplot from EPIC CNV**#

SampleSheet_DFSP <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_97pt_Meth_immune_CNV.csv", row.names = 1)

SampleSheet_DFSP <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_all.csv")
SampleSheet_DFSP <- SampleSheet_DFSP %>% filter(Main == "Yes",
                                                !(Specimen.Nature %in% c("Metastasis", "Recurrence (Post-imatinib)")),
                                                Specimen.Class == "Tumour"
                                                )
colnames(SampleSheet_DFSP)[2] <- "Tumor_Sample_Barcode"
SampleSheet_DFSP_Meth <- read.csv("D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_96pt_Meth_immune.csv")


colnames(SampleSheet_DFSP)[1] <- "Sample.ID"
colnames(SampleSheet_DFSP)[1] <- "Tumor_Sample_Barcode"

gistic.dfsp = readGistic(gisticDir = "D:/Sarcoma/Methylation_profiling/DFSP/CNV/GISTIC/All97pts")
gistic.dfsp.pdgfb <- readGistic(gisticDir = "D:/Sarcoma/Methylation_profiling/DFSP/CNV/GISTIC/PDGFB")
gistic.dfsp.pdgfd <- readGistic(gisticDir = "D:/Sarcoma/Methylation_profiling/DFSP/CNV/GISTIC/PDGFD")
gistic.dfsp.pairedClassic <- readGistic(gisticDir = "D:/Sarcoma/Methylation_profiling/DFSP/CNV/GISTIC/PairedClassic")
gistic.dfsp.pairedFST <- readGistic(gisticDir = "D:/Sarcoma/Methylation_profiling/DFSP/CNV/GISTIC/PairedFST")
gistic.dfsp.AllFST <- readGistic(gisticDir = "D:/Sarcoma/Methylation_profiling/DFSP/CNV/GISTIC/AllFST")
gistic.dfsp.ClassicOnly <- readGistic(gisticDir = "D:/Sarcoma/Methylation_profiling/DFSP/CNV/GISTIC/ClassicOnly")
gistic.dfsp.Metastasis <- readGistic(gisticDir = "D:/Sarcoma/Methylation_profiling/DFSP/CNV/GISTIC/Metastasis")

gisticChromPlot(gistic = gistic.dfsp.pdgfb, markBands = "all", fdrCutOff = 0.25,
                ref.build = "hg19")

coGisticChromPlot(gistic1 = gistic.dfsp.ClassicOnly, gistic2 = gistic.dfsp.Metastasis, 
                  g1Name = "Classic", g2Name = "Metastasis", type = 'Del',
                  txtSize = 0.8, rugTickSize = 0.2, mutGenesTxtSize = 2.6,
                  fdrCutOff = 0.25, ref.build = "hg19")
coGisticChromPlot(gistic1 = gistic.dfsp.pdgfb, gistic2 = gistic.dfsp.pdgfd, 
                  g1Name = "PDGFB", g2Name = "PDGFD", type = 'Del', 
                  fdrCutOff = 0.05)

features <- c("CN.subgroup.WardD2", "Meth.subgroup.probesel", "Molecular.subtype", 
              "Histology.Subtype", "Sex" ,"Age.Cat",
              "Site.Category.Broad", "Local.recurrence","Metastasis")

phenoData$Tumor_Sample_Barcode <- phenoData$Sample.ID

gisticOncoPlot(gistic = gistic.dfsp, clinicalData = phenoData, 
               clinicalFeatures = features, sortByAnnotation = TRUE, top = 30,
               gene_mar	= 10, barcode_mar = 10, sepwd_genes = 0.5, sepwd_samples = 0.25,
               annotationFontSize = 1.2, legendFontSize = 1.2, fontSize = 0.3,
               showTumorSampleBarcodes = FALSE, removeNonAltered = FALSE,
               colors = colors)

colors<- c("red","blue")
names(colors) <- c("Amp","Del")
gisticOncoPlot(gistic = gistic.dfsp.pdgfd, removeNonAltered = FALSE)

cnv.df <- data.frame(t(gistic.dfsp@cnMatrix))
cnv.df$Sample.ID <- rownames(cnv.df)
SampleSheet_DFSP_CNV <- merge(SampleSheet_DFSP, cnv.df, by = "Tumor_Sample_Barcode", all.x = TRUE)
write.csv(SampleSheet_DFSP_CNV, "D:/Sarcoma/Methylation_profiling/DFSP/phenoData_DFSP_97pt_Meth_immune_CNV.csv")

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
#GSEA analysis requires a ranked gene list, which contains three features:
#numeric vector: fold change or other type of numerical variable
#named vector: every number has a name, the corresponding gene ID
#sorted vector: number should be sorted in decreasing order

d = data.frame(gistic.dfsp.pdgfb@gene.summary)
## assume 1st column is ID
## 2nd column is FC or other type of numerical variable

## feature 1: numeric vector
geneList = d[,2]

## feature 2: named vector
names(geneList) = as.character(d[,1])

## feature 3: decreasing order
geneList = sort(geneList, decreasing = TRUE)

gene_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = names(geneList), columns = "ENTREZID", keytype = "SYMBOL")
names(geneList) <- gene_ids[,2]
geneList <- geneList[!is.na(names(geneList))]

ego.cc <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              keyType      = "ENTREZID",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = TRUE)

ego.bp <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                keyType      = "ENTREZID",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = TRUE)

ego.mf <- gseGO(geneList     = geneList,
                OrgDb        = org.Hs.eg.db,
                ont          = "MF",
                keyType      = "SYMBOL",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = TRUE)

dotplot(ego.cc, showCategory=10) + ggtitle("dotplot for GSEA")

#Survival analysis
library(survival)
library(survminer)

# Convert survival months into numeric variables
SampleSheet_DFSP_CNV$OS.time <- as.numeric(SampleSheet_DFSP_CNV$OS.time)
SampleSheet_DFSP_CNV$RFS.time <- as.numeric(SampleSheet_DFSP_CNV$RFS.time)

# Produce survival object
survival_obj <- Surv(SampleSheet_DFSP_CNV$OS.time, SampleSheet_DFSP_CNV$OS.status == "1:DECEASED")
survival_obj <- Surv(SampleSheet_DFSP_CNV$RFS.time, SampleSheet_DFSP_CNV$RFS.status== "1:Recurrence")

km_curves <- survfit(survival_obj ~ AP_4.11q24.1, data = SampleSheet_DFSP_CNV)
DP_12.22q13.33 AP_4.11q24.1
SampleSheet_DFSP_CNV$AP_4.11q24.1

surv_plot <- ggsurvplot(km_curves, data = SampleSheet_DFSP_CNV, conf.int = FALSE, pval = TRUE,
                        pval.method = TRUE, risk.table = TRUE)
surv_plot


fisher.test(table(SampleSheet_DFSP_CNV$AP_4.11q24.1, SampleSheet_DFSP_CNV$FST))

ddfsp <- read.maf(maf = dfsp.maf.filtered)

, 
                 gisticAllLesionsFile = all.lesions,
                 gisticAmpGenesFile = amp.genes,
                 gisticDelGenesFile = del.genes,
                 gisticScoresFile = scores.gis,
                 clinicalData = dfsp.clin.path)

dfsp2 <- read.maf(maf = dfsp.maf, 
                 gisticAllLesionsFile = all.lesions.fst,
                 gisticAmpGenesFile = amp.genes.fst,
                 gisticDelGenesFile = del.genes.fst,
                 gisticScoresFile = scores.gis.fst,
                 clinicalData = dfsp.clin.path)

dfsp.gistic <- readGistic(gisticAllLesionsFile = all.lesions,
                          gisticAmpGenesFile = amp.genes,
                          gisticDelGenesFile = del.genes,
                          gisticScoresFile = scores.gis)
dfsp.classic <- read.maf(maf = dfsp.maf.classic, 
                         gisticAllLesionsFile = all.lesions.classic,
                         gisticAmpGenesFile = amp.genes.classic,
                         gisticDelGenesFile = del.genes.classic,
                         gisticScoresFile = scores.gis.classic,
                         clinicalData = dfsp.clin.path)

dfsp.classic <- read.maf(maf = dfsp.maf.classic, 
                         clinicalData = dfsp.clin.path)
dfsp.fst <- read.maf(maf = dfsp.maf.fst, 
                     clinicalData = dfsp.clin.path)
dfsp.paired <- read.maf(maf = dfsp.maf.paired, 
                        clinicalData = dfsp.clin.path)

dfsp.fst <- read.maf(maf = dfsp.maf.fst, 
                     gisticAllLesionsFile = all.lesions.fst,
                     gisticAmpGenesFile = amp.genes.fst,
                     gisticDelGenesFile = del.genes.fst,
                     gisticScoresFile = scores.gis.fst,
                     clinicalData = dfsp.clin.path)

dfsp.pdgfb <- read.maf(maf = dfsp.maf.pdgfb, 
                     gisticAllLesionsFile = all.lesions.pdgfb,
                     gisticAmpGenesFile = amp.genes.pdgfb,
                     gisticDelGenesFile = del.genes.pdgfb,
                     gisticScoresFile = scores.gis.pdgfb,
                     clinicalData = dfsp.clin.path)

dfsp.pdgfd <- read.maf(maf = dfsp.maf.pdgfd, 
                       gisticAllLesionsFile = all.lesions.pdgfd,
                       gisticAmpGenesFile = amp.genes.pdgfd,
                       gisticDelGenesFile = del.genes.pdgfd,
                       gisticScoresFile = scores.gis.pdgfd,
                       clinicalData = dfsp.clin.path)

dfsp.gistic <- readGistic(gisticAllLesionsFile = all.lesions,
                                gisticAmpGenesFile = amp.genes,
                                gisticDelGenesFile = del.genes,
                                gisticScoresFile = scores.gis)

dfsp.gistic.pdgfb <- readGistic(gisticAllLesionsFile = all.lesions.pdgfb,
                               gisticAmpGenesFile = amp.genes.pdgfb,
                               gisticDelGenesFile = del.genes.pdgfb,
                               gisticScoresFile = scores.gis.pdgfb)

dfsp.gistic.pdgfd <- readGistic(gisticAllLesionsFile = all.lesions.pdgfd,
                                gisticAmpGenesFile = amp.genes.pdgfd,
                                gisticDelGenesFile = del.genes.pdgfd,
                                gisticScoresFile = scores.gis.pdgfd)

dfsp.gistic.fst <- readGistic(gisticAllLesionsFile = all.lesions.fst,
                              gisticAmpGenesFile = amp.genes.fst,
                              gisticDelGenesFile = del.genes.fst,
                              gisticScoresFile = scores.gis.fst)

dfsp.gistic.classic <- readGistic( gisticAllLesionsFile = all.lesions.classic,
                                   gisticAmpGenesFile = amp.genes.classic,
                                   gisticDelGenesFile = del.genes.classic,
                                   gisticScoresFile = scores.gis.classic)
  
  
getSampleSummary(dfsp) #Shows gene summary.
getGeneSummary(dfsp) #shows clinical data associated with samples
getClinicalData(dfsp) #Shows all fields in MAF
getFields(dfsp) #Writes maf summary to an output file with basename laml.
write.mafSummary(maf = dfsp, basename = 'dfsp')

plotmafSummary(maf = dfsp, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)


genes <- c("CDKN2A","CDKN2B","TP53", "TERT", "MYC", fst.df)
del.genes.table <- read.delim(del.genes.fst, sep ="\t", header = TRUE)
genes <- c(del.genes$X22q13.31[-(1:3)][1:6],del.genes$X22q13.33[-(1:3)][1:30])
genes <- c(fst.df[2:11], del.genes.table$X22q13.31[-(1:3)][1:6], "CDKN2A","CDKN2B","TP53")
genes <- c(fst.df2, del.genes.table$X22q13.31[-(1:3)][1:6], "CDKN2A","CDKN2B","TP53")


cols = RColorBrewer::brewer.pal(n = 6, name = 'Paired')
hist_cols <- cols[1:2]
names(hist_cols) <- c("Classic", "FST")
hist_cols <- list(Histology.Subtype = hist_cols)
mol_cols <- cols[3:4]
names(mol_cols) <- c("PDGFB", "PDGFD")
mol_cols <- list(Molecular.Subtype = mol_cols)
ori_cols <- cols[5:6]
names(ori_cols) <- c("Primary", "Metastasis")
ori_cols <- list(Origin = ori_cols)

oncoplot(maf = dfsp2, top = 20, 
         clinicalFeatures = c("Histology.Subtype", "Molecular.Subtype", "Origin"), 
         genes = genes,  removeNonMutated = FALSE,
         annotationColor = c(hist_cols, mol_cols, ori_cols),
         sortByAnnotation = TRUE, pathways = NULL, showTumorSampleBarcodes = TRUE,
         anno_height = 1, legend_height = 4, gene_mar = 10, barcode_mar = 6,
         showTitle = FALSE)



mol_cols <- c("#B2DF8A", "#33A02C")
"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C"
"#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00"
cols <- list(Histology.Subtype = hist_cols, Molecular.Subtype = mol_cols)
cnv_cols <- c("#1F78B4", "#E31A1C")
names(cnv_cols) <- c("Del", "Amp")

gisticOncoPlot(gistic = dfsp.gistic, 
               clinicalData = getClinicalData(x = dfsp), 
               clinicalFeatures = c("Histology.Subtype", "Molecular.Subtype", "Origin"),
               sortByAnnotation = TRUE, top = 10, removeNonAltered = FALSE, 
               colors = cnv_cols,
               annotationColor = c(hist_cols, mol_cols, ori_cols),
               gene_mar = 10, barcode_mar = 6, SampleNamefontSize = 1, annotationFontSize = 1.5,
               showTumorSampleBarcodes = TRUE)

gisticChromPlot(gistic = dfsp.gistic.pdgfb, ref.build = "hg38", txtSize = 1.5, cytobandTxtSize = 1.5)
gisticChromPlot(gistic = dfsp.gistic.fst, ref.build = "hg38", markBands = "all",
                txtSize = 1.5, cytobandTxtSize = 1, mutGenesTxtSize = 1.5)

gisticChromPlot(gistic = dfsp.gistic.pdgfd, ref.build = "hg38", txtSize = 1.5, cytobandTxtSize = 1.5)
gisticChromPlot(gistic = dfsp.gistic.classic)

gisticChromPlot(gistic1 = dfsp.pdgfb, gistic2 = dfsp.fst, g1Name = "Classic", g2Name = "FST", type = 'Amp')
coGisticChrom

fst.vs.classic <- mafCompare(m1 = dfsp.fst, m2 = dfsp.classic, 
                             m1Name = 'FST', m2Name = 'Classic', 
                             minMut = 1, useCNV = TRUE,
                             pathways = FALSE)
print(fst.vs.classic)
fst.vs.classic$results[adjPval<0.05]
forestPlot(mafCompareRes = fst.vs.classic, pVal = 0.1)

genes <- print(fst.vs.classic)$results$Hugo_Symbol[1:10]
coOncoplot(m1 = dfsp.classic, m2 = dfsp.fst, m1Name = 'Classic', m2Name = 'FST', genes = genes, removeNonMutated = FALSE)

pdgfb.vs.pdgfd <- mafCompare(m1 = dfsp.pdgfb, m2 = dfsp.pdgfd, 
                             m1Name = 'PDGFB', m2Name = 'PDGFD', 
                             minMut = 2, useCNV = TRUE,
                             pathways = FALSE)
print(pdgfb.vs.pdgfd)
coOncoplot(m1 = dfsp.pdgfb, m2 = dfsp.pdgfd, m1Name = 'PDGFB', m2Name = 'PDGFD', genes = genes, removeNonMutated = FALSE)

pws = pathways(maf = dfsp, plotType = 'treemap')
PlotOncogenicPathways(maf = dfsp, removeNonMutated = FALSE, pathways = c("MYC", "TP53"),  showTumorSampleBarcodes = TRUE)
PlotOncogenicPathways()

fst.ce <- clinicalEnrichment(maf = dfsp, clinicalFeature = "Histology.Subtype", 
                             annotationDat = getClinicalData(x = dfsp),
                             useCNV = TRUE,
                             pathways = FALSE)
fst.ce$groupwise_comparision[p_value <0.05]
fst.df <- fst.ce$pairwise_comparision[fdr<0.05]$Hugo_Symbol
fst.df2 <- fst.ce2$pairwise_comparision[fdr<0.05]$Hugo_Symbol
plotEnrichmentResults(enrich_res = fst.ce, pVal = 0.05, 
                      geneFontSize = 1.5, 
                      annoFontSize = 0.8, 
                      legendFontSize = 1.8)



#Mutational signature
library("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)
BSgenome.Hsapiens.UCSC.hg38
dfsp.tnm = trinucleotideMatrix(maf = dfsp, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")
library('NMF')
dfsp.sign = estimateSignatures(mat = dfsp.tnm, nTry = 6)
plotCophenetic(res = dfsp.sign)
dfsp.sig = extractSignatures(mat = dfsp.tnm, n = 3)
dfsp.v3.cosm = compareSignatures(nmfRes = dfsp.sig, sig_db = "SBS")
library('pheatmap')
pheatmap::pheatmap(mat = dfsp.v3.cosm$cosine_similarities, cluster_rows = FALSE, 
                   main = "cosine similarity against validated signatures")
maftools::plotSignatures(nmfRes = dfsp.sig, title_size = 1.2, sig_db = "SBS")
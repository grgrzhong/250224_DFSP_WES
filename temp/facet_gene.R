
library(GenomicRanges)
library(rtracklayer)
library(VariantAnnotation)
library(here)
library(tidyverse)

# Simple script to convert CNV_facets segments to gene-level CNVs

# Define file paths
vcf_file <- "/home/zhonggr/projects/250224_DFSP_WES/data/mfs/cnv/facets/SARC-022-T/SARC-022-T.vcf.gz"

gtf_file <- "/home/zhonggr/projects/250224_DFSP_WES/data/reference/gene_annotation/gencode.v43.basic.annotation.gtf"

output_file <- "/home/zhonggr/projects/250224_DFSP_WES/data/mfs/cnv/facets/SARC-022-T/SARC-022-T.gene_level_cnvs.txt"


vcf <- readVcf(vcf_file)

# Step 1: Read the VCF file using system command (to avoid complex VCF parsing)
cat("Reading VCF file...\n")
vcf_cmd <- paste("gunzip -c", vcf_file, "| grep -v '^##' | grep -v '^#'")

vcf_data <- read.table(
  text = system(vcf_cmd, intern = TRUE), sep = "\t", stringsAsFactors = FALSE
)
test <- vcf_data |> as_tibble()

colnames(vcf_data) <- c(
  "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"
)

# Extract CNV segments
cat("Extracting CNV segments...\n")
segments <- data.frame(
  chrom = vcf_data$CHROM,
  start = vcf_data$POS,
  stringsAsFactors = FALSE
)

# Extract END positions
segments$end <- sapply(vcf_data$INFO, function(info) {
  end_match <- regexpr("END=[0-9]+", info)
  if (end_match > 0) {
    end_str <- regmatches(info, end_match)
    return(as.numeric(sub("END=", "", end_str)))
  } else {
    return(NA)
  }
})

# Extract total copy number (TCN_EM)
segments$total_cn <- sapply(vcf_data$INFO, function(info) {
  cn_match <- regexpr("TCN_EM=[0-9]+", info)
  if (cn_match > 0) {
    cn_str <- regmatches(info, cn_match)
    return(as.numeric(sub("TCN_EM=", "", cn_str)))
  } else {
    return(NA)
  }
})

# Step 2: Read gene annotations (only gene-level features)
cat("Reading gene annotations...\n")
genes <- import.gff(gtf_file, format="gtf", feature.type="gene")

# Convert to GRanges
cat("Setting up genomic ranges...\n")
segment_gr <- GRanges(
  seqnames = segments$chrom,
  ranges = IRanges(start = segments$start, end = segments$end),
  total_cn = segments$total_cn
)

# Step 3: Find overlaps between genes and segments
cat("Finding overlaps between genes and CNV segments...\n")
overlaps <- findOverlaps(genes, segment_gr)

# Create a data frame with gene-CNV overlaps
cat("Creating gene-level CNV annotations...\n")
gene_cnvs <- data.frame(
  gene_name = genes$gene_name[queryHits(overlaps)],
  gene_id = genes$gene_id[queryHits(overlaps)],
  chromosome = as.character(seqnames(genes)[queryHits(overlaps)]),
  gene_start = start(genes)[queryHits(overlaps)],
  gene_end = end(genes)[queryHits(overlaps)],
  copy_number = segment_gr$total_cn[subjectHits(overlaps)]
)

# Calculate overlap percentages
gene_cnvs$overlap_percent <- sapply(1:nrow(gene_cnvs), function(i) {
  gene_start <- gene_cnvs$gene_start[i]
  gene_end <- gene_cnvs$gene_end[i]
  seg_start <- start(segment_gr)[subjectHits(overlaps)[i]]
  seg_end <- end(segment_gr)[subjectHits(overlaps)[i]]

  overlap_start <- max(gene_start, seg_start)
  overlap_end <- min(gene_end, seg_end)
  overlap_size <- overlap_end - overlap_start + 1
  gene_size <- gene_end - gene_start + 1

  return((overlap_size / gene_size) * 100)
})

# Filter for significant overlaps
significant_cnvs <- gene_cnvs[gene_cnvs$overlap_percent >= 50, ]

# Determine CNV type
significant_cnvs$cnv_type <- "NEUTRAL" # copy number = 2: Normal diploid state in humans

significant_cnvs$cnv_type[significant_cnvs$copy_number > 2] <- "GAIN" # Copy number gain/amplification

significant_cnvs$cnv_type[significant_cnvs$copy_number >= 5] <- "HLAMP" # High-level amplification

significant_cnvs$cnv_type[significant_cnvs$copy_number < 2] <- "LOSS" #  Heterozygous deletion/loss

significant_cnvs$cnv_type[significant_cnvs$copy_number == 0] <- "DLOSS" # Homozygous/deep deletion

# Aggregate by gene (taking the mean of copy numbers and the most common CNV type)
gene_summary <- aggregate(
  significant_cnvs[, c("copy_number", "overlap_percent")],
  by = list(gene_name = significant_cnvs$gene_name),
  FUN = function(x) mean(x, na.rm = TRUE)
)

# Get the most common CNV type for each gene
gene_types <- aggregate(
  data.frame(count = rep(1, nrow(significant_cnvs))),
  by = list(
    gene_name = significant_cnvs$gene_name, 
    cnv_type = significant_cnvs$cnv_type
  ),
  FUN = sum
)

# Find the dominant CNV type
dominant_types <- by(gene_types, gene_types$gene_name, function(x) {
  x[which.max(x$count), "cnv_type"]
})
dominant_df <- data.frame(
  gene_name = names(dominant_types),
  cnv_type = unlist(dominant_types),
  stringsAsFactors = FALSE
)

# Join the results
final_gene_cnvs <- merge(gene_summary, dominant_df, by = "gene_name")

# Write results to file
cat("Writing results to file...\n")
write.table(
  final_gene_cnvs, output_file, sep="\t", quote=FALSE, row.names=FALSE
)

cat("Done! Found", nrow(final_gene_cnvs), "genes with CNV alterations.\n")

final_gene_cnvs[1:5, ]
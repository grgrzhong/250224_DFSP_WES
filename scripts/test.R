# 1. Create a specific FST vs non-FST comparison
# Combine your post_fst and fs_dfsp groups as "FST" samples
fst_samples <- clinical_data %>% 
  filter(sample_group %in% c("post_fst", "fs_dfsp")) %>%
  pull(Tumor_Sample_Barcode)

# All other samples as "non-FST"
non_fst_samples <- clinical_data %>% 
  filter(sample_group %in% c("u_dfsp", "pre_fst")) %>%
  pull(Tumor_Sample_Barcode)

# 2. Create MAF objects for each group
fst_maf <- subsetMaf(maf_filtered, tsb = fst_samples, verbose = FALSE)
non_fst_maf <- subsetMaf(maf_filtered, tsb = non_fst_samples, verbose = FALSE)

# 3. Compare the two groups to find enriched mutations
fst_comparison <- mafCompare(
  m1 = fst_maf, 
  m2 = non_fst_maf,
  m1Name = "FST", 
  m2Name = "non-FST",
  minMut = 2
)

# 4. Extract significant results (adjust p-value threshold as needed)
fst_enriched_genes <- fst_comparison$results %>%
  filter(adjPval < 0.1 & OR > 1) %>%
  arrange(adjPval)

# 5. Save the results
write.csv(fst_enriched_genes, here("results/fst_enriched_mutations.csv"), row.names = FALSE)

# 6. Visualize the results
pdf(here("figures/oncoplot/fst_enriched_forestplot.pdf"), width = 8, height = 6)
forestPlot(mafCompareRes = fst_comparison, pVal = 0.1)
dev.off()

# 7. Create an oncoplot with the enriched genes
if(nrow(fst_enriched_genes) > 0) {
  enriched_genes <- fst_enriched_genes$Hugo_Symbol[1:min(20, nrow(fst_enriched_genes))]
  
  CompareMafPlot(
    maf = maf_filtered,
    genes = enriched_genes,
    clinical_features = "sample_group",
    annotation_colors = sample_group_colors,
    sort_by_annotation = TRUE,
    show_sample_names = FALSE,
    title = "FST-Enriched Mutations",
    font_size = 0.7,
    width = 10,
    height = 8,
    fig_dir = "figures/oncoplot",
    fig_name = "fst_enriched_genes"
  )
}

### Clinical Enrichment Analysis --------------
## Additional Analyses for Finding FST-Enriched Changes
# Add FST status to clinical data
maf_filtered@clinical.data$fst_status <- ifelse(
  maf_filtered@clinical.data$sample_group %in% c("post_fst", "fs_dfsp"), 
  "FST", "non-FST"
)

# Perform clinical enrichment analysis
fst_ce <- clinicalEnrichment(
  maf = maf_filtered, 
  clinicalFeature = "fst_status"
)

# Extract significant results
fst_enriched <- fst_ce$groupwise_comparision[p_value < 0.05]

# Visualize the results
plotEnrichmentResults(
  enrich_res = fst_ce, 
  pVal = 0.05,
  geneFontSize = 1.5,
  annoFontSize = 0.8
)

# Save the results
write.csv(
  fst_ce$groupwise_comparision, 
  here("results/fst_clinical_enrichment.csv"), 
  row.names = FALSE
)

# Pathway analysis for FST samples ----------------
# Check for enriched pathways in FST samples
fst_pathways <- pathways(
  maf = fst_maf,
  plotType = 'barplot'
)

# Plot oncogenic pathways specific to FST
PlotOncogenicPathways(
  maf = fst_maf, 
  removeNonMutated = TRUE
)

# Copy number variations ------------
# If you have GISTIC data for your samples (which appears to be the case from your script), you can analyze copy number differences:
# If you have GISTIC results for FST and non-FST samples separately
# Load or prepare these first, then use coGisticChromPlot
coGisticChromPlot(
  gistic1 = fst_gistic, 
  gistic2 = non_fst_gistic,
  g1Name = "FST", 
  g2Name = "non-FST", 
  type = 'both',  # or 'Amp' or 'Del'
  ref.build = "hg38"
)
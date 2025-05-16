#!/bin/bash
# Variant Filtering Pipeline for FFPE Whole Exome Sequencing Data
# Uses Mutect2 and bcftools with FFPE artifact, oxidative damage, RepeatMasker, and blacklist filtering

set -euo pipefail

# Configuration variables
REFERENCE="reference.fa"
TUMOR_BAM="tumor.bam"
NORMAL_BAM="normal.bam"
NORMAL_SAMPLE="normal_sample"
PON="pon.vcf.gz"
GNOMAD="af-only-gnomad.vcf.gz"
SMALL_EXAC="small_exac_common.vcf"
INTERVALS="intervals.list"
REPEAT_MASKER="repeatmasker.bed"
BLACKLIST="blacklist.bed"
OUTPUT_PREFIX="ffpe_filtered"
THREADS=8
JAVA_OPTS="-Xmx8G -XX:ParallelGCThreads=$THREADS"

# Create header files needed for repeatmasker and blacklist annotation
echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > vcf.rm.header
echo -e "##INFO=<ID=Blacklist,Number=1,Type=String,Description=\"Variant overlaps blacklisted region\">" > vcf.bl.header

# 1. Variant calling with Mutect2
echo "Step 1/8: Calling variants with Mutect2..."
gatk --java-options "$JAVA_OPTS" Mutect2 \
    -R "$REFERENCE" \
    -I "$TUMOR_BAM" \
    -I "$NORMAL_BAM" \
    -normal "$NORMAL_SAMPLE" \
    --germline-resource "$GNOMAD" \
    --panel-of-normals "$PON" \
    --f1r2-tar-gz f1r2.tar.gz \
    -O "${OUTPUT_PREFIX}_unfiltered.vcf.gz"

# 2. Learn orientation bias for FFPE artifacts
echo "Step 2/8: Learning read orientation bias..."
gatk LearnReadOrientationModel \
    -I f1r2.tar.gz \
    -O read-orientation-model.tar.gz

# 3. Calculate contamination
echo "Step 3/8: Calculating contamination..."
gatk GetPileupSummaries \
    -I "$TUMOR_BAM" \
    -V "$SMALL_EXAC" \
    -L "$INTERVALS" \
    -O tumor.pileups.table

gatk GetPileupSummaries \
    -I "$NORMAL_BAM" \
    -V "$SMALL_EXAC" \
    -L "$INTERVALS" \
    -O normal.pileups.table

gatk CalculateContamination \
    -I tumor.pileups.table \
    -matched normal.pileups.table \
    -O contamination.table \
    --tumor-segmentation segments.table

# 4. Filter variants with FilterMutectCalls
echo "Step 4/8: Filtering variants with FilterMutectCalls..."
gatk FilterMutectCalls \
    -R "$REFERENCE" \
    -V "${OUTPUT_PREFIX}_unfiltered.vcf.gz" \
    --ob-priors read-orientation-model.tar.gz \
    --contamination-table contamination.table \
    --tumor-segmentation segments.table \
    -O "${OUTPUT_PREFIX}_mutect_filtered.vcf.gz"

# 5. Annotate with RepeatMasker and blacklist regions
echo "Step 5/8: Annotating with RepeatMasker and blacklist regions..."
# First annotate with RepeatMasker
bcftools annotate \
    --header-lines vcf.rm.header \
    --annotations "$REPEAT_MASKER" \
    --columns CHROM,FROM,TO,RepeatMasker \
    "${OUTPUT_PREFIX}_mutect_filtered.vcf.gz" \
    -Oz -o "${OUTPUT_PREFIX}_repmask_annotated.vcf.gz"

# Then annotate with blacklist
bcftools annotate \
    --header-lines vcf.bl.header \
    --annotations "$BLACKLIST" \
    --columns CHROM,FROM,TO,Blacklist \
    "${OUTPUT_PREFIX}_repmask_annotated.vcf.gz" \
    -Oz -o "${OUTPUT_PREFIX}_annotated.vcf.gz"
bcftools index "${OUTPUT_PREFIX}_annotated.vcf.gz"

# 6. Apply quality filters and region filters
echo "Step 6/8: Applying quality and region filters..."
# Apply quality filters
bcftools filter -i 'QUAL>=30 && INFO/DP>=20 && INFO/AF>=0.05' \
    "${OUTPUT_PREFIX}_annotated.vcf.gz" \
    -Oz -o "${OUTPUT_PREFIX}_quality_filtered.vcf.gz"

# Filter out variants in RepeatMasker regions
bcftools view -i 'INFO/RepeatMasker=="."' \
    "${OUTPUT_PREFIX}_quality_filtered.vcf.gz" \
    -Oz -o "${OUTPUT_PREFIX}_no_repeats.vcf.gz"

# Filter out variants in blacklist regions
bcftools view -i 'INFO/Blacklist=="."' \
    "${OUTPUT_PREFIX}_no_repeats.vcf.gz" \
    -Oz -o "${OUTPUT_PREFIX}_no_blacklist.vcf.gz"

# Filter out common variants (gnomAD population frequency > 0.01)
bcftools filter -i 'INFO/gnomAD_AF<=0.01 || INFO/gnomAD_AF="."' \
    "${OUTPUT_PREFIX}_no_blacklist.vcf.gz" \
    -Oz -o "${OUTPUT_PREFIX}_rare_variants.vcf.gz"

# 7. FFPE artifact and oxidative damage-specific filtering
echo "Step 7/8: Applying FFPE artifact and oxidative damage filters..."
bcftools filter -i '
    # Keep non-FFPE and non-oxidative damage variants
    ((REF!="C" | (ALT!="T" & ALT!="A")) & (REF!="G" | (ALT!="A" & ALT!="T"))) |
    # Keep FFPE C>T variants with DP_5prime_end<5 AND AF>=0.10
    (REF="C" & ALT="T" & INFO/DP_5prime_end<5 & FORMAT/AF>=0.10) |
    # Keep FFPE G>A variants with DP_5prime_end<5 AND AF>=0.10
    (REF="G" & ALT="A" & INFO/DP_5prime_end<5 & FORMAT/AF>=0.10) |
    # Keep oxidative damage C>A variants with AF>=0.10
    (REF="C" & ALT="A" & FORMAT/AF>=0.10) |
    # Keep oxidative damage G>T variants with AF>=0.10
    (REF="G" & ALT="T" & FORMAT/AF>=0.10)
' "${OUTPUT_PREFIX}_rare_variants.vcf.gz" \
    -Oz -o "${OUTPUT_PREFIX}_ffpe_filtered.vcf.gz"

# 8. Final output and statistics
echo "Step 8/8: Generating final output and statistics..."
# Rename the final filtered file
mv "${OUTPUT_PREFIX}_ffpe_filtered.vcf.gz" "${OUTPUT_PREFIX}_final.vcf.gz"
bcftools index "${OUTPUT_PREFIX}_final.vcf.gz"

# Generate statistics
bcftools stats "${OUTPUT_PREFIX}_final.vcf.gz" > "${OUTPUT_PREFIX}_final_stats.txt"
bcftools +tstv "${OUTPUT_PREFIX}_final.vcf.gz" > "${OUTPUT_PREFIX}_tstv_ratio.txt"

# Create a filtered VCF that keeps all variants but marks filtered ones in the FILTER field
bcftools view "${OUTPUT_PREFIX}_annotated.vcf.gz" | \
    awk -F'\t' -v OFS='\t' 'BEGIN {
                   print "##FILTER=<ID=LowQual,Description=\"Low quality: QUAL<30 or DP<20 or AF<0.05\">";
                   print "##FILTER=<ID=CommonVariant,Description=\"Variant with gnomAD AF>0.01\">";
                   print "##FILTER=<ID=RepeatRegion,Description=\"Variant in RepeatMasker region\">";
                   print "##FILTER=<ID=BlacklistRegion,Description=\"Variant in blacklisted region\">";
                   print "##FILTER=<ID=FFPEArtifact,Description=\"Potential FFPE artifact (C>T or G>A at read end with AF<0.10)\">";
                   print "##FILTER=<ID=OxidativeDamage,Description=\"Potential oxidative damage variant (C>A or G>T with AF<0.10)\">";
                  }
        /^#/ {print; next}
        {
            filter = "PASS";
            if ($6 < 30 || !match($8, /DP=[0-9]+/) || (match($8, /DP=([0-9]+)/, m) && m[1] < 20)) filter = "LowQual";
            if (match($8, /AF=([0-9.]+)/, m) && m[1] < 0.05) filter = (filter == "PASS" ? "LowQual" : filter ";LowQual");
            if (match($8, /gnomAD_AF=([0-9.]+)/, m) && m[1] > 0.01) filter = (filter == "PASS" ? "CommonVariant" : filter ";CommonVariant");
            if (match($8, /RepeatMasker!="."/) || match($8, /RepeatMasker=1/)) filter = (filter == "PASS" ? "RepeatRegion" : filter ";RepeatRegion");
            if (match($8, /Blacklist!="."/) || match($8, /Blacklist=1/)) filter = (filter == "PASS" ? "BlacklistRegion" : filter ";BlacklistRegion");
            
            # FFPE artifacts (C>T or G>A) with AF<0.10 or at read end
            if (($4 == "C" && $5 == "T") || ($4 == "G" && $5 == "A")) {
                if ((match($8, /DP_5prime_end=([0-9]+)/, m) && m[1] >= 5) || (match($9, /AF=([0-9.]+)/, m) && m[1] < 0.10)) {
                    filter = (filter == "PASS" ? "FFPEArtifact" : filter ";FFPEArtifact");
                }
            }
            
            # Oxidative damage artifacts (C>A or G>T) with AF<0.10
            if (($4 == "C" && $5 == "A") || ($4 == "G" && $5 == "T")) {
                if (match($9, /AF=([0-9.]+)/, m) && m[1] < 0.10) {
                    filter = (filter == "PASS" ? "OxidativeDamage" : filter ";OxidativeDamage");
                }
            }
            
            $7 = filter;
            print
        }' | bgzip > "${OUTPUT_PREFIX}_fully_annotated.vcf.gz"

bcftools index "${OUTPUT_PREFIX}_fully_annotated.vcf.gz"

# Clean up temporary header files
rm vcf.rm.header vcf.bl.header

echo "Pipeline completed successfully!"
echo "Final filtered variants: ${OUTPUT_PREFIX}_final.vcf.gz"
echo "Fully annotated variants with filter reasons: ${OUTPUT_PREFIX}_fully_annotated.vcf.gz"
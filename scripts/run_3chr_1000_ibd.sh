#!/bin/bash
# ===============================================
# Phase 1 & 2 — 1000 Genomes subset, chr20, chr21, chr22
# 1000 samples, 3 small chromosomes, Beagle phasing + Refined IBD
# All paths relative to the cse284-final-project directory
# ===============================================
set -euo pipefail

PROJDIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJDIR"

DATARAW="data_raw"
DATAQC="data_qc"
RESULTS="results"
BEAGLE_JAR="beagle.27Feb25.75f.jar"
REFINEDIBD_JAR="refined-ibd.17Jan20.102.jar"

mkdir -p "$DATARAW" "$DATAQC" "$RESULTS"

# Step 1: Download VCFs, indices, and sample panel (if not already present)
for chr in 20 21 22; do
    vcf="$DATARAW/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    [ -f "$vcf" ] || wget -P "$DATARAW" \
      "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    [ -f "${vcf}.tbi" ] || wget -P "$DATARAW" \
      "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"
done

PANEL="$DATARAW/integrated_call_samples_v3.20130502.ALL.panel"
[ -f "$PANEL" ] || wget -P "$DATARAW" \
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

# Step 2: Select 1000 random sample IDs (skip header, extract sample column)
SAMPLES="$DATARAW/final_samples_1000.txt"
# tail -n +2 "$PANEL" | cut -f1 | shuf | head -n 1000 > "$SAMPLES"
# echo "Selected $(wc -l < "$SAMPLES") sample IDs"

# Step 3: Subset VCFs to selected samples & normalize to biallelic
for chr in 20 21 22; do
    vcf="$DATARAW/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    subset="$DATAQC/chr${chr}_subset.vcf.gz"
    biallelic="$DATAQC/chr${chr}_subset_biallelic.vcf.gz"

    echo "--- chr${chr}: subsetting samples ---"
    bcftools view -S "$SAMPLES" --force-samples "$vcf" -Oz -o "$subset"
    bcftools index -f "$subset"

    echo "--- chr${chr}: normalizing to biallelic ---"
    bcftools norm -m -both -Oz -o "$biallelic" "$subset"
    bcftools index -f "$biallelic"
done

# Step 4: Remove sites with missing genotypes, then deduplicate markers
# Refined IBD / Beagle require clean, unique markers.
for chr in 20 21 22; do
    biallelic="$DATAQC/chr${chr}_subset_biallelic.vcf.gz"
    nomiss="$DATAQC/chr${chr}_subset_biallelic_nomiss.vcf.gz"
    dedup="$DATAQC/chr${chr}_subset_biallelic_nomiss_dedup.vcf.gz"

    echo "--- chr${chr}: removing sites with missing genotypes ---"
    bcftools view -g ^miss "$biallelic" -Oz -o "$nomiss"
    bcftools index -f "$nomiss"

    echo "--- chr${chr}: removing duplicate markers ---"
    bcftools norm -d all -Oz -o "$dedup" "$nomiss"
    bcftools index -f "$dedup"
done

# Step 5: Phase with Beagle
for chr in 20 21 22; do
    dedup="$DATAQC/chr${chr}_subset_biallelic_nomiss_dedup.vcf.gz"
    outprefix="$RESULTS/chr${chr}_phased"

    echo "--- chr${chr}: phasing with Beagle ---"
    java -Xmx8g -jar "$BEAGLE_JAR" \
      gt="$dedup" \
      out="$outprefix" \
      nthreads=8
done

# Step 6: Run Refined IBD
# Using permissive thresholds to capture shorter segments on small chromosomes.
for chr in 20 21 22; do
    phased="$RESULTS/chr${chr}_phased.vcf.gz"
    outprefix="$RESULTS/chr${chr}_refinedibd"

    echo "--- chr${chr}: running Refined IBD ---"
    java -Xmx16g -jar "$REFINEDIBD_JAR" \
      gt="$phased" \
      out="$outprefix" \
      nthreads=8 \
      lod=1.0 \
      length=0.5
done

# Step 7: Summarize number of detected segments per chromosome
SUMMARY="$RESULTS/refinedibd_segment_counts.txt"
: > "$SUMMARY"

for chr in 20 21 22; do
    ibdfile="$RESULTS/chr${chr}_refinedibd.ibd.gz"
    count=$(gunzip -c "$ibdfile" | wc -l | awk '{print $1}')
    echo "chr${chr}  ${count}" >> "$SUMMARY"
done

echo "Segment counts written to $SUMMARY"

# Step 8: Concatenate all segment outputs into one combined file
COMBINED="$RESULTS/refinedibd_3chr_all.ibd"
: > "$COMBINED"

for chr in 20 21 22; do
    gunzip -c "$RESULTS/chr${chr}_refinedibd.ibd.gz" >> "$COMBINED"
done

gzip -f "$COMBINED"
echo "Combined Refined IBD segments written to ${COMBINED}.gz"

# Step 9: Generate a simple pairwise summary (pair, number of segments, total shared cM)
PAIRWISE="$RESULTS/refinedibd_3chr_pairwise_summary.tsv"

gunzip -c "${COMBINED}.gz" | \
awk '{
    a=$1; b=$3;
    if (a>b) {t=a; a=b; b=t}
    key=a"\t"b;
    n[key]+=1;
    sum[key]+=$9;
    if ($9>mx[key]) mx[key]=$9;
}
END {
    print "ID1\tID2\tn_segments\ttotal_shared_cm\tmax_segment_cm";
    for (k in n) {
        print k "\t" n[k] "\t" sum[k] "\t" mx[k];
    }
}' > "$PAIRWISE"

echo "Pairwise summary written to $PAIRWISE"

echo "All steps complete. Refined IBD segment files and pairwise summaries are ready."
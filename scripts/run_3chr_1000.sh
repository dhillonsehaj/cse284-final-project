#!/bin/bash
# ===============================================
# Phase 1 & 2 — 1000 Genomes subset, chr20, chr21, chr22
# 1000 samples, 3 small chromosomes, PLINK2 KING analysis
# All paths relative to the cse284-final-project directory
# ===============================================
set -euo pipefail

PROJDIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJDIR"

DATARAW="data_raw"
DATAQC="data_qc"
RESULTS="results"

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
tail -n +2 "$PANEL" | cut -f1 | shuf | head -n 1000 > "$SAMPLES"
echo "Selected $(wc -l < "$SAMPLES") sample IDs"

# Step 3: Subset VCFs to selected samples & normalize to biallelic
for chr in 20 21 22; do
    vcf="$DATARAW/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    subset="$DATAQC/chr${chr}_subset.vcf.gz"
    biallelic="$DATAQC/chr${chr}_subset_biallelic.vcf.gz"

    echo "--- chr${chr}: subsetting samples ---"
    bcftools view -S "$SAMPLES" --force-samples "$vcf" -Oz -o "$subset"
    bcftools index "$subset"

    echo "--- chr${chr}: normalizing to biallelic ---"
    bcftools norm -m -both -Oz -o "$biallelic" "$subset"
    bcftools index "$biallelic"
done

# Step 4: Convert to PLINK2 binary format
for chr in 20 21 22; do
    plink2 \
      --vcf "$DATAQC/chr${chr}_subset_biallelic.vcf.gz" \
      --make-bed \
      --out "$DATAQC/chr${chr}_subset"
done

# Step 5: Apply QC filters
for chr in 20 21 22; do
    plink2 \
      --bfile "$DATAQC/chr${chr}_subset" \
      --geno 0.02 --mind 0.02 --maf 0.05 \
      --make-bed \
      --out "$DATAQC/chr${chr}_qc"
done

# Step 6: Assign unique variant IDs per chromosome, then merge
for chr in 20 21 22; do
    plink2 \
      --bfile "$DATAQC/chr${chr}_qc" \
      --set-missing-var-ids '@:#:\$r:\$a' \
      --new-id-max-allele-len 500 \
      --make-bed \
      --out "$DATAQC/chr${chr}_qc_ids"
done

MERGELIST="$DATAQC/merge_list.txt"
echo -e "$DATAQC/chr21_qc_ids\n$DATAQC/chr22_qc_ids" > "$MERGELIST"

plink2 \
  --bfile "$DATAQC/chr20_qc_ids" \
  --pmerge-list "$MERGELIST" bfile \
  --make-bed \
  --out "$DATAQC/merged_3chr_1000"

# Step 7: Run KING to detect relatives (PLINK2)
plink2 \
  --bfile "$DATAQC/merged_3chr_1000" \
  --make-king-table \
  --out "$RESULTS/king_3chr_1000"

# Step 8: Run method-of-moments IBD estimation (PLINK 1.9)
plink \
  --bfile "$DATAQC/merged_3chr_1000" \
  --genome \
  --out "$RESULTS/genome_3chr_1000"

# Step 9: List putative relatives (3rd-degree or closer, kinship > 0.0442)
awk 'NR>1 && $NF>0.0442 {print $1, $2, $3, $4, $NF}' "$RESULTS/king_3chr_1000.kin0" > "$RESULTS/relatives_detected.txt"
echo "Pairs of putative relatives written to $RESULTS/relatives_detected.txt"

echo "All steps complete. KING table, IBD estimates, and relatives list are ready."

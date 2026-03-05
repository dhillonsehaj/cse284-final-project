# CSE 284 Final Project — Benchmarking IBD Detection Methods for Relative Finding

## Project Overview

This project benchmarks three identity-by-descent (IBD) detection methods for relative finding: **PLINK** (genotype-based), **GERMLINE** (haplotype-based), and **Refined IBD / Beagle** (probabilistic haplotype-based). Using a subset of the 1000 Genomes Project, we evaluate each method on multiple quantitative and computational metrics to assess their sensitivity, consistency, and efficiency.

Specifically, we compare:
- Number of detected IBD pairs and segments
- Total shared segment length per pair and segment length distributions
- Agreement between methods (pairwise overlap and Jaccard similarity)
- Genomic overlap of detected segments (for segment-based approaches)
- Computational performance (runtime and memory usage)
- Inferred relatedness via kinship estimates and total IBD sharing

## Dataset

### Source: 1000 Genomes Project

The data comes from the [1000 Genomes Project](https://www.internationalgenome.org/), a large international initiative characterizing human genetic variation across diverse populations. The dataset includes whole-genome sequencing data from over 2,500 individuals representing 26 populations worldwide, spanning African, European, East Asian, South Asian, and Admixed American ancestries. All samples were processed using standardized pipelines and rigorous quality control procedures.

### Subset Used

To limit computational load, we use **chromosome 22** from the Phase 3 v5b release:
- **File**: `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
- **Intended subset**: 350 individuals selected from the sample manifest
- **Actual subset**: 230 individuals (after filtering against samples present in the VCF)

### Data Processing and Quality Control

1. **Sample selection**: 350 samples of interest were selected from the 1000 Genomes sample info file (`20130606_sample_info.txt`). After cross-referencing with the VCF header, 230 samples were retained.

2. **VCF subsetting**: The chr22 VCF was subset to the 230 samples using `bcftools view --samples-file`.

3. **Normalization**: Multiallelic variants were split into biallelic records using `bcftools norm -m -both` to avoid downstream issues.

4. **PLINK conversion**: The normalized VCF was converted to PLINK binary format (`.bed/.bim/.fam`) using PLINK2.

5. **Quality control filters** (applied via PLINK2):
   - Remove SNPs with >2% missing genotypes (`--geno 0.02`)
   - Remove individuals with >2% missing genotypes (`--mind 0.02`)
   - Remove variants with minor allele frequency <5% (`--maf 0.05`)

6. **LD pruning**: Skipped due to duplicate variant IDs in the dataset and small sample size.

**Final QC'd dataset**: 230 samples × 118,164 biallelic SNPs on chromosome 22.

### File Structure

```
data_raw/           Raw and intermediate data files
  20130606_sample_info.txt          1000 Genomes sample manifest
  selected_samples.txt              Initial 350-sample selection
  final_samples_filtered.txt        230 samples present in VCF
  vcf_samples.txt                   Sample IDs from VCF header

data_qc/            Quality-controlled PLINK files
  chr22_qc.bed / .bim / .fam       Final QC'd binary genotype data
  chr22_subset.bed / .bim / .fam    Pre-QC binary genotype data

results/            Analysis outputs
  plink_king.kin0                   PLINK2 KING kinship estimates
  chr22_genome.genome               PLINK 1.9 IBD proportions (Z0/Z1/Z2/PI_HAT)
  plink_results_summary.txt         Summary of PLINK metrics and results
  results_guide.txt                 Column definitions and format guide

scripts/            Analysis scripts
```

---

## Method 1: PLINK

### Approach

PLINK estimates pairwise relatedness directly from unphased genotype data, without requiring haplotype phasing. We ran two complementary analyses:

1. **PLINK2 KING-robust kinship** (`--make-king-table`): Computes the KING-robust kinship estimator for all sample pairs. This method is robust to population structure and provides a single kinship coefficient per pair.

2. **PLINK 1.9 method-of-moments IBD** (`--genome`): Estimates the proportions of the genome shared IBD=0, IBD=1, and IBD=2 (Z0, Z1, Z2) for all sample pairs using a method-of-moments approach, along with the overall IBD proportion (PI_HAT = Z1/2 + Z2).

### Commands

```bash
# KING-robust kinship (PLINK2 v2.0.0-a.7 M1)
plink2 --bfile data_qc/chr22_qc --make-king-table --out results/plink_king

# Method-of-moments IBD proportions (PLINK 1.9.0-b.7.11)
plink --bfile data_qc/chr22_qc --genome --out results/chr22_genome
```

### Results

**26,335 pairwise comparisons** were computed (all pairs among 230 samples).

KING kinship summary:
| Statistic | Value |
|-----------|-------|
| Mean kinship | -0.1550 |
| Max kinship | 0.0932 |
| Min kinship | -0.4720 |
| Pairs with kinship > 0.0442 (3rd degree+) | 255 |
| Pairs with kinship > 0.0884 (2nd degree+) | 0 |
| Pairs with kinship > 0.177 (1st degree+) | 0 |

IBD proportions summary:
| Statistic | Value |
|-----------|-------|
| Mean PI_HAT | 0.0384 |
| Max PI_HAT | 0.3012 |
| Pairs with PI_HAT > 0.125 | 4,530 |
| Pairs with PI_HAT > 0.25 | 277 |

No close relatives (1st or 2nd degree) were detected by KING kinship. The 255 pairs exceeding the 3rd-degree threshold are consistent with population substructure among 1000 Genomes samples rather than true recent relatedness.

### Metrics for Benchmarking

PLINK provides the following **pair-level metrics** that are directly comparable to aggregated outputs from GERMLINE and Refined IBD:

| Metric | Source | Column |
|--------|--------|--------|
| Kinship coefficient | `plink_king.kin0` | KINSHIP |
| P(IBD=0) | `chr22_genome.genome` | Z0 |
| P(IBD=1) | `chr22_genome.genome` | Z1 |
| P(IBD=2) | `chr22_genome.genome` | Z2 |
| Total IBD sharing | `chr22_genome.genome` | PI_HAT |
| Number of IBD pairs | Derived from thresholds on KINSHIP or PI_HAT | — |
| IBS distance | `chr22_genome.genome` | DST |

**Note**: PLINK does not produce segment-level output (individual IBD segments with genomic coordinates). Segment-level comparisons (segment counts, length distributions, genomic overlap, Jaccard similarity) will only be evaluated between GERMLINE and Refined IBD.

### Computational Performance

| Analysis | Runtime | Memory Reserved | Threads |
|----------|---------|-----------------|---------|
| PLINK2 KING | < 1 second | 8,192 MiB | 8 |
| PLINK 1.9 --genome | < 1 second | 8,192 MB | 8 |

---

## Method 2: GERMLINE

*To be completed.*

---

## Method 3: Refined IBD (Beagle)

*To be completed.*
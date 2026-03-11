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

To limit computational load, we use **chromosomes 20, 21, and 22** from the Phase 3 v5b release:
- **Files**: `ALL.chr{20,21,22}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
- **Sample panel**: `integrated_call_samples_v3.20130502.ALL.panel`
- **Intended subset**: 1,000 individuals randomly selected from the sample panel
- **Actual subset**: 1,000 individuals (all present in the VCFs)

### Data Processing and Quality Control

The full pipeline is implemented in `scripts/run_3chr_1000.sh` and can be reproduced by running:

```bash
conda activate ibd_env
bash scripts/run_3chr_1000.sh
```

Steps performed by the script:

1. **Download**: VCFs, tabix indices, and sample panel are downloaded from the 1000 Genomes FTP (skipped if already present).

2. **Sample selection**: 1,000 sample IDs are randomly selected from the panel file.

3. **VCF subsetting**: Each chromosome VCF is subset to the selected samples using `bcftools view -S`.

4. **Normalization**: Multiallelic variants are split into biallelic records using `bcftools norm -m -both`.

5. **PLINK conversion**: Normalized VCFs are converted to PLINK binary format (`.bed/.bim/.fam`) using PLINK2.

6. **Quality control filters** (applied via PLINK2):
   - Remove SNPs with >2% missing genotypes (`--geno 0.02`)
   - Remove individuals with >2% missing genotypes (`--mind 0.02`)
   - Remove variants with minor allele frequency <5% (`--maf 0.05`)

7. **Variant ID assignment**: Missing variant IDs (`.`) are replaced with unique identifiers using `--set-missing-var-ids` to avoid merge conflicts.

8. **Chromosome merging**: The three per-chromosome QC'd files are merged into a single dataset using `--pmerge-list`.

**Final QC'd dataset**: 1,000 samples × 401,845 biallelic SNPs across chromosomes 20, 21, and 22.

### File Structure

```
data_raw/           Raw and intermediate data files
  integrated_call_samples_v3.20130502.ALL.panel   1000 Genomes sample panel
  final_samples_1000.txt            1,000 randomly selected sample IDs
  ALL.chr{20,21,22}...vcf.gz        Downloaded VCFs (gitignored)

data_qc/            Quality-controlled PLINK files
  chr{20,21,22}_qc.bed/.bim/.fam    Per-chromosome QC'd binary genotype data
  chr{20,21,22}_qc_ids.bed/.bim/.fam Per-chromosome with unique variant IDs
  merged_3chr_1000.bed/.bim/.fam    Merged 3-chromosome QC'd dataset

results/            Analysis outputs
  king_3chr_1000.kin0               PLINK2 KING kinship estimates (1000 samples)
  genome_3chr_1000.genome           PLINK 1.9 --genome IBD estimates (Z0/Z1/Z2/PI_HAT)
  relatives_detected.txt            Pairs with kinship > 0.0442

scripts/            Analysis and pipeline scripts
  run_3chr_1000.sh                  End-to-end pipeline: download, QC, merge, KING
  analyze_germline2.py              GERMLINE2 post-processing and analysis
  analyze_refinedibd.py             Refined IBD post-processing and analysis
```

---

## Method 1: PLINK

### Approach

PLINK estimates pairwise relatedness directly from unphased genotype data, without requiring haplotype phasing. We use **PLINK2 KING-robust kinship** (`--make-king-table`), which computes the KING-robust kinship estimator for all sample pairs. This method is robust to population structure and provides a single kinship coefficient per pair.

### Commands

The full pipeline (download, QC, merge, KING, --genome) is automated in `scripts/run_3chr_1000.sh`. The key commands are:

```bash
# Per-chromosome QC
plink2 --vcf data_qc/chr${chr}_subset_biallelic.vcf.gz --make-bed --out data_qc/chr${chr}_subset
plink2 --bfile data_qc/chr${chr}_subset --geno 0.02 --mind 0.02 --maf 0.05 --make-bed --out data_qc/chr${chr}_qc

# Assign unique variant IDs and merge chromosomes
plink2 --bfile data_qc/chr${chr}_qc --set-missing-var-ids '@:#:$r:$a' --new-id-max-allele-len 500 --make-bed --out data_qc/chr${chr}_qc_ids
plink2 --bfile data_qc/chr20_qc_ids --pmerge-list data_qc/merge_list.txt bfile --make-bed --out data_qc/merged_3chr_1000

# KING-robust kinship (PLINK2)
plink2 --bfile data_qc/merged_3chr_1000 --make-king-table --out results/king_3chr_1000

# Method-of-moments IBD estimation (PLINK 1.9)
plink --bfile data_qc/merged_3chr_1000 --genome --out results/genome_3chr_1000
```

### Results

**499,500 pairwise comparisons** were computed (all pairs among 1,000 samples across chromosomes 20, 21, and 22).

- Samples selected: 1,000
- Samples kept after QC: 1,000 (0 removed by `--mind 0.02`)
- Variants after QC: chr20: 172,665 / chr21: 116,095 / chr22: 113,085
- Total merged variants: 401,845

#### KING-robust Kinship (PLINK2)

| Relationship category | Kinship threshold | Pairs detected |
|-----------------------|-------------------|----------------|
| 1st degree | > 0.177 | 1 |
| 2nd degree | 0.0884 – 0.177 | 1 |
| 3rd degree | 0.0442 – 0.0884 | 321 |
| Unrelated | < 0.0442 | 499,177 |
| **Total related pairs** | **> 0.0442** | **323** |

Top kinship value: **0.2446** (pair NA20320–NA20321, 1st-degree relatives). All 323 related pairs are written to `results/relatives_detected.txt`.

#### Method-of-Moments IBD (PLINK 1.9 `--genome`)

The `--genome` command estimates IBD state probabilities (Z0, Z1, Z2) and the proportion of genome shared IBD (PI_HAT) for every pair using a method-of-moments approach. These metrics provide a finer-grained view of relatedness compared to KING's single kinship coefficient.

| Relationship category | PI_HAT threshold | Pairs detected |
|-----------------------|------------------|----------------|
| 1st degree | ≥ 0.5 | 1 |
| 2nd degree | 0.25 – 0.5 | 36 |
| 3rd degree | 0.125 – 0.25 | 29,428 |
| Unrelated | < 0.0625 | 414,300 |

**Top pair**: NA20320–NA20321 with PI_HAT = 0.5049, Z0 = 0.0172, Z1 = 0.9558, Z2 = 0.0270. The very high Z1 and near-zero Z0/Z2 are consistent with a **parent-child relationship** (expected: Z0=0, Z1=1, Z2=0).

Both methods agree on the top pair (NA20320–NA20321). However, the --genome method detects more pairs at lower relatedness thresholds, while KING's robust estimator is more conservative and better handles population structure.

### Metrics for Benchmarking

PLINK provides the following **pair-level metrics** that are directly comparable to aggregated outputs from GERMLINE and Refined IBD:

| Metric | Source | Column |
|--------|--------|--------|
| Kinship coefficient | `king_3chr_1000.kin0` | KINSHIP |
| P(IBD=0) | `genome_3chr_1000.genome` | Z0 |
| P(IBD=1) | `genome_3chr_1000.genome` | Z1 |
| P(IBD=2) | `genome_3chr_1000.genome` | Z2 |
| Proportion IBD | `genome_3chr_1000.genome` | PI_HAT |
| Number of IBD pairs | Derived from thresholds on KINSHIP or PI_HAT | — |

**Note**: PLINK does not produce segment-level output (individual IBD segments with genomic coordinates). Segment-level comparisons (segment counts, length distributions, genomic overlap, Jaccard similarity) will only be evaluated between GERMLINE and Refined IBD.

### Computational Performance

| Analysis | Runtime | Threads |
|----------|---------|---------|
| Full pipeline (download + QC + merge + KING + --genome) | ~7 min | 8 |
| KING kinship only | ~2 seconds | 8 |
| `--genome` only | ~1.5 seconds | 8 |

---

## Method 2: GERMLINE2

### Approach

GERMLINE2 is a haplotype-based method for detecting **identity-by-descent (IBD) segments**, defined as genomic regions shared between individuals due to inheritance from a recent common ancestor. The software identifies long shared haplotypes by first locating short exact matches between haplotypes (seeds) and then extending them across neighboring markers while allowing a limited number of mismatches. Segments exceeding a user-specified minimum genetic length are reported as candidate IBD segments.

GERMLINE2 is designed to scale efficiently to large genomic datasets by using a hash-based search strategy that identifies candidate matches before extending them into full segments. The algorithm outputs genomic coordinates for each detected segment, enabling detailed analysis of IBD segment lengths, counts, and genomic distribution.

Unlike PLINK’s allele-frequency–based approach, GERMLINE2 operates on **phased haplotype data** and therefore provides explicit **segment-level information** rather than only pairwise relatedness coefficients.

For this project we used the open-source implementation available from the Gusev Lab:

- **GERMLINE2**
- Repository: https://github.com/gusevlab/germline2

The source code was downloaded and compiled locally using:
git clone https://github.com/gusevlab/germline2
cd germline2
make

This produces the executable `g2`, which performs IBD detection on phased haplotype input files.

### Data Preparation

GERMLINE2 requires phased haplotypes in SHAPEIT/IMPUTE hap/sample format, along with a sample file containing sample identifiers and a genetic map specifying marker positions and genetic distances.

The starting dataset consisted of phased 1000 Genomes Project Phase 3 VCF files for chromosomes 20, 21, and 22. These files were already phased, so no additional phasing with Beagle was required. We analyzed a subset of 1000 individuals listed in:

data_raw/final_samples_1000.txt

#### Subsetting the chromosome VCFs

The phased chromosome VCFs were subset to the 1000 target individuals using bcftools view:
bcftools view -S data_raw/final_samples_1000.txt -Oz \
  -o data_raw/chr20_1000.vcf.gz \
  data_raw/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

bcftools view -S data_raw/final_samples_1000.txt -Oz \
  -o data_raw/chr21_1000.vcf.gz \
  data_raw/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

bcftools view -S data_raw/final_samples_1000.txt -Oz \
  -o data_raw/chr22_1000.vcf.gz \
  data_raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

These subset VCFs were then indexed with bcftools index.

#### Converting VCF to haplotype format

The subset VCFs were converted to GERMLINE2 input format using bcftools convert:

bcftools convert --hapsample data_raw/chr20_1000_g2 data_raw/chr20_1000.vcf.gz
bcftools convert --hapsample data_raw/chr21_1000_g2 data_raw/chr21_1000.vcf.gz
bcftools convert --hapsample data_raw/chr22_1000_g2 data_raw/chr22_1000.vcf.gz

This produced:

- data_raw/chr20_1000_g2.hap.gz
- data_raw/chr20_1000_g2.samples
- data_raw/chr21_1000_g2.hap.gz
- data_raw/chr21_1000_g2.samples
- data_raw/chr22_1000_g2.hap.gz
- data_raw/chr22_1000_g2.samples

Because GERMLINE2 was run on uncompressed haplotype files, the .hap.gz files were decompressed to .hap before running g2.

#### Constructing the genetic map

GERMLINE2 requires a genetic map with three columns:

physical_position cM/Mb cumulative_cM

For this analysis, map files were stored in data_raw/:
- data_raw/chr20_g2.map
- data_raw/chr21_g2.map
- data_raw/chr22_g2.map

The cumulative centimorgan column was required to increase with genomic position. The chromosome 22 map, for example, spanned approximately 60.35 cM, which is a realistic scale for that chromosome. These map files were used directly in the final GERMLINE2 runs.

### Commands

GERMLINE2 was run in two configurations across chromosomes 20, 21, and 22.

A key parameter change during debugging was setting -d 0 instead of -d 1. Using -d 1 removed too many seeds and led to empty output files. Setting -d 0 allowed segments to be detected successfully.

#### Diploid segment detection

This run identifies IBD segments between individuals without distinguishing haplotypes. Final parameters were:
- -m 1.0
- -g 2
- -d 0
- -f 0.05

Example loop:
for chr in 20 21 22; do
  /usr/bin/time -v ./g2 \
    -m 1.0 \
    -g 2 \
    -d 0 \
    -f 0.05 \
    /mnt/c/Users/Luvly/Documents/cse284-final-project/data_raw/chr${chr}_1000_g2.hap \
    /mnt/c/Users/Luvly/Documents/cse284-final-project/data_raw/chr${chr}_1000_g2.samples \
    /mnt/c/Users/Luvly/Documents/cse284-final-project/data_raw/chr${chr}_g2.map \
    /mnt/c/Users/Luvly/Documents/cse284-final-project/results/GERMLINE2/chr${chr}_1000_g2_out \
    2> /mnt/c/Users/Luvly/Documents/cse284-final-project/results/GERMLINE2/chr${chr}_1000_g2_runtime.log
done

Parameter descriptions:
| Parameter | Description |
|-----------|-------------|
| -m | Minimum segment length in cM |
| -g | Maximum allowed gaps |
| -d | Dynamic seed threshold |
| -f | Minor allele frequency filter |

#### Haploid mode (IBD state estimation)

To estimate IBD state probabilities, GERMLINE2 was also run in haploid mode, which restricts matches to a single haplotype and appends .0 or .1 to sample identifiers.

Final haploid parameters matched the diploid minimum length:
- -h
- -m 1.0
- -g 2
- -d 0
- -f 0.05

Example loop:
for chr in 20 21 22; do
  /usr/bin/time -v ./g2 \
    -h \
    -m 1.0 \
    -g 2 \
    -d 0 \
    -f 0.05 \
    /mnt/c/Users/Luvly/Documents/cse284-final-project/data_raw/chr${chr}_1000_g2.hap \
    /mnt/c/Users/Luvly/Documents/cse284-final-project/data_raw/chr${chr}_1000_g2.samples \
    /mnt/c/Users/Luvly/Documents/cse284-final-project/data_raw/chr${chr}_g2.map \
    /mnt/c/Users/Luvly/Documents/cse284-final-project/results/GERMLINE2/chr${chr}_1000_g2_out_hap \
    2> /mnt/c/Users/Luvly/Documents/cse284-final-project/results/GERMLINE2/chr${chr}_1000_g2_hap_runtime.log
done

In haploid mode, output identifiers include .0 or .1 suffixes indicating which haplotype participated in the match.

### Post-processing and Analysis Script

To summarize and compare GERMLINE2 results with other methods, we implemented a Python script:
scripts/analyze_germline2.py

This script parses GERMLINE2 output files and computes pairwise IBD statistics comparable to those produced by PLINK.

Input parsing

The script reads the following inputs:
| Input | Description |
|------|-------------|
| GERMLINE2 diploid output | Segment-level IBD results between individuals |
| GERMLINE2 haploid output | Haplotypic IBD segments used to infer IBD states |
| Genetic maps | Used to estimate chromosome lengths in centimorgans |
| Sample file | Used to construct the list of individuals |

The script extracts sample identifiers from the .samples file and constructs pairwise sample summaries across the three chromosomes.

#### Segment aggregation

For each pair of individuals, the script aggregates GERMLINE2 segment output to compute summary statistics:
| Metric | Description |
|------|-------------|
| Number of shared segments | Count of detected IBD segments |
| Total shared length | Sum of shared segment lengths |
| Maximum segment length | Longest detected segment |
| Mean segment length | Average length of segments |

These aggregated values allow GERMLINE2 results to be summarized at the pair level, making them comparable to PLINK’s pairwise relatedness metrics.

#### IBD state estimation

When haploid mode output is provided, the script estimates the probabilities of the three IBD states:
| State | Interpretation |
|------|---------------|
| P(IBD = 0) | No shared haplotypes |
| P(IBD = 1) | One haplotype shared |
| P(IBD = 2) | Two haplotypes shared |

These probabilities are estimated from haplotype-level segment overlap and normalized by the approximate chromosome lengths obtained from the genetic maps. The script then computes a kinship coefficient:
φ = 0.5 × IBD2 + 0.25 × IBD1

#### Output files

The script produces the following files:
| File | Description |
|------|-------------|
| pairwise_ibd.csv | Pair-level summary statistics from GERMLINE2 segments |
| pairwise_ibd_by_chr.csv | Pair-level summary statistics by chromosome |
| pairwise_ibd012.csv | Estimated IBD state probabilities and kinship |
| segments_summary.csv | Distribution statistics of detected segments |
| segments_summary_by_chr.csv | Segment summaries by chromosome |
| runtime_summary.csv | Parsed runtime and memory summaries |
| runtime_summary_by_mode.csv | Runtime summary grouped by mode |

It also generates figures including segment length and pairwise sharing distributions.

#### Running the analysis script

The analysis was executed from the project root as:
python3 scripts/analyze_germline2.py \
  --prefix results/GERMLINE2/chr20_1000_g2_out \
  --prefix results/GERMLINE2/chr21_1000_g2_out \
  --prefix results/GERMLINE2/chr22_1000_g2_out \
  --hap_prefix results/GERMLINE2/chr20_1000_g2_out_hap \
  --hap_prefix results/GERMLINE2/chr21_1000_g2_out_hap \
  --hap_prefix results/GERMLINE2/chr22_1000_g2_out_hap \
  --map data_raw/chr20_g2.map \
  --map data_raw/chr21_g2.map \
  --map data_raw/chr22_g2.map \
  --sample data_raw/chr22_1000_g2.samples \
  --runtime_log results/GERMLINE2/chr20_1000_g2_runtime.log \
  --runtime_log results/GERMLINE2/chr21_1000_g2_runtime.log \
  --runtime_log results/GERMLINE2/chr22_1000_g2_runtime.log \
  --hap_runtime_log results/GERMLINE2/chr20_1000_g2_hap_runtime.log \
  --hap_runtime_log results/GERMLINE2/chr21_1000_g2_hap_runtime.log \
  --hap_runtime_log results/GERMLINE2/chr22_1000_g2_hap_runtime.log \
  --outdir results
  
| Argument | Description |
|----------|-------------|
| --prefix | GERMLINE2 diploid segment output |
| --hap_prefix | GERMLINE2 haploid output |
| --map | Genetic map file |
| --sample | Sample identifier file |
| --runtime_log | Diploid GERMLINE2 runtime log |
| --hap_runtime_log | Haploid GERMLINE2 runtime log |
| --outdir | Directory for output summaries |

### Results

GERMLINE2 outputs individual IBD segments with the following structure:
| Column | Description |
|-------|-------------|
| ID1 | Individual 1 or haplotype 1 |
| ID2 | Individual 2 or haplotype 2 |
| P0 | Segment start position (bp) |
| P1 | Segment end position (bp) |
| cM | Segment length |
| #words | Number of seed matches |
| #gaps | Allowed mismatches |

Example diploid segment:
HG02356 HG00596 19581331 19844492 1.04501 12 7

Example haploid segment:
NA18999.0 HG03800.0 16793045 17616889 1.0664 32 18

These outputs indicate shared haplotype segments of roughly 1 cM or greater, consistent with the -m 1.0 threshold used in the final runs.

#### Metrics for Benchmarking

GERMLINE2 produces segment-level outputs that can be aggregated to generate pairwise statistics comparable to PLINK.
| Metric | Source | Description |
|-------|-------|-------------|
| Segment count per pair | GERMLINE2 output | Number of detected IBD segments |
| Segment length | GERMLINE2 output | Length of each segment |
| Total shared length | Aggregated segments | Sum of segment lengths |
| P(IBD ≥ 1) | Shared length / chromosome length | Probability of sharing |
| P(IBD = 0,1,2) | Haploid mode segments | Estimated IBD states |

Segment-level statistics also enable comparison with other haplotype-based methods such as Refined IBD.

#### Computational Performance

The final GERMLINE2 runs were performed on the 1000-individual subset across chromosomes 20, 21, and 22. Runtime logs were collected using /usr/bin/time -v for both diploid and haploid runs. These logs were later parsed into runtime_summary.csv and runtime_summary_by_mode.csv to compare elapsed time, CPU usage, and memory consumption across modes.

In general, haploid mode produced substantially more output segments than diploid mode, as expected from the larger number of haplotype-to-haplotype comparisons.

Observed Segment Counts

Using final parameters (-m 1.0, -g 2, -d 0, -f 0.05), the diploid and haploid outputs showed the following chromosome-level counts:
| Chromosome | Diploid segments | Haploid segments |
|-----------|------------------|------------------|
| chr20 | 1296 | 167,823 |
| chr21 | 1064 | 110,666 |
| chr22 | 591 | 70,822 |

These counts are biologically plausible and show the expected pattern of more detected segments in haploid mode than diploid mode.

#### Kinship Classification

Kinship estimates were derived from pairwise_ibd012.csv using the thresholds:
| Kinship coefficient | Relationship |
|---------------------|-------------|
| > 0.177 | 1st degree |
| 0.088–0.177 | 2nd degree |
| 0.044–0.088 | 3rd degree |
| < 0.044 | Unrelated |

A simple classification script was used to assign relationship labels and count pair types. In the analyzed subset, all detected pairs were classified as unrelated, which is consistent with the design of the 1000 Genomes Project, where close relatives are generally excluded.

---

## Method 3: Refined IBD (Beagle)

### Approach

Refined IBD is a haplotype-based method for detecting **identity-by-descent (IBD) segments** using a probabilistic model of haplotype sharing. Unlike PLINK, which estimates relatedness directly from unphased genotype data, Refined IBD operates on **phased haplotypes** and identifies genomic segments shared between individuals due to inheritance from a recent common ancestor.

The Refined IBD pipeline consists of two major steps:

1. **Haplotype phasing** using **Beagle v5**
2. **IBD segment detection** using **Refined IBD**

The resulting segment-level output contains genomic coordinates, estimated segment lengths in centimorgans (cM), and log-odds (LOD) scores representing confidence in the detected IBD segment.

Because Refined IBD reports **explicit segment coordinates**, it enables detailed comparison with other haplotype-based approaches such as **GERMLINE2**, while aggregated pairwise statistics allow comparison with **PLINK’s pairwise relatedness estimates**.

Since no external recombination map was supplied, the default approximation used by both Beagle and Refined IBD was:

```
1 cM ≈ 1 Mb
```

This assumption was used consistently during downstream analysis.

---

### Input Data Preparation

The Refined IBD pipeline was run on the same **1,000 individuals** selected from the 1000 Genomes Project dataset.

To improve signal detection and increase genomic coverage, we analyzed **chromosomes 20, 21, and 22**.

The full analysis pipeline is implemented in:

```
scripts/run_3chr_1000_ibd.sh
```

For each chromosome, the pipeline performs the following preprocessing steps:

1. **Subset samples**

```
bcftools view -S final_samples_1000.txt
```

2. **Normalize variants to biallelic form**

```
bcftools norm -m -both
```

3. **Remove variants containing missing genotypes**

```
bcftools view -g ^miss
```

4. **Remove duplicate markers**

```
bcftools norm -d all
```

Removing duplicate markers is required because Beagle terminates when duplicate variants are present.

After filtering and deduplication, the marker counts used for phasing were:

| Chromosome | Markers   |
| ---------- | --------- |
| chr20      | 1,811,167 |
| chr21      | 1,104,408 |
| chr22      | 1,102,381 |

---

### Phasing with Beagle

Genotypes were phased using **Beagle v5.5**.

Example command:

```bash
java -Xmx8g -jar beagle.27Feb25.75f.jar \
  gt=data_qc/chr22_subset_biallelic_nomiss_dedup.vcf.gz \
  out=results/chr22_phased \
  nthreads=8
```

Equivalent commands were run for chromosomes **20** and **21**.

The resulting phased genotype files were:

```
results/chr20_phased.vcf.gz
results/chr21_phased.vcf.gz
results/chr22_phased.vcf.gz
```

Phasing runtime was relatively short:

| Chromosome | Runtime    |
| ---------- | ---------- |
| chr20      | 51 seconds |
| chr21      | 33 seconds |
| chr22      | 31 seconds |

Total Beagle runtime across all three chromosomes was approximately **2 minutes**.

---

### Running Refined IBD

Refined IBD was then applied to each phased chromosome to detect shared haplotype segments.

Example command:

```bash
java -Xmx16g -jar refined-ibd.17Jan20.102.jar \
  gt=results/chr22_phased.vcf.gz \
  out=results/chr22_refinedibd \
  nthreads=8 \
  lod=1.0 \
  length=0.5
```

#### Parameters

| Parameter  | Description                                                  |
| ---------- | ------------------------------------------------------------ |
| `lod`      | Minimum log-odds score required for reporting an IBD segment |
| `length`   | Minimum segment length threshold (cM)                        |
| `nthreads` | Number of CPU threads                                        |

These thresholds were selected to capture shorter segments expected in a dataset composed largely of unrelated individuals.

---

#### Output Format

Refined IBD outputs segment-level results in compressed files:

```
results/chr20_refinedibd.ibd.gz
results/chr21_refinedibd.ibd.gz
results/chr22_refinedibd.ibd.gz
```

Each row describes a detected IBD segment:

| Column     | Description                 |
| ---------- | --------------------------- |
| Sample1    | First individual            |
| Hap1       | Haplotype index             |
| Sample2    | Second individual           |
| Hap2       | Haplotype index             |
| Chromosome | Chromosome number           |
| Start      | Segment start position (bp) |
| End        | Segment end position (bp)   |
| LOD        | Log-odds score              |
| Length     | Segment length (cM)         |

Example segment:

```
NA18923 1 NA18507 2 22 28470074 29055755 13.28 0.586
```

This indicates a shared haplotype segment of approximately **0.586 cM** on chromosome 22.

The per-chromosome outputs were merged into a single dataset:

```
results/refinedibd_3chr_all.ibd.gz
```

Additional summary files generated by the pipeline include:

```
results/refinedibd_segment_counts.txt
results/refinedibd_3chr_pairwise_summary.tsv
```

---

### Post-processing and Analysis

To summarize Refined IBD results and enable comparison with other methods, we implemented a Python analysis script:

```
scripts/analyze_refinedibd.py
```

The script produces several outputs:

```
results/refinedibd_pairwise_ibd.csv
results/refinedibd_pairwise_ibd012.csv
results/refinedibd_segments_summary.csv
refinedibd_segment_lengths.png
refinedibd_shared_cm_hist.png
```

These files summarize segment-level and pair-level IBD statistics.

---

#### Segment-Level Results

Across chromosomes 20–22, Refined IBD detected:

| Chromosome | Segments |
| ---------- | -------- |
| chr20      | 380      |
| chr21      | 112      |
| chr22      | 90       |
| **Total**  | **582**  |

The majority of segments were close to the minimum detection threshold of **0.5 cM**, as shown in the segment length distribution.

Observed segment lengths ranged approximately from:

```
0.50 cM – 0.85 cM
```

with most segments clustered between **0.50 and 0.60 cM**.

This pattern is expected because:

* the minimum segment threshold was set to **0.5 cM**
* the dataset contains mostly **unrelated individuals**
* only **three chromosomes** were analyzed.

---

#### Pairwise IBD Sharing

After aggregating segments by pair of individuals:

* **528 pairs** shared at least one IBD segment
* Mean number of segments per pair: **~1.1**
* Maximum segments for a single pair: **11**

Total shared genomic length per pair ranged approximately from:

```
0.5 cM – 2.3 cM
```

Most pairs shared **less than 1 cM** of total IBD.

The distribution of shared genomic length is shown in the generated histogram:

```
refinedibd_shared_cm_hist.png
```

---

### Computational Performance

Refined IBD required substantially more computation than PLINK due to haplotype modeling.

#### Beagle Phasing

| Chromosome | Runtime    |
| ---------- | ---------- |
| chr20      | 51 seconds |
| chr21      | 33 seconds |
| chr22      | 31 seconds |

#### Refined IBD Detection

| Chromosome | Runtime               |
| ---------- | --------------------- |
| chr20      | 17 minutes 37 seconds |
| chr21      | 9 minutes 43 seconds  |
| chr22      | 9 minutes 35 seconds  |

Total Refined IBD runtime across the three chromosomes was approximately **37 minutes**.

Most runtime was spent constructing the haplotype model used during IBD inference.

---

### Metrics for Benchmarking

Refined IBD provides both **segment-level** and **pair-level** statistics that can be compared with other methods.

| Metric                      | Description                       |
| --------------------------- | --------------------------------- |
| Segment count per pair      | Number of detected IBD segments   |
| Segment length distribution | Distribution of segment sizes     |
| Total shared genomic length | Sum of segment lengths per pair   |
| Maximum segment length      | Longest shared segment            |
| P(IBD ≥1)                   | Fraction of chromosome shared     |
| P(IBD = 0,1,2)              | Estimated IBD state probabilities |
| Kinship estimate            | Derived from IBD probabilities    |

These metrics enable:

* **segment-level comparison with GERMLINE2**
* **pairwise relatedness comparison with PLINK**
* **runtime and scalability benchmarking across methods**

---

## Benchmarking Plots

Five comparison plots are generated to evaluate agreement and differences across the three IBD detection methods. All scripts are in `scripts/` and output figures to `plots/`.

### How to Run

```bash
cd cse284-final-project
python3 scripts/plot1_upset_overlap.py
python3 scripts/plot2_z0_z1_scatter.py
python3 scripts/plot3_segment_length.py
python3 scripts/plot4_kinship_scatter.py
python3 scripts/plot5_jaccard.py
```

### Plot Descriptions

| Plot | Script | Description |
|------|--------|-------------|
| **Plot 1** – UpSet Overlap | `plot1_upset_overlap.py` | UpSet plot showing which IBD pairs are detected by each method and their overlaps. |
| **Plot 2** – Z0/Z1 & Shared cM | `plot2_z0_z1_scatter.py` | Two-panel: PLINK Z0 vs Z1 scatter with relationship reference points (left), GERMLINE2 vs Refined IBD total shared cM histogram (right). |
| **Plot 3** – Segment Length Distribution | `plot3_segment_length.py` | Overlaid histograms of IBD segment lengths (cM and Mb) for GERMLINE2 and Refined IBD. PLINK does not produce segment-level data. |
| **Plot 4** – Kinship Scatter | `plot4_kinship_scatter.py` | Three-panel pairwise kinship comparison: PLINK vs GERMLINE2, PLINK vs Refined IBD, GERMLINE2 vs Refined IBD. |
| **Plot 5** – Segment Overlap (Jaccard) | `plot5_jaccard.py` | Per-chromosome IBD coverage bar chart, segment position map on chr20, and Jaccard similarity summary for GERMLINE2 vs Refined IBD. |

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

GERMLINE2 requires phased haplotypes in **SHAPEIT/IMPUTE `.haps` format**, along with a `.sample` file containing sample identifiers and a genetic map specifying marker positions and genetic distances.

The starting dataset consisted of chromosome 22 genotype data derived from the 1000 Genomes Project samples after quality control.

#### Generating the no-missing VCF

GERMLINE2 requires haplotype input without missing genotype calls. To satisfy this requirement, we filtered our dataset to remove any variants containing missing genotypes across the samples.

Starting from the filtered dataset:

`chr22_subset_biallelic.vcf`

which contains only biallelic SNPs from chromosome 22, we removed sites with missing genotype calls using **bcftools**.
bcftools view -g ^miss \
chr22_subset_biallelic.vcf \
-o chr22_subset_biallelic_nomiss.vcf

This filtering step removes any variants with missing genotype values in any sample.

The resulting file `chr22_subset_biallelic_nomiss.vcf` contains **1,110,240 SNPs across 230 individuals**, with complete genotype information for every marker.

---

#### Converting the VCF to haplotype format

The filtered VCF was converted to SHAPEIT-style haplotype and sample files using `bcftools convert`:
bcftools convert --hapsample chr22_g2 chr22_subset_biallelic_nomiss.vcf

This produced two files:

| File | Description |
|-----|-------------|
| `chr22_g2.haps` | Phased haplotypes |
| `chr22_g2.sample` | Sample identifiers |

---

#### Constructing the genetic map

GERMLINE2 requires a genetic map with three columns:
physical_position cm/Mb cumulative_cM

Because a recombination map was used during initial testing, we approximated genetic distance using physical position in megabases:
**cM ≈ physical_position / 1,000,000**

The map was generated using:
awk '{pos=$3; cm=pos/1000000.0; print pos,1,cm}' chr22_g2.haps > chr22_g2.map

The resulting map spans approximately **35.2 cM** across chromosome 22.

---

### Commands

GERMLINE2 was run in two configurations.

#### Diploid segment detection

This run identifies IBD segments between individuals without distinguishing haplotypes.
./g2 -m 1.0 -g 2 -d 1 -f 0.05 \
data_raw/chr22_g2.haps \
data_raw/chr22_g2.sample \
data_raw/chr22_g2.map \
data_raw/chr22_g2_out 

Parameter descriptions:

| Parameter | Description |
|-----------|-------------|
| `-m` | Minimum segment length (cM) |
| `-g` | Maximum allowed gaps |
| `-d` | Dynamic seed threshold |
| `-f` | Minor allele frequency filter |

---

#### Haploid mode (IBD state estimation)

To estimate **IBD state probabilities**, we also ran GERMLINE2 in haploid mode, which restricts matches to the same haplotype:
./g2 -h -m 0.2 -g 6 -d 1 -f 0.05 \
data_raw/chr22_g2.haps \
data_raw/chr22_g2.sample \
data_raw/chr22_g2.map \
data_raw/chr22_g2_out_hap

In haploid mode, output identifiers include `.0` or `.1` suffixes indicating which haplotype participated in the match.

---

### Post-processing and Analysis Script

To summarize and compare GERMLINE2 results with other methods, we implemented a Python script:

`analyze_germline2.py`

This script parses GERMLINE2 output files and computes **pairwise IBD statistics** comparable to those produced by PLINK.

---

#### Input parsing

The script reads the following inputs:

| Input | Description |
|------|-------------|
| GERMLINE2 diploid output | Segment-level IBD results between individuals |
| GERMLINE2 haploid output | Haplotypic IBD segments used to infer IBD states |
| Genetic map | Used to estimate chromosome length in centimorgans |
| Sample file | Used to construct the list of individuals |

The script extracts sample identifiers from the `.sample` file and constructs all possible **pairwise sample combinations**.

---

#### Segment aggregation

For each pair of individuals, the script aggregates GERMLINE2 segment output to compute summary statistics:

| Metric | Description |
|------|-------------|
| Number of shared segments | Count of detected IBD segments |
| Total shared length | Sum of shared segment lengths |
| Maximum segment length | Longest detected segment |
| Mean segment length | Average length of segments |

These aggregated values allow GERMLINE2 results to be summarized at the **pair level**, making them comparable to PLINK’s pairwise relatedness metrics.

---

#### IBD state estimation

When haploid mode output is provided, the script estimates the probabilities of the three IBD states:

| State | Interpretation |
|------|---------------|
| **P(IBD = 0)** | No shared haplotypes |
| **P(IBD = 1)** | One haplotype shared |
| **P(IBD = 2)** | Two haplotypes shared |

These probabilities are estimated from haplotype-level segment overlap and normalized by the approximate chromosome length obtained from the genetic map.

---

#### Output files

The script produces the following files:

| File | Description |
|-----|-------------|
| `pairwise_ibd.csv` | Pair-level summary statistics from GERMLINE2 segments |
| `pairwise_ibd012.csv` | Estimated IBD state probabilities |
| `segment_summary.csv` | Distribution statistics of detected segments |

These outputs provide a consistent format for comparing GERMLINE2 results with PLINK and Refined IBD.

---

#### Running the analysis script

The analysis can be executed using:
python3 analyze_germline2.py \
--prefix data_raw/chr22_g2_out \
--hap_prefix data_raw/chr22_g2_out_hap \
--map data_raw/chr22_g2.map \
--sample data_raw/chr22_g2.sample \
--outdir results

| Argument | Description |
|----------|-------------|
| `--prefix` | GERMLINE2 diploid segment output |
| `--hap_prefix` | GERMLINE2 haploid output |
| `--map` | Genetic map file |
| `--sample` | Sample identifier file |
| `--outdir` | Directory for output summaries |

---

### Results

GERMLINE2 outputs individual **IBD segments** with the following structure:

| Column | Description |
|------|-------------|
| ID1 | Individual 1 |
| ID2 | Individual 2 |
| P0 | Segment start position (bp) |
| P1 | Segment end position (bp) |
| cM | Segment length |
| #words | Number of seed matches |
| #gaps | Allowed mismatches |

Example segment:
NA18632.0 NA18530.1 28342752 29035817 0.693 264 10
This indicates a shared haplotype segment approximately **0.69 Mb long** on chromosome 22.

---

### Metrics for Benchmarking

GERMLINE2 produces segment-level outputs that can be aggregated to generate pairwise statistics comparable to PLINK.

| Metric | Source | Description |
|------|-------|-------------|
| Segment count per pair | GERMLINE2 output | Number of detected IBD segments |
| Segment length | GERMLINE2 output | Length of each segment |
| Total shared length | Aggregated segments | Sum of segment lengths |
| P(IBD ≥1) | Shared length / chromosome length | Probability of sharing |
| P(IBD = 0,1,2) | Haploid mode segments | Estimated IBD states |

Segment-level statistics also enable comparisons with other haplotype-based methods such as **Refined IBD**.

---

### Computational Performance

GERMLINE2 processed the chromosome 22 dataset efficiently despite the large number of markers.

| Analysis | Runtime | Input Size |
|---------|--------|-----------|
| GERMLINE2 diploid detection | ~11 seconds | 1.11M SNPs |
| GERMLINE2 haploid detection | ~9 seconds | 1.11M SNPs |

The hash-based seed-and-extend algorithm enables fast scanning of large haplotype datasets.

---

### Upcoming Action Items

#### Incorporating a recombination map

In the current analysis, genetic distance was approximated using physical position in megabases. As a next step, we plan to obtain a **high-resolution recombination map for chromosome 22** (e.g., from the HapMap or 1000 Genomes genetic maps) so that segment lengths can be measured using accurate centimorgan distances rather than approximate physical distances.

Using a proper recombination map should improve the biological accuracy of detected segment lengths and provide more reliable IBD state estimates.

---

#### Parameter tuning

GERMLINE2 performance depends strongly on parameter selection. Future work will explore:

- **Minimum match length (-m)**  
  Evaluate thresholds between **0.1–2.0 Mb**

- **Allowed gaps (-g)**  
  Test values between **2–8**

- **Minor allele frequency filter (-f)**  
  Evaluate filtering of rare variants

The goal is to balance **sensitivity and specificity** of IBD detection.

---

#### Expanding genomic scope

Current results use only **chromosome 22**. Future analyses may include:

- Additional chromosomes  
- Whole-genome analysis  
- Cross-chromosome consistency checks  

---

#### Increasing sample coverage

The dataset currently includes **230 individuals**. Future analyses may evaluate:

- Additional individuals from the 1000 Genomes dataset  
- Population-specific subsets  
- Effects of sample size on IBD detection sensitivity  

---

#### Cross-method comparison

The final stage will compare results from three approaches:

| Method | Type |
|------|------|
| **PLINK** | Allele-frequency relatedness estimation |
| **GERMLINE2** | Haplotype segment detection |
| **Refined IBD** | Probabilistic haplotype inference |

Comparison metrics will include:

- Pairwise relatedness estimates  
- Segment length distributions  
- Total shared genomic length  
- Agreement between methods  
- Computational runtime and scalability

---

## Method 3: Refined IBD (Beagle)

## Approach

Refined IBD is a haplotype-based method for detecting **identity-by-descent (IBD) segments** using a probabilistic model of haplotype sharing. The method evaluates whether two haplotypes share genomic regions inherited from a recent common ancestor by computing likelihood scores for candidate shared segments.

Because Refined IBD operates on **phased haplotypes**, genotype data must first be phased before IBD detection can be performed. In this project, phasing was performed using **Beagle**, which infers haplotype phase using localized haplotype clustering.

After phasing, Refined IBD scans the phased haplotypes and identifies shared segments between individuals. Each detected segment includes genomic coordinates, a log-odds (LOD) score representing confidence in the IBD call, and an estimated segment length in centimorgans (cM).

This segment-level output enables direct comparison with other haplotype-based methods such as **GERMLINE2**, while aggregated pairwise statistics allow comparison with **PLINK’s pairwise relatedness estimates**.

---

## Phasing with Beagle

Genotypes were phased using **Beagle v5**.

The chromosome 22 VCF containing the selected 1000 Genomes samples was used as input.

```bash
java -Xmx8g -jar beagle.27Feb25.75f.jar \
  gt=data_raw/chr22_subset.vcf.gz \
  out=results/beagle/chr22_phased \
  nthreads=8
```

This produced the phased VCF file:

```
results/beagle/chr22_phased.vcf.gz
```

Phased genotypes use the `|` separator to indicate haplotype phase.

Example genotype:

```
0|1
```

This indicates that the reference allele occurs on one haplotype and the alternate allele on the other.

---

## Running Refined IBD

IBD segments were detected using the Refined IBD implementation by Browning & Browning.

```bash
java -Xmx16g -jar refined-ibd.17Jan20.102.jar \
  gt=results/beagle/chr22_phased.vcf.gz \
  out=results/refinedibd/chr22_refinedibd \
  nthreads=8 \
  lod=1.0 \
  length=0.5
```

### Parameters

| Parameter  | Description                                                  |
| ---------- | ------------------------------------------------------------ |
| `lod`      | Minimum log-odds score required for reporting an IBD segment |
| `length`   | Minimum segment length in centimorgans                       |
| `nthreads` | Number of CPU threads used                                   |

Because no recombination map was provided, Refined IBD used the default approximation:

```
1 cM ≈ 1 Mb
```

This assumption was used consistently in downstream analysis.

---

## Output Format

Refined IBD produces detected segments in the compressed file:

```
results/refinedibd/chr22_refinedibd.ibd.gz
```

Each row represents a detected IBD segment and contains:

| Column     | Description                   |
| ---------- | ----------------------------- |
| Sample1    | First individual              |
| Hap1       | Haplotype index               |
| Sample2    | Second individual             |
| Hap2       | Haplotype index               |
| Chromosome | Chromosome number             |
| Start      | Segment start position (bp)   |
| End        | Segment end position (bp)     |
| LOD        | Log-odds score supporting IBD |
| Length     | Segment length (cM)           |

Example output:

```
NA18923 1 NA18507 2 22 28470074 29055755 13.28 0.586
```

This indicates a shared haplotype segment approximately **0.586 cM long** on chromosome 22.

---

## Post-processing and Analysis Script

To summarize Refined IBD results and make them comparable with GERMLINE2 and PLINK outputs, we implemented a Python analysis script:

```
scripts/analyze_refinedibd.py
```

The script parses the Refined IBD segment file and computes both **segment-level statistics** and **pairwise relatedness metrics**.

---

## Segment Aggregation

First, the script aggregates detected segments across all individuals to compute global summary statistics.

Computed metrics include:

| Metric                           | Description                     |
| -------------------------------- | ------------------------------- |
| Number of segments               | Total detected IBD segments     |
| Segment length distribution      | Distribution of segment lengths |
| Mean / median segment length     | Average segment size            |
| Minimum / maximum segment length | Range of detected segments      |
| Mean / maximum LOD score         | Confidence of detected segments |

The script also generates a histogram of segment lengths:

```
refinedibd_segment_lengths.png
```

---

## Pairwise IBD Aggregation

To enable comparison with other methods, segment-level results are aggregated at the **pair level**.

For each pair of individuals the script computes:

| Metric                 | Description                    |
| ---------------------- | ------------------------------ |
| Number of segments     | Count of shared segments       |
| Total shared length    | Sum of segment lengths (cM)    |
| Maximum segment length | Longest detected segment       |
| Mean segment length    | Average segment size           |
| Mean / maximum LOD     | Average confidence of segments |

These statistics are written to:

```
results/refinedibd_pairwise_ibd.csv
```

This file provides a consistent pairwise format that can be compared directly with:

* GERMLINE2 aggregated segment statistics
* PLINK pairwise relatedness estimates

---

## IBD State Estimation

To approximate IBD state probabilities, the analysis script uses the haplotype indices provided in the Refined IBD output.

Segments are grouped into two haplotype coverage sets:

* segments involving **haplotype 0**
* segments involving **haplotype 1**

Using these intervals, the script estimates the following probabilities:

| State          | Interpretation        |
| -------------- | --------------------- |
| **P(IBD = 0)** | No shared haplotypes  |
| **P(IBD = 1)** | One haplotype shared  |
| **P(IBD = 2)** | Two haplotypes shared |

The estimates are computed using genomic coverage in centimorgan space:

```
P(IBD ≥1) = union(shared intervals) / chromosome_length
P(IBD2) = overlap(hap0_intervals, hap1_intervals) / chromosome_length
P(IBD1) = P(IBD ≥1) − P(IBD2)
P(IBD0) = 1 − P(IBD ≥1)
```

A kinship estimate can also be derived:

```
φ = 0.5 × P(IBD2) + 0.25 × P(IBD1)
```

These results are written to:

```
results/refinedibd_pairwise_ibd012.csv
```

This format is comparable to:

* **PLINK Z0/Z1/Z2 probabilities**
* **GERMLINE2 haplotype-based IBD state estimates**

---

## Results

Refined IBD detected a small number of IBD segments among the analyzed individuals.

| Metric                  | Value            |
| ----------------------- | ---------------- |
| Total segments detected | 10               |
| Segment length range    | 0.509 – 0.692 cM |
| Mean segment length     | ~0.59 cM         |
| Number of samples       | 230              |

The small number of detected segments is expected because the 1000 Genomes samples consist primarily of **unrelated individuals**, and the analysis was restricted to **a single chromosome**.

---

## Computational Performance

Refined IBD processed the chromosome 22 dataset efficiently.

| Step                  | Runtime     | Threads |
| --------------------- | ----------- | ------- |
| Beagle phasing        | ~10 seconds | 8       |
| Refined IBD detection | ~1 minute   | 8       |

Most of the runtime was spent constructing the haplotype model used during IBD inference.

---

## Metrics for Benchmarking

Refined IBD provides both **segment-level** and **pair-level** metrics for benchmarking.

| Metric                      | Description                       |
| --------------------------- | --------------------------------- |
| Segment count per pair      | Number of detected IBD segments   |
| Segment length distribution | Distribution of segment sizes     |
| Total shared genomic length | Sum of segment lengths            |
| P(IBD ≥1)                   | Fraction of chromosome shared     |
| P(IBD = 0,1,2)              | Estimated IBD state probabilities |
| Kinship estimate            | Derived from IBD probabilities    |

These metrics allow direct comparison with:

| Method          | Output type                       |
| --------------- | --------------------------------- |
| **PLINK**       | Pairwise relatedness coefficients |
| **GERMLINE2**   | Segment-level haplotype matches   |
| **Refined IBD** | Probabilistic haplotype segments  |

Together, these methods provide complementary perspectives on genetic relatedness and enable benchmarking of **allele-frequency–based vs haplotype-based IBD detection approaches**.

# CSE 284 Final Project — Reproducibility Guide (PLINK, GERMLINE2, Refined IBD)

This README is a run guide for reproducing the analysis. It focuses on:
- how to install dependencies,
- how to run each method,
- where outputs are written,
- and quick sanity checks to verify your run matches ours.

It intentionally avoids extended result interpretation (that is in the report).

---

## 1) Dataset used in this project

All analyses in this repository are run on:
- **1000 individuals** listed in:
  - `data_raw/final_samples_1000.txt`
- **chromosomes 20, 21, 22** only
- 1000 Genomes Phase 3 VCFs:
  - `ALL.chr20.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
  - `ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
  - `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`

> Note: `scripts/run_3chr_1000.sh` now reuses `data_raw/final_samples_1000.txt` if it exists (for reproducibility).

---

## 2) Project structure

```text
cse284-final-project/
  data_raw/                 input files and sample list
  data_qc/                  subset/QC intermediate files
  results/
    king_3chr_1000.kin0     PLINK KING output
    genome_3chr_1000.genome PLINK --genome output
    refinedibd_3chr_all.ibd.gz
    refinedibd_3chr_pairwise_summary.tsv
    GERMLINE2/
      chr*_1000_g2_out
      chr*_1000_g2_out_hap
      pairwise_ibd.csv
      pairwise_ibd012.csv
  scripts/
    run_3chr_1000.sh
    run_3chr_1000_ibd.sh
    analyze_germline2.py
    analyze_refinedibd.py
    plot1_upset_overlap.py
    plot2_z0_z1_scatter.py
    plot3_segment_length.py
    plot4_kinship_scatter.py
    plot5_jaccard.py
  plots/                    generated benchmark figures
```

### Main output files by method

- **PLINK (KING + --genome)**
  - `results/king_3chr_1000.kin0`: pairwise KING kinship estimates (`KINSHIP`) used for related-pair thresholding.
  - `results/genome_3chr_1000.genome`: pairwise `Z0/Z1/Z2` and `PI_HAT` values from `--genome`.
  - `results/relatives_detected.txt`: filtered list of pairs reported as related by the pipeline.

- **GERMLINE2**
  - `results/GERMLINE2/chr*_1000_g2_out`: diploid segment calls per chromosome.
  - `results/GERMLINE2/chr*_1000_g2_out_hap`: haplotype-mode segment calls per chromosome.
  - `results/GERMLINE2/pairwise_ibd.csv`: pair-level totals (segment counts and shared cM).
  - `results/GERMLINE2/pairwise_ibd012.csv`: pair-level IBD state summary (`P_IBD0/P_IBD1/P_IBD2`) and kinship proxy.

- **Beagle + Refined IBD**
  - `results/chr*_refinedibd.ibd.gz`: per-chromosome Refined IBD segment calls.
  - `results/refinedibd_3chr_all.ibd.gz`: merged segment calls across chr20/21/22.
  - `results/refinedibd_3chr_pairwise_summary.tsv`: pairwise summary table used by downstream comparisons/plots.

---

## 3) Install dependencies

You can use any equivalent setup; below is one reproducible option.

### 3.1 Core tools

- `bcftools`
- `plink2`
- `plink` (1.9)
- `python3` + `pip`
- Java (for Beagle / Refined IBD)

### 3.2 Python packages

```bash
python3 -m pip install pandas numpy matplotlib upsetplot scipy
```

### 3.3 External binaries/jars

- GERMLINE2 executable (`g2`) from: https://github.com/gusevlab/germline2
- Beagle jar: `beagle.27Feb25.75f.jar`
- Refined IBD jar: `refined-ibd.17Jan20.102.jar`

Place Beagle/Refined-IBD jars in the project root (same level as `scripts/`).

---

## 4) Run Method 1: PLINK (KING + --genome)

From project root:

```bash
bash scripts/run_3chr_1000.sh
```

This pipeline:
1. Downloads chr20/21/22 VCFs + panel if missing
2. Uses `data_raw/final_samples_1000.txt`
3. Subsets + normalizes VCFs
4. Builds QC PLINK files
5. Runs KING and PLINK `--genome`

### Main outputs

- `results/king_3chr_1000.kin0`
- `results/genome_3chr_1000.genome`
- `results/relatives_detected.txt`

### Quick sanity checks

```bash
# Related pairs by KING threshold > 0.0442
awk 'NR>1 && $8>0.0442{n++} END{print n+0}' results/king_3chr_1000.kin0
# expected: 323

# Top --genome pair should include NA20320 / NA20321 with PI_HAT ~0.5049
grep -E 'NA20320.*NA20321|NA20321.*NA20320' results/genome_3chr_1000.genome
```

---

## 5) Run Method 2: GERMLINE2

GERMLINE2 was run in diploid + haploid mode, then post-processed for pair-level summaries.

### 5.1 Prepare GERMLINE2 input format (if re-running from VCF)

```bash
bcftools convert --hapsample data_raw/chr20_1000_g2 data_raw/chr20_1000.vcf.gz
bcftools convert --hapsample data_raw/chr21_1000_g2 data_raw/chr21_1000.vcf.gz
bcftools convert --hapsample data_raw/chr22_1000_g2 data_raw/chr22_1000.vcf.gz

gunzip -f data_raw/chr20_1000_g2.hap.gz
gunzip -f data_raw/chr21_1000_g2.hap.gz
gunzip -f data_raw/chr22_1000_g2.hap.gz
```

Map files used:
- `data_raw/chr20_g2.map`
- `data_raw/chr21_g2.map`
- `data_raw/chr22_g2.map`

### 5.2 Run GERMLINE2 (diploid + haploid)

Assuming `g2` is available in your shell path:

```bash
for chr in 20 21 22; do
  g2 -m 1.0 -g 2 -d 0 -f 0.05 \
    data_raw/chr${chr}_1000_g2.hap \
    data_raw/chr${chr}_1000_g2.samples \
    data_raw/chr${chr}_g2.map \
    results/GERMLINE2/chr${chr}_1000_g2_out

done

for chr in 20 21 22; do
  g2 -h -m 1.0 -g 2 -d 0 -f 0.05 \
    data_raw/chr${chr}_1000_g2.hap \
    data_raw/chr${chr}_1000_g2.samples \
    data_raw/chr${chr}_g2.map \
    results/GERMLINE2/chr${chr}_1000_g2_out_hap

done
```

### 5.3 Post-process GERMLINE2 outputs

```bash
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
  --outdir results/GERMLINE2
```

### Main outputs

- `results/GERMLINE2/pairwise_ibd.csv`
- `results/GERMLINE2/pairwise_ibd012.csv`
- `results/GERMLINE2/chr*_1000_g2_out`
- `results/GERMLINE2/chr*_1000_g2_out_hap`

### Quick sanity checks

```bash
awk 'END{print NR-1}' results/GERMLINE2/pairwise_ibd.csv
# expected: 2905

cat results/GERMLINE2/chr*_1000_g2_out | wc -l
# expected: 2951
```

---

## 6) Run Method 3: Beagle + Refined IBD

From project root:

```bash
bash scripts/run_3chr_1000_ibd.sh
```

This pipeline:
1. Uses `data_raw/final_samples_1000.txt`
2. Subsets + normalizes chr20/21/22
3. Removes missing-genotype sites + duplicate markers
4. Phases with Beagle
5. Runs Refined IBD (`lod=1.0`, `length=0.5`)
6. Produces merged segment file + pairwise summary

### Main outputs

- `results/chr20_refinedibd.ibd.gz`
- `results/chr21_refinedibd.ibd.gz`
- `results/chr22_refinedibd.ibd.gz`
- `results/refinedibd_3chr_all.ibd.gz`
- `results/refinedibd_3chr_pairwise_summary.tsv`

### Quick sanity checks

```bash
awk 'END{print NR-1}' results/refinedibd_3chr_pairwise_summary.tsv
# expected: 524

gzip -dc results/refinedibd_3chr_all.ibd.gz | wc -l
# expected: 545

# sample consistency check (should be 0 outside final_samples_1000.txt)
comm -23 <(gzip -dc results/refinedibd_3chr_all.ibd.gz | awk '{print $1; print $3}' | sort -u) <(sort data_raw/final_samples_1000.txt) | wc -l
# expected: 0
```

---

## 7) Metrics used for benchmarking

### PLINK
- `KINSHIP` (from KING)
- `Z0`, `Z1`, `Z2`, `PI_HAT` (from `--genome`)

### GERMLINE2
- segment counts and lengths
- per-pair total shared cM
- derived IBD state probabilities and kinship proxy (`pairwise_ibd012.csv`)

### Refined IBD
- segment counts and lengths
- per-pair total shared cM
- per-pair summary from `refinedibd_3chr_pairwise_summary.tsv`

---

## 8) Generate benchmark figures

From project root:

```bash
python3 scripts/plot1_upset_overlap.py
python3 scripts/plot2_z0_z1_scatter.py
python3 scripts/plot3_segment_length.py
python3 scripts/plot4_kinship_scatter.py
python3 scripts/plot5_jaccard.py
```

Output files:
- `plots/plot1_upset_overlap.png`
- `plots/plot2_z0_z1_scatter.png`
- `plots/plot3_segment_length.png`
- `plots/plot4_kinship_scatter.png`
- `plots/plot5_jaccard.png`

---

## 9) Quick verification checklist

After running pipelines/scripts, verify:

1. `final_samples_1000.txt` is used (not a newly randomized cohort)
2. PLINK KING related pair count = **323**
3. GERMLINE2 pair count (`pairwise_ibd.csv`) = **2905**
4. Refined IBD pair count (`refinedibd_3chr_pairwise_summary.tsv`) = **524**
5. Refined IBD total segments (`refinedibd_3chr_all.ibd.gz`) = **545**
6. All 5 plots are generated under `plots/`

If these checks match, your run is consistent with the repository’s reference outputs.

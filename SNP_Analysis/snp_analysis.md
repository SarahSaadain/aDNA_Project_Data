# aDNA Consensus Sequence and SNP Calling

This part of the pipeline processes ancient DNA (aDNA) sequencing data to:

1. **Call SNPs**
2. **Generate consensus sequences**
3. **Output summary files for downstream analyses (e.g., haplotypes, MAF, visualization)**

This implementation is tailored for **aDNA**, taking into account damage, low coverage, and short fragment lengths.


## Overview

We use **ANGSD** (Analysis of Next Generation Sequencing Data) to create consensus sequences and call SNPs in one streamlined process.

### Script: `execute_angsd_create_consensus_sequence`

For each sorted BAM file:

* A consensus FASTA sequence is generated.
* Minor allele frequencies (MAF) are computed.
* SNPs are identified.

### **Parallel Execution**

Multiple BAM files are processed in parallel using Python’s `multiprocessing.Pool`.

## ANGSD Settings and Rationale

```bash
-angsd \
  -out <output> \
  -i <sorted.bam> \
  -ref <reference.fa> \
  -doFasta 2 \
  -doCounts 1 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -GL 2 \
  -SNP_pval 1e-6 \
  -minMapQ 30 \
  -minQ 20 \
  -uniqueOnly 1 \
  -remove_bads 1 \
  -baq 1 \
  -C 50
```

### ✅ Explanation of Key Parameters

| Option            | Purpose                                              | Rationale (aDNA-specific)                                            |
| ----------------- | ---------------------------------------------------- | -------------------------------------------------------------------- |
| `-doFasta 2`      | Emit consensus using EBD model                       | EBD (Extended Basecall Distribution) is robust for low-coverage aDNA |
| `-doCounts 1`     | Count per-base nucleotide frequencies                | Required for allele inference                                        |
| `-doMajorMinor 1` | Infer major and minor alleles                        | Needed for MAF and consensus                                         |
| `-doMaf 1`        | Output minor allele frequencies                      | For SNP analysis                                                     |
| `-GL 2`           | Use SAMtools genotype likelihood model               | Good default for aDNA with damage                                    |
| `-SNP_pval 1e-6`  | P-value threshold for SNP calling                    | Conservative filtering to reduce false positives in damaged DNA      |
| `-minMapQ 30`     | Only consider reads with mapping quality ≥ 30        | Excludes poorly mapped reads common in aDNA                          |
| `-minQ 20`        | Only consider base calls with quality ≥ 20           | Filters low-confidence bases                                         |
| `-uniqueOnly 1`   | Exclude multi-mapping reads                          | Ensures specificity in mapping                                       |
| `-remove_bads 1`  | Remove reads flagged as "bad"                        | Further improves data integrity                                      |
| `-baq 1`          | Adjust base qualities around indels                  | Mitigates false SNPs at indel boundaries                             |
| `-C 50`           | Recalibrate mapping quality for excessive mismatches | Critical for aDNA which often has mismatches due to degradation      |


## Output Files

For each input BAM file:

* `<sample>_consensus.fa.gz`: Consensus sequence (compressed FASTA)
* `<sample>.mafs.gz`: SNPs and allele frequencies
* `<sample>.saf.idx` (if extended): For downstream selection or diversity analysis

## Downstream Use Cases

* Haplotype inference (via consensus)
* Population diversity (π, Tajima’s D)
* Introgression and admixture analyses
* aDNA authenticity checks (e.g., excess C>T transitions)

## Notes

* BAM files must be sorted and indexed.
* The reference genome (`-ref`) must be indexed (`.fai`).
* Ensure adapter trimming and deduplication is completed **before** running ANGSD.

## SNP and Minor Allele Frequency (MAF) Visualizations

This section describes the plots automatically generated from ANGSD `.mafs.gz` files, used to assess allele frequency distributions and reference allele concordance in ancient DNA (aDNA) data.

These visualizations were generated using the `plot_snp_by_individual.R` script.

### 1. **MAF Density Distribution** (`*_maf_density.png`)

* **Description**: Kernel density plot of minor allele frequencies (`knownEM` column).
* **Purpose**: Visualizes how common or rare minor alleles are across all polymorphic sites.
* **Interpretation**:

  * A peak near 0 or 0.5 indicates many low-frequency or heterozygous variants.
  * Skewed distributions may reflect damage, contamination, or demographic events.

### 2. **Histogram of MAF Values** (`*_maf_histogram.png`)

* **Description**: Histogram showing the number of SNPs in each MAF bin.
* **Purpose**: Complements the density plot by giving raw counts of minor allele frequencies.
* **Interpretation**:

  * Useful for spotting discrete peaks or gaps in allele frequency distributions.
  * Helps differentiate true variation from sequencing noise or post-mortem damage.

### 3. **SNP Density per Scaffold** (`*_snp_density_per_scaffold.png`)

* **Description**: Barplot of the number of SNPs per scaffold.
* **Purpose**: Identifies uneven variant distribution across the genome.
* **Interpretation**:

  * High-density scaffolds may represent better-covered or more variable regions.
  * Scaffolds with no SNPs may be poorly covered or invariant.

### 4. **Reference vs. Major Allele Heatmap** (`*_ref_vs_major.png`)

* **Description**: Heatmap showing counts of combinations between reference bases and inferred major alleles.
* **Purpose**: Detects divergence between the sample and the reference genome.
* **Interpretation**:

  * Cells along the diagonal (ref = major) suggest high similarity or reference bias.
  * Off-diagonal counts reveal fixed differences or strong minor alleles becoming major in the sample.
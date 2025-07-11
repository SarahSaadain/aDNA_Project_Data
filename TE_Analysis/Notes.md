# Working Directory

/mnt/data5/sarah/TE_Analysis

# Check candidates vs NCBI and DFAM Database/API

Use: 
* [run_ncbi_dfam_analysis.sh](run_ncbi_dfam_analysis.sh) to call
* [heck_fasta_against_ncbi.py](check_fasta_against_ncbi.py)
* [check_fasta_against_dfam.py](check_fasta_against_dfam.py)

```bash
nohup bash run_ncbi_dfam_analysis.sh > ncbi_dfam_analysis.log 2>&1 &
```

# Analysis of gaps by species accross individuals

To check mathces accross individuals, I want to compare all gaps with each other.

Run [run_combine_gaps_fasta_per_species.py](run_combine_gaps_fasta_per_species.py) to collect all gaps fasta files and combine them by species. During combining the individual ID/name is added to the sequence header.

Modified GD to start with analyzing the gaps fasta file: [genomedelta_main_gaps_analysis.sh](../GenomeDelta/genomedelta_main_gaps_analysis.sh)

This creates in GenomeDeltaResults a <species> folder.

Run it with combined reads to get clusters across all individuals:

```bash

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dfun/Dfun_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dfun --prefix Dfun  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dfun/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dimm/Dimm_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dimm --prefix Dimm  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dimm/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dsim/Dsim_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dsim --prefix Dsim  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dsim/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dbus/Dbus_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dbus --prefix Dbus  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dbus/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dhis/Dhis_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dhis --prefix Dhis  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dhis/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Drep/Drep_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Drep --prefix Drep  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Drep/genomedelta.nohup.out 2>&1 &

```

After GD is finished, we can run NCBI Blast and Dfam checks.

# Analyze clusters

```bash
python run_analyse_clusters.py GenomeDeltaResult/Dhis GenomeDeltaResult/Dbus GenomeDeltaResult/Dimm GenomeDeltaResult/Dfun GenomeDeltaResult/Dsim GenomeDeltaResult/Drep --min 5 --only-dfam-hits
```

Sample output:

```
Species	Cluster	Sequences	Unique_Individuals	DFAM_Hit	DFAM_Tandem
Dhis	cluster_17090.fasta	10	5	yes	no
Dhis	cluster_17183.fasta	5	5	no	yes
Dhis	cluster_21602.fasta	5	5	no	yes
Dhis	cluster_6857.fasta	5	5	no	no
Dhis	cluster_19652.fasta	30	5	no	no
Dhis	cluster_23764.fasta	7	5	yes	no
Dhis	cluster_30201.fasta	5	5	yes	yes
Dhis	cluster_14980.fasta	5	5	yes	no
Dhis	cluster_14995.fasta	5	5	yes	no
Dhis	cluster_21778.fasta	5	5	no	yes
Dhis	cluster_11349.fasta	8	5	no	no
Dhis	cluster_5841.fasta	5	5	no	yes
Dhis	cluster_27049.fasta	5	5	no	no
```

# TE Candidate Processing with DeviaTE

[run_gd_result_to_deviate.py](run_gd_result_to_deviate.py)

```bash
nohup python -u run_gd_result_to_deviate_candidates.py > deviate_candidates.log 2>&1 &
```

```bash
nohup python -u run_gd_result_to_deviate_all_gaps.py > deviate_all_gaps.log 2>&1 &
```

## Processing Workflow (Per Species)

For each species (e.g. `DSIM`, `DFUN`, ...), the script performs:

1. **Gather Input Files**

   * Scans `GenomeDeltaResult/*/file-GD-candidates.fasta`
   * Groups them by species (first 4 letters of folder name, uppercased)

2. **Combine Sequences**

   * All FASTA files for a species are merged into one
   * Each sequence header is prefixed with its folder name (indivisual):

     ```
     >Dsim01|original_header
     ```

3. **Deduplicate**

   * Runs `fastp` with deduplication enabled
   * Produces a single deduplicated `.fastq` file per species

4. **Run deviaTE**

   * Passes the deduplicated `.fastq` to `deviaTE`:

     ```bash
     deviaTE --input species_deduplicated_candidates.fastq
     ```

5. **Organize Output**

   * All results for a species are saved under:

     ```
     deviate/species/
     ├── species_combined_candidates.fasta
     ├── species_deduplicated_candidates.fastq
     └── deviaTE results
     ```

## Notes

* Skips steps if output files already exist
* Stops execution if any expected input file is missing
* Creates output folders automatically

## Example Output Structure

```
deviate/
├── Dsim/
│   ├── Dsim_combined_candidates.fasta
│   ├── Dsim_deduplicated_candidates.fastq
│   └── deviaTE output
├── Dfun/
│   └── ...
```

## Single copy gene normalization with DeviaTE

https://github.com/W-L/deviaTE/tree/master?tab=readme-ov-file#normalization-methods
https://github.com/W-L/deviaTE/tree/master?tab=readme-ov-file#special-use-case-drosophila

## Single copy gene with Busco

Not required for drosophila:

Here's a breakdown of the most relevant files/directories:

* **`busco_sequences/`**
  👉 This directory contains the **actual sequences** (in FASTA format) of the **complete**, **duplicated**, **fragmented**, and **missing** BUSCO genes found in your input genome or transcriptome.

* **`full_table.tsv`**
  👉 This is a detailed table listing **all BUSCOs in the database**, with status for each (Complete, Duplicated, Fragmented, or Missing), and the **sequence IDs** of matches in your reference genome. This is ideal if you want to extract the coordinates or further process the hits.

* **`missing_busco_list.tsv`**
  👉 Just a list of BUSCO IDs not found in your input.

* **`short_summary.txt`** and **`.json`**
  👉 These provide a summary of how many BUSCOs were found as single-copy, duplicated, etc.


## Get single copy gene for all speceis

```bash
python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dbus/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dfun/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Drep/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dhis/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dimm/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dsim/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna

```
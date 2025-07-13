# Working Directory

/mnt/data5/sarah/TE_Analysis

# Check candidates vs NCBI and DFAM Database/API

Use: 
* [run_ncbi_dfam_analysis.sh](run_ncbi_dfam_analysis.sh) to call
* [check_fasta_against_ncbi.py](check_fasta_against_ncbi.py)
* [check_fasta_against_dfam.py](check_fasta_against_dfam.py)


```bash
nohup bash run_ncbi_dfam_analysis.sh > ncbi_dfam_analysis.log 2>&1 &
```

## NCBI
This will perform a request against NCBI Blast using the python library. The results are stored in `NCBI_Analysis`

## DFAM
This will perform a request against DFAM API using HTTP (python request library). This is a multi step process:

1. send request with sequence to DFAM API. This returns a unique ID for the job.
2. request status with the unique ID. Once the job is done, this will return the result. Otherwise, a status message is returned.

DFAM mostly takes more time.

# Analysis of gaps by species accross individuals

To check matches accross individuals, I want to compare all gaps with each other.

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

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Brac/Brac_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Brac --prefix Brac  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Brac/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bgca/Bgca_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bgca --prefix Bgca  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bgca/genomedelta.nohup.out 2>&1 &

```

After GD is finished, we can run NCBI Blast and Dfam checks.

# Analyze clusters

Related to the combined analysis, i want to know if there are clusters with multiple individuals. Also i want to know which clusters have a DFAM hit and which clusters dont.

Run [run_analyse_clusters.py](run_analyse_clusters.py) to get a list of all clusters for one or many Result folders. It also compares it to DFAM hits and checks if the relevant cluster has a hit. If the DFAM hit is not found, it will say `error`. So its good to make sure the DFAM has finished before analysing the data.

Use optional parameters `--min 5 --only-dfam-hits`to filter for clusters with a minimum number of individuals and only show clusters with dfam hits.

If this is used on an individual GD folder, the number of individuals is always 1!

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

Columns: 

- Species: species ID
- Cluster: cluster file
- Sequences: # of sequences in the cluster
- Unique_Individuals: # number of unique individuals in the cluster 
- DFAM_Hit: cluster has a DFAM hit
- DFAM_Tandem: cluster has a DFAM tandem repeat

```bash
python run_analyse_clusters.py GenomeDeltaResult/Dhis GenomeDeltaResult/Dbus GenomeDeltaResult/Dimm GenomeDeltaResult/Dfun GenomeDeltaResult/Dsim GenomeDeltaResult/Drep --min 5 --only-dfam-hits
```

# TE Candidate Processing with DeviaTE

https://github.com/W-L/deviaTE 

Run [run_gd_result_to_deviate.py](run_gd_result_to_deviate.py) to map sequences against known TEs.

Initially 

```bash
nohup python -u run_gd_result_to_deviate_candidates.py > deviate_candidates.log 2>&1 &
```

```bash
nohup python -u run_gd_result_to_deviate_all_gaps.py > deviate_all_gaps.log 2>&1 &
```

## Processing Workflow of run_gd_result_to_deviate_candidates (Per Species)

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

## Processing Workflow of run_gd_result_to_deviate_all_gaps (Per Species)

For each species GD result folder, the script performs:

1. **Gather Input Files**

   * Scans `GenomeDeltaResult/*/file-GD.fasta` (= all gaps)

2. **Deduplicate**

   * Runs `fastp` with deduplication enabled

3. **Run deviaTE**

   * Passes the deduplicated `.fastq` to `deviaTE`:

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

For Dmel DeviaTE provides out of the box SCG. See here:

https://github.com/W-L/deviaTE/tree/master?tab=readme-ov-file#normalization-methods
https://github.com/W-L/deviaTE/tree/master?tab=readme-ov-file#special-use-case-drosophila

This is what i used to run it.

## Single copy gene with Busco

Not required for drosophila, here we use out of the box single copy gene isntead!

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

Run [get_single_copy_gene.py](get_single_copy_gene.py) to extract first single copy gene.

```bash
python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dbus/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dfun/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Drep/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dhis/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dimm/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dsim/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna

```

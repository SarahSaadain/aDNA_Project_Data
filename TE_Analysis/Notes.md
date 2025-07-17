# GenomeDelta

Below steps rely on having GenomeDelta outputs for the different individuals.

Example how to run GD normally:
* [Drosi GenomeDelta Notes](<../Drosi_NHM_2024/resources/run _genome_delta.md>)
* [Dsim GenomeDelta Notes](<../Dsim/resources/run _genome_delta.md>)

# Working Directory

All steps below are done from `/mnt/data5/sarah/TE_Analysis` if not stated otherwise.

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

## Step 1: Combine gaps from all individuals into a big fasta file

Run [run_combine_gaps_fasta_per_species.py](run_combine_gaps_fasta_per_species.py) to collect **all gaps** fasta files (`file-GD.fasta` from individuals) and combine them by species. During combining the individual ID/name is added to the sequence header.

## Step 2: Run the gaps with GD to get clusters

Modified GD to start with analyzing the gaps fasta file (no bam processing ... 
Makes consensus sequences of all found clusters across all individuals of the same species. Directly start with blast, then MSA, building clusters, consensus of cluster, ...): [genomedelta_main_gaps_analysis.sh](../GenomeDelta/genomedelta_main_gaps_analysis.sh)

This creates in GenomeDeltaResults a `<species>` (e.g. `Dsim`) folder (original runs create individuals as folders -> e.g. `Dsim01`).

Run it with combined reads (created by [run_combine_gaps_fasta_per_species.py](run_combine_gaps_fasta_per_species.py)) and provide them via new parameter `--gap_fasta` to get clusters across all individuals. It will do the same as GD, but in this case we provide the gaps fasta file ourselves. Since we give our own gap sequences, we dont need bam and ref genome files (used before to determine gaps).

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dfun/Dfun_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dfun --prefix Dfun  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dfun/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dimm/Dimm_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dimm --prefix Dimm  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dimm/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dsim/Dsim_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dsim --prefix Dsim  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dsim/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dbus/Dbus_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dbus --prefix Dbus  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dbus/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dhis/Dhis_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dhis --prefix Dhis  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dhis/genomedelta.nohup.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main_gaps_analysis.sh --gap_fasta /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Drep/Drep_combined_gaps.fasta --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Drep --prefix Drep  > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Drep/genomedelta.nohup.out 2>&1 &

```

## Step 3: Analyze with DFAM/NCBI Blast

After GD is finished, we can run NCBI Blast and Dfam checks.

## Step 4: Analyze with run_analyse_clusters.py

After Dfam we can use run_analyse_clusters.py to get a list of all clusters, their hits, regions and more. See below.

# Analyze clusters (run_analyse_clusters.py)

Related to the combined analysis, i want to know if there are clusters with multiple individuals. Also i want to know which clusters have a DFAM hit and which clusters dont.

Run [run_analyse_clusters.py](run_analyse_clusters.py) to get a list of all clusters for one or many Result folders. It also compares it to DFAM hits and checks if the relevant cluster has a hit. If the DFAM hit is not found, it will say `error`. So its good to make sure the DFAM has finished before analysing the data.

## run_analyse_clusters.py parameters

| **Parameter**       | **Type**   | **Description**                                                                        | **Default**           |
| ------------------- | ---------- | -------------------------------------------------------------------------------------- | --------------------- |
| `species_dirs`      | Positional | One or more input directories for species (e.g., `GenomeDelta/Dbus01`).                | Required              |
| `--min_individuals` | `int`      | Minimum number of unique individuals required per cluster to include it in the output. | `0`                   |
| `--min_ratio`       | `float`    | Minimum **sequence-to-individual** ratio required per cluster to include it.           | `0.00`                |
| `--dfam_dir`        | `str`      | Directory path to DFAM analysis results for checking repeat annotations.               | `'DFAM_Analysis'`     |
| `--only-dfam-hits`  | `flag`     | If set, **only** include clusters that have a DFAM match. Acts as a filter.            | `False` (not applied) |


If this is used on an individual GD folder, the number of individuals is always 1!

## Columns for run_analyse_clusters.py

| **Column Name**         | **Explanation**                                                                                                                                                                                                                                      |
| ----------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Species**             | Species code or name (e.g., `Dhis`, `Dbus`). Represents the organism from which the sequences come.                                                                                                                                                  |
| **Cluster**             | The filename of the FASTA file for a sequence cluster (e.g., `cluster_17090.fasta`). Each cluster likely groups homologous or related sequences.                                                                                                     |
| **Sequences**           | Total number of sequences in the FASTA cluster file. Count of lines starting with `>` in the FASTA.                                                                                                                                                  |
| **Unique\_Individuals** | Number of unique individuals (e.g., `Dhis01`, `Dbus04`) contributing sequences to the cluster. This is based on individual IDs extracted from sequence headers.                                                                                      |
| **Ratio\_Seq\_Indiv**   | Ratio of sequences to individuals: `Sequences / Unique_Individuals`. A value >1 means some individuals contribute multiple sequences.                                                                                                                |
| **DFAM\_Hit**           | Indicates whether the cluster matches a DFAM entry (likely a transposable element or repeat). `yes` or `no`.                                                                                                                                         |
| **DFAM\_Tandem**        | Indicates if the DFAM hit is a **tandem repeat** (repeats occurring adjacent to each other). `yes` or `no`.                                                                                                                                          |
| **Individuals**         | Comma-separated list of all unique individuals contributing to the cluster.                                                                                                                                                                          |
| **Regions**             | Comma-separated list of genomic regions (coordinates) corresponding to sequences in the cluster. Format: `scaffold:start-end`. If multiple sequences from an individual map to overlapping regions, these may be merged to a single coordinate span. **Note:** This region is the base for the GenomeDelta cluster but it is not identical. The consensus sequence from this cluster can have a different length than the region, but it will definately lie in the mentioned region. This is due to how GenomeDelta handles the creation of the consensus sequence. If one sequence of the cluster is bigger than the others, during alignment this can lead to `-` characters in the MSA file. If there are more `-` characters than NT during comparison of one Base, the Base will be omitted. If that happens on the beginning and/or on the end, the consensus sequence will be shorter than the source sequences. Du to this, it is not possible to give the exact region for the TE hit from Dfam. |
| **TE\_Descriptions** | Quoted description(s) of transposable element(s) hit (e.g., `"FW2_DM is a non-LTR retrotransposon"`).                                                     |
| **TE\_Coordinates**  | Coordinates in the consensus sequence of the TE hit, in `start-end` format. Adjusted for strand (reverse strand reports coordinates in descending order). |
| **TE\_Strand**       | Strand of the TE hit: `"+"` or `"-"`.                                                                                                                     |
| **TE\_Bit\_Score**   | Bit score of the DFAM match — reflects match quality.                                                                                                     |

## Run Drosi NHM data with run_analyse_clusters.py

```bash
python run_analyse_clusters.py GenomeDeltaResult/Dhis GenomeDeltaResult/Dbus GenomeDeltaResult/Dimm GenomeDeltaResult/Dfun GenomeDeltaResult/Dsim GenomeDeltaResult/Drep --min_individuals 5 --only-dfam-hits
```

### Example output for run_analyse_clusters.py

(Needs updating)

```
Species	Cluster	Sequences	Unique_Individuals	Ratio_Seq_Indiv	DFAM_Hit	DFAM_Tandem	Individuals	Regions	TE_Descriptions	TE_Coordinates	TE_Strand	TE_Bit_Score
Dhis	cluster_17090.fasta	10	5	2.0	yes	no	Dhis01,Dhis02,Dhis03,Dhis04,Dhis05	OY282582.1:29821006-29822677,OY282582.1:45901083-45904485,OY282582.1:45914821-45916819,OY282582.1:54477633-54479975,OY282582.1:58195753-58198013	"M4DM is a nonautonomous DNA transposon."	564-675	+	13.4
Dhis	cluster_17090.fasta	10	5	2.0	yes	no	Dhis01,Dhis02,Dhis03,Dhis04,Dhis05	OY282582.1:29821006-29822677,OY282582.1:45901083-45904485,OY282582.1:45914821-45916819,OY282582.1:54477633-54479975,OY282582.1:58195753-58198013	"M4DM is a nonautonomous DNA transposon."	1312-1818	+	146.4
Dhis	cluster_23764.fasta	7	5	1.4	yes	no	Dhis01,Dhis02,Dhis03,Dhis04,Dhis05	OY729166.1:28943535-28947857,OY729166.1:29588341-29593850	"DNAREP1_DM is a putative nonautonomous DNA transposable sequence."	1668-1769	-	57.6
Dhis	cluster_23764.fasta	7	5	1.4	yes	no	Dhis01,Dhis02,Dhis03,Dhis04,Dhis05	OY729166.1:28943535-28947857,OY729166.1:29588341-29593850	"BS retrotransposon."	1833-2111	-	37.8
Dhis	cluster_30201.fasta	5	5	1.0	yes	yes	Dhis01,Dhis02,Dhis03,Dhis04,Dhis05	OY729167.1:7714116-7716601	"DNAREP1_DM is a putative nonautonomous DNA transposable sequence."	1241-1363	-	43.2
Dhis	cluster_14980.fasta	5	5	1.0	yes	no	Dhis01,Dhis02,Dhis03,Dhis04,Dhis05	OY282582.1:33590097-33592374	"DNAREP1_DM is a putative nonautonomous DNA transposable sequence."	536-913	-	117.9
```

## Run Bger mapped against my ref genome with run_analyse_clusters.py

```bash
python run_analyse_clusters.py GenomeDeltaResult/Brac --min_individuals 1 --only-dfam-hits
```

## Run Bger mapped against GCA... with run_analyse_clusters.py

```bash
python run_analyse_clusters.py GenomeDeltaResult/Bgca --min_individuals 1 --only-dfam-hits
```
# TE Candidate Processing with DeviaTE

https://github.com/W-L/deviaTE 

## DeviaTE for candidates 

Run [run_gd_result_to_deviate_candidates.py](run_gd_result_to_deviate_candidates.py) to map sequences against known TEs. The candidates are combined by species

```bash
nohup python -u run_gd_result_to_deviate_candidates.py --scg_dir . --library_base TE_library_combined.fasta > deviate_candidates.log 2>&1 &
```

### Processing Workflow of run_gd_result_to_deviate_candidates.py (Per Species)

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

3. **Single Copy Gene**

   * Use `<species>_scg_<busco_id>_nt.fa` if availble. 
   * The SCG sequences are added to library provided by `--library_base`
   * The SCG sequences are added to fasta file for processing in DeviaTE

4. **Deduplicate**

   * Runs `fastp` with deduplication enabled
   * Produces a single deduplicated `.fastq` file per species

5. **Run deviaTE**

   * Passes the deduplicated `.fastq` to `deviaTE`:

## DeviaTE for individuals and all gaps 

Run [run_gd_result_to_deviate_individual_all_gaps.py](run_gd_result_to_deviate_individual_all_gaps.py) to map all gaps of individuals against known TEs. The candidates are combined by species

```bash
nohup python -u run_gd_result_to_deviate_individual_all_gaps.py > deviate_individual_all_gaps.log --scg_dir . --library_base TE_library_combined.fasta 2>&1 &
```

### Processing Workflow of run_gd_result_to_deviate_individual_all_gaps.py (Per Species)

For each species GD result folder, the script performs:

1. **Gather Input Files**

   * Scans `GenomeDeltaResult/*/file-GD.fasta` (= all gaps)

2. **Single Copy Gene**

   * Use `<species>_scg_<busco_id>_nt.fa` if availble. 
   * The SCG sequences are added to library provided by `--library_base`
   * The SCG sequences are added to fasta file for processing in DeviaTE

3. **Deduplicate**

   * Runs `fastp` with deduplication enabled

4. **Run deviaTE**

   * Passes the deduplicated `.fastq` to `deviaTE`:

## Single copy gene normalization with DeviaTE

For Dmel DeviaTE provides out of the box SCG. See here:

https://github.com/W-L/deviaTE/tree/master?tab=readme-ov-file#normalization-methods
https://github.com/W-L/deviaTE/tree/master?tab=readme-ov-file#special-use-case-drosophila

## Single copy gene with Busco

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

By adding `--num_genes n`it is possible to define up to n single copy genes that should be used. They are extracted into a file calles `<prefix>_scg.fasta`

```bash
python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dbus/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna --prefix Dbus

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dfun/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna --prefix Dfun

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Drep/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna --prefix Drep

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dhis/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna --prefix Dhis

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dimm/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna --prefix Dimm

python /mnt/data5/sarah/TE_Analysis/get_single_copy_gene.py /mnt/data5/sarah/TE_Analysis/BUSCO_Analysis_Dsim/BUSCO_Analysis/ /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna --prefix Dsim
```

Above commands create these files (`<species>_scg_<busco_id>_nt.fa`):

* Dbus_scg_11at7214_nt.fa
* Dfun_scg_11at7214_nt.fa 
* Dhis_scg_11at7214_nt.fa 
* Dimm_scg_11at7214_nt.fa  
* Drep_scg_11at7214_nt.fa
* Dsim_scg_11at7214_nt.fa

This files was manually creates based on the SCG in for Dmel `deviate_transposon_sequence_set_v10.2.fa` 
* Dmel_scg_deviate_nt.fa 

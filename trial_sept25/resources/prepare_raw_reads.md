# Download

# MD5
## Expected MD5


## MD5 Check



# Extract Reads

```bash
nohup tar -xvzf 235K3TLT3_5_R19504_20250910.tar.gz > unpack.log 2>&1 &
```

# Pre Checking

How to get a list of all relevant `fastq.gz` files. This also includes sub folders

```bash
find . -type f -name "*_R*.fastq.gz"
```

The original data also contains some undetermined reads. to exclude them add `! -iname "*undetermined*"`

```bash
find . -type f -name "*_R*.fastq.gz" ! -iname "*undetermined*"
```

By adding `| wc -l` we can get the number of files

```bash
find . -type f -name "*_R*.fastq.gz" ! -iname "*undetermined*" | wc -l
 ```

# FastQ File Renaming and Moving

## Step 1 - Collect fastq.gz

### Purpose
This scripts process sequencing data files organized in the following folder structure:

```
<sequencing\_run>/demultiplexed/<sequence\_id/>
```

Each `sequence_id` folder contains `.fastq.gz` files (paired-end reads: R1 and R2).

First we want to get all fastq.gz files into the same folder. Since each lane has the same file names for the sequences, we need to adapt the filenames so that they are unique. To achive this, we add the foldername of the tar.gz unzipped folder to the filename.

### Objective

- Move all `R1` and `R2` `.fastq.gz` files from their nested folders to the current working directory.
- Rename files to ensure unique filenames by appending the lane folder name
- Only move files matching `*_R1_*.fastq.gz` and `*_R2_*.fastq.gz` (exclude `I1` and `I2` index files).

#### Example

Original file path:

```
/232V3GLT3_1_R19096_20250626/demultiplexed/353205/353205_S28_R2_001.fastq.gz
```

Renamed and moved to current directory as:
```
./353205_S28_R2_001_232V3GLT3_1_R19096_20250626.fastq.gz
```

### How to Use

1. Run the **dry-run test script** to see what files would be moved and renamed, and count the total files:
2. The dry-run script prints:
   * Source file path
   * New filename it would move to
   * Total count of files matching the criteria
3. Once verified, run the other script to perform actual moves.


### Script

#### Dry Run
Run this script from the folder containing the unzipped tar.gz folders. 
This "dry run" script will not perform any changes, it will only show the user how the moving and renaming would be executed.

The script can simply be copied into the terminal as seen below.

```bash
#!/bin/bash

count=0

echo "Simulated file moves (only R1/R2, with full path context):"
echo "-----------------------------------------------------------"

for run_dir in ./*; do
    if [[ -d "$run_dir" ]]; then
        run_id=$(basename "$run_dir")
        demux_dir="$run_dir/demultiplexed"
        if [[ -d "$demux_dir" ]]; then
            for seq_dir in "$demux_dir"/3*; do
                seq_id=$(basename "$seq_dir")
                for fq in "$seq_dir"/*_R[12]_*.fastq.gz; do
                    if [[ -f "$fq" ]]; then
                        ((count++))
                        base_name=$(basename "$fq" .fastq.gz)
                        new_name="${base_name}.fastq.gz"
                        echo "$fq -> ./$(basename "$new_name")"
                    fi
                done
            done
        fi
    fi
done


echo "-----------------------------------------------------------"
echo "Total R1/R2 .fastq.gz files to be moved: $count"

```

### Real run

```bash
#!/bin/bash

count=0

echo "Simulated file moves (only R1/R2, with full path context):"
echo "-----------------------------------------------------------"

for run_dir in ./*; do
    if [[ -d "$run_dir" ]]; then
        run_id=$(basename "$run_dir")
        demux_dir="$run_dir/demultiplexed"
        if [[ -d "$demux_dir" ]]; then
            for seq_dir in "$demux_dir"/3*; do
                seq_id=$(basename "$seq_dir")
                for fq in "$seq_dir"/*_R[12]_*.fastq.gz; do
                    if [[ -f "$fq" ]]; then
                        ((count++))
                        base_name=$(basename "$fq" .fastq.gz)
                        new_name="${base_name}.fastq.gz"
                        #echo "$fq -> ./$(basename "$new_name")"
                        mv "$fq" "./$new_name"
                        echo "Moved: $fq -> ./$new_name"
                    fi
                done
            done
        fi
    fi
done


echo "-----------------------------------------------------------"
echo "Total R1/R2 .fastq.gz files to be moved: $count"

```

## Step 2 - Renaming the read files for the pipeline

### Objective
For the pipeline we want to perform additional renaming of the files:

1. We want to know the individual. This can be extracted from the `sequence id`. The relevant values are in the rename file in Google Docs.

## Renaming files

### Add indiviual and protocol

```bash
python /mnt/data5/sarah/aDNA/resources/rename.py /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/rename.csv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected --test
```

## Create folders and parent folders
```bash
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tDcar/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tDbus/raw/{ref_genome,reads} 
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tDfun/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tDimm/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tDrep/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tDmer/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tDwil/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tDana/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tBger/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tTni/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tOcap/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tLlit/raw/{ref_genome,reads}
mkdir -p /mnt/data2/sarah/aDNA_snakemake/tPvul/raw/{ref_genome,reads}
```

## Move to pipeline folder

The target folder is `<species>/raw/reads`.

```bash
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tDcar* /mnt/data2/sarah/aDNA_snakemake/tDcar/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tDbus* /mnt/data2/sarah/aDNA_snakemake/tDbus/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tDfun* /mnt/data2/sarah/aDNA_snakemake/tDfun/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tDimm* /mnt/data2/sarah/aDNA_snakemake/tDimm/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tDrep* /mnt/data2/sarah/aDNA_snakemake/tDrep/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tDmer* /mnt/data2/sarah/aDNA_snakemake/tDmer/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tDwil* /mnt/data2/sarah/aDNA_snakemake/tDwil/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tDana* /mnt/data2/sarah/aDNA_snakemake/tDana/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tBger* /mnt/data2/sarah/aDNA_snakemake/tBger/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tTni* /mnt/data2/sarah/aDNA_snakemake/tTni/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tOcap* /mnt/data2/sarah/aDNA_snakemake/tOcap/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tLlit* /mnt/data2/sarah/aDNA_snakemake/tLlit/raw/reads
mv /mnt/data2/sarah/sequencing_data/trial_sequencing_sept25/reads_collected/tPvul* /mnt/data2/sarah/aDNA_snakemake/tPvul/raw/reads
```

## Move ref genomes

```bash
cp /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna* /mnt/data2/sarah/aDNA_snakemake/tDbus/raw/ref_genome/
cp /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna* /mnt/data2/sarah/aDNA_snakemake/tDfun/raw/ref_genome/
cp /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna* /mnt/data2/sarah/aDNA_snakemake/tDimm/raw/ref_genome/
cp /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna* /mnt/data2/sarah/aDNA_snakemake/tDrep/raw/ref_genome/
cp /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/Bgerracon1.fasta* /mnt/data2/sarah/aDNA_snakemake/tBger/raw/ref_genome/
```

## Get missing ref genomes
### Dmer
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/961/210/405/GCA_961210405.1_dmerc_wildtype_assembly/GCA_961210405.1_dmerc_wildtype_assembly_genomic.fna.gz
```

### Dana
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/639/315/GCF_017639315.1_ASM1763931v2/GCF_017639315.1_ASM1763931v2_genomic.fna.gz
gunzip GCF_017639315.1_ASM1763931v2_genomic.fna.gz
```


### Llit
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/330/335/GCA_964330335.1_xgLitLitt3.hap1.1/GCA_964330335.1_xgLitLitt3.hap1.1_genomic.fna.gz
gunzip GCA_964330335.1_xgLitLitt3.hap1.1_genomic.fna.gz
```

### Dwil
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/018/902/025/GCF_018902025.1_UCI_dwil_1.1/GCF_018902025.1_UCI_dwil_1.1_genomic.fna.gz
gunzip GCF_018902025.1_UCI_dwil_1.1_genomic.fna.gz
```

### Tni
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/590/095/GCF_003590095.1_tn1/GCF_003590095.1_tn1_genomic.fna.gz
gunzip GCF_003590095.1_tn1_genomic.fna.gz
```

### Pvul
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/932/274/485/GCF_932274485.2_xgPatVulg1.2/GCF_932274485.2_xgPatVulg1.2_genomic.fna.gz
gunzip GCF_932274485.2_xgPatVulg1.2_genomic.fna.gz
```
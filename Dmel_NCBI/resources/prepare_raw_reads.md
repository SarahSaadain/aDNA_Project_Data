# Working directory on roco

/mnt/data2/sarah/Dmel_NCBI_Pipeline/

# Download

## SRR Numbers for prefetch

[Dmel_SRR_List.txt](<Dmel_SRR_List.txt >)

Note: Had to download latest version of SRA Toolit because the one from conad was not the latest. It had some issues with certificates. The latest (`/mnt/data2/sarah/Dmel_NCBI_Pipeline/sratoolkit.3.2.1-ubuntu64`) version worked. 

Run Download:

```bash
nohup sratoolkit.3.2.1-ubuntu64/bin/prefetch --option-file Dmel_SRR_List.txt > srr_prefetch.log 2>&1 &
```

it seems that SRR23876576 is corrupt:

2025-07-05T20:02:56 fasterq-dump.2.9.6 err: column undefined while opening cursor within virtual database module - failed to resolve column 'READ' idx '1'
2025-07-05T20:02:56 fasterq-dump.2.9.6 err: cmn_iter.c cmn_iter_range().VCursorOpen() -> RC(rcVDB,rcCursor,rcOpening,rcColumn,rcUndefined)

# Converting SRR to FASTQ

```bash
cd /mnt/data2/sarah/Dmel_NCBI_Pipeline/SRR_FILES

for sra_path in SRR*/SRR*.sra; do
  sra_file=$(basename "$sra_path" .sra)
  sra_dir=$(dirname "$sra_path")
  log_file="/mnt/data2/sarah/Dmel_NCBI_Pipeline/SRR_FILES/${sra_dir}/fasterq_${sra_file}.log"
  
  nohup fasterq-dump "$sra_path" \
    --outdir "$sra_dir" \
    --split-files \
    --threads 8 > "$log_file" 2>&1 &
done
```

What This Does: 
- Iterates over each .sra file.
- Runs fasterq-dump using:
    --outdir "$sra_dir" → outputs to the same directory as the .sra
    --split-files → outputs *_1.fastq and *_2.fastq if paired-end
    --threads 8 → adjusts threads (you can increase if your system has more CPUs)
- Logs output to: SRRxxxxxxx/fasterq_SRRxxxxxxx.log

Files will look like below:

```
./SRR23876580/SRR23876580.sra_1.fastq.gz
```

# FASTQ to FASTQ.GZ

Use below command to zip all files.

```bash
nohup find . -type f -name "*.fastq" -exec gzip {} \; > gzip_log.txt 2>&1 &
```

Command to remove all .fastq and *.sra files (after gzip is done):

```bash
find . -type f -name "*.fastq" -delete
find . -type f -name "*.sra" -delete
```

# Pre Checking

How to get a list of all relevant `fastq.gz` files. This also includes sub folders

```bash
find . -type f -name "*.fastq.gz"
```

By adding `| wc -l` we can get the number of files

```bash
find . -type f -name "*.fastq.gz" | wc -l
 ```

# File Renaming and Moving

## Step 1 - Collect fastq.gz

### Purpose
This scripts process sequencing data files organized in the following folder structure:

```
./<SRR FOLDER>/
```

First we want to get all fastq.gz files into the same folder. 

### Objective

- Move all `.fastq.gz` files from their nested folders to the current working directory.

### Move files

```bash
find . -type f -name "*.fastq.gz" -exec mv {} . \;
```

## Step 2 - Renaming the read files for the pipeline

### Objective
For the pipeline we want to perform additional renaming of the files:

1. We want to know the individual and protocol (not relevant here, but can be used for extra info in the pipeline). The relevant values are in the rename file in Google Docs.

This will get us from

```
SRR23876580.sra_1.fastq.gz
```
to 
```
Dmel01_NCBI_SRR23876586_R1.fastq.gz
```

## Renaming files

### Add indiviual and protocol

```bash
python /mnt/data5/sarah/aDNA/resources/rename.py rename_SRR.csv . --test
```

## Move to pipeline folder

The target folder is `<species>/raw/reads`.

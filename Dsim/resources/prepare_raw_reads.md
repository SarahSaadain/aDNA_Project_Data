# MD5
## Expected MD5

| File                                 | Expected MD5                       |
| ------------------------------------ | ---------------------------------- |
| `22TLJTLT4_1_R18987_20250612.tar.gz` | `20c32fbb97af170fb1a43ecbfcaa0f69` |
| `22TLJTLT4_2_R18987_20250612.tar.gz` | `8bf3e6e02ce710c5091d8b7f014ad540` |
| `22TLJTLT4_3_R18987_20250612.tar.gz` | `8f91beda9e7f305682de5e20660dc48d` |
| `22TLJTLT4_4_R18987_20250612.tar.gz` | `2485e39702d89ffb3dff7f702c78152e` |
| `22TLJTLT4_5_R18987_20250612.tar.gz` | `77ad7c2c7f5f888c2a329dea1f0f4281` |
| `22TLJTLT4_6_R18987_20250612.tar.gz` | `87cfbdc15a0f15d4c795ef206a235d4f` |
| `22TLJTLT4_7_R18987_20250612.tar.gz` | `8c4d437078a44e5b50430e976039835f` |
| `22TLJTLT4_8_R18987_20250612.tar.gz` | `d2833805841f8dba565305e378d5dd65` |

## MD5 Check

```bash
nohup sh -c 'md5sum 22TLJTLT4_*_R18987_20250612.tar.gz > md5sums.txt' &> md5sum.log &

nohup sh -c 'md5sum 22TLJTLT4_1_R18987_20250612.tar.gz > md5sum_1.txt' &> md5sum_1.log &
nohup sh -c 'md5sum 22TLJTLT4_2_R18987_20250612.tar.gz > md5sum_2.txt' &> md5sum_2.log &
nohup sh -c 'md5sum 22TLJTLT4_3_R18987_20250612.tar.gz > md5sum_3.txt' &> md5sum_3.log &
nohup sh -c 'md5sum 22TLJTLT4_4_R18987_20250612.tar.gz > md5sum_4.txt' &> md5sum_4.log &
nohup sh -c 'md5sum 22TLJTLT4_5_R18987_20250612.tar.gz > md5sum_5.txt' &> md5sum_5.log &
nohup sh -c 'md5sum 22TLJTLT4_6_R18987_20250612.tar.gz > md5sum_6.txt' &> md5sum_6.log &
nohup sh -c 'md5sum 22TLJTLT4_7_R18987_20250612.tar.gz > md5sum_7.txt' &> md5sum_7.log &
nohup sh -c 'md5sum 22TLJTLT4_8_R18987_20250612.tar.gz > md5sum_8.txt' &> md5sum_8.log &
```

Check if running:

ps aux | grep md5sum


# Extract Reads


```bash
nohup tar -xvzf 22TLJTLT4_1_R18987_20250612.tar.gz > extract_1.log 2>&1 &
nohup tar -xvzf 22TLJTLT4_2_R18987_20250612.tar.gz > extract_2.log 2>&1 &
nohup tar -xvzf 22TLJTLT4_3_R18987_20250612.tar.gz > extract_3.log 2>&1 &
nohup tar -xvzf 22TLJTLT4_4_R18987_20250612.tar.gz > extract_4.log 2>&1 &
nohup tar -xvzf 22TLJTLT4_5_R18987_20250612.tar.gz > extract_5.log 2>&1 &
nohup tar -xvzf 22TLJTLT4_6_R18987_20250612.tar.gz > extract_6.log 2>&1 &
nohup tar -xvzf 22TLJTLT4_7_R18987_20250612.tar.gz > extract_7.log 2>&1 &
nohup tar -xvzf 22TLJTLT4_8_R18987_20250612.tar.gz > extract_8.log 2>&1 &
```

Check if unzipping is still running:

```bash
ps aux | grep tar.gz
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
./22TLJTLT4_8_R18987_20250612/demultiplexed/347154/347154_S28_R1_001.fastq.gz
```

Renamed and moved to current directory as:
```
./347154_S28_R1_001_22TLJTLT4_8_R18987_20250612.fastq.gz
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
                        new_name="${base_name}_${run_id}.fastq.gz"
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
                        new_name="${base_name}_${run_id}.fastq.gz"
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

1. We want to know the individual and protocol (not relevant here, but can be used for extra info in the pipeline). This can be extracted from the `sequence id`. The relevant values are in the rename file in Google Docs.
2. We want to shorten the filname by
    * replacing the `lane folder` with `lane#` (example `22TLJTLT4_8_R18987_20250612` -> `lane8`)
    * removing `_001`as it is always the same.

This will get us from

```
347154_S28_R1_001_22TLJTLT4_8_R18987_20250612
```

to 
```
Dfun09_NHM_353205_S28_R2_lane1.fastq.gz
```

## Renaming files

### Add indiviual and protocol

```bash
python /mnt/data5/sarah/aDNA/resources/rename.py /mnt/data5/sarah/aDNA/Dsim/resources/rename_step1_runID_to_individual.csv /mnt/data5/sarah/aDNA/sequencing_data/Dsim/reads --test
```

### Folder to lane

```bash
python /mnt/data5/sarah/aDNA/resources/rename.py /mnt/data5/sarah/aDNA/Dsim/resources/rename_step2_folder_to_lane.csv /mnt/data5/sarah/aDNA/sequencing_data/Dsim/reads --test
```

### Remove 001 in name

```bash
python /mnt/data5/sarah/aDNA/resources/rename.py /mnt/data5/sarah/aDNA/Dsim/resources/rename_step3_remove_001.csv /mnt/data5/sarah/aDNA/sequencing_data/Dsim/reads --test
```

## Move to pipeline folder

The target folder is `<species>/raw/reads`.

```bash
mv /mnt/data5/sarah/aDNA/sequencing_data/Dsim/reads/Dsim* /mnt/data5/sarah/aDNA/Dsim/raw/reads
```
# Pre Checking

How to get a list of all relevant `fastq.gz` files. This also includes sub folders

```bash
find . -type f -name "*.fq.gz"
```

By adding `| wc -l` we can get the number of files

```bash
find . -type f -name "*.fq.gz" | wc -l
 ```

 864 read files in total

# FastQ File Renaming and Moving

## Step 1 - Collect fastq.gz files in Bger/raw/reads

### Script

#### Dry Run
Run this script from Bger/raw/reads. It will consider the sets (set1, set2) in the original folder.
This "dry run" script will not perform any changes, it will only show the user how the moving and renaming would be executed

```bash
#!/bin/bash

count=0

echo "Simulated file moves (only R1/R2, with full path context):"
echo "-----------------------------------------------------------"

for set_dir in /mnt/data4/sarah/BGI_aDNA2026/2026_01_JAN_F25A430002286-01_LIBquhzR_20260126140548_Philipp/Bger/*; do
    if [[ -d "$set_dir" ]]; then
        for fq in "$set_dir"/*.fq.gz; do
            if [[ -f "$fq" ]]; then
                ((count++))
                base_name=$(basename "$fq" .fq.gz)
                echo "$fq -> ./$(basename "$base_name").fastq.gz"
            fi
        done
    fi
done

echo "-----------------------------------------------------------"
echo "Total R1/R2 .fastq.gz files to be moved: $count"

```

### Real run

```bash
#!/bin/bash

count=0

echo "Moving fastq.gz files ..."

for set_dir in /mnt/data4/sarah/BGI_aDNA2026/2026_01_JAN_F25A430002286-01_LIBquhzR_20260126140548_Philipp/Bger/*; do
    if [[ -d "$set_dir" ]]; then
        for fq in "$set_dir"/*.fq.gz; do
            if [[ -f "$fq" ]]; then
                ((count++))
                base_name=$(basename "$fq" .fq.gz)
                echo "$fq -> ./$(basename "$base_name").fastq.gz"
                mv "$fq" "./$(basename "$base_name").fastq.gz"
            fi
        done
    fi
done

echo "-----------------------------------------------------------"
echo "Total files moved: $count"
```

## Renaming files

### Add indiviual and protocol

```bash
python rename.py rename_2026_step1.csv . --test
```

### Folder to lane

```bash
python resources/rename.py Bger/resources/rename_step2_folder_to_lane.csv Bger/raw/reads/ --test
```

### Remove 001 in name

```bash
python resources/rename.py Bger/resources/rename_step3_remove_001.csv Bger/raw/reads/ --test
```
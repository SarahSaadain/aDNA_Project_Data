-rwxr-xr-x 1 roco roco 190M Aug 25 20:04 file-1.fastq.gz
-rwxr-xr-x 1 roco roco 326M Aug 25 20:04 file-2.fastq.gz
-rwxr-xr-x 1 roco roco 256M Aug 25 20:04 file-3.fastq.gz
-rw-r--r-- 1 roco roco 1,2G Aug 25 20:05 file-4.fastq.gz
-rwxr-xr-x 1 roco roco 2,7G Aug 25 20:04 file-4_1.fastq.gz
-rwxr-xr-x 1 roco roco 3,1G Aug 25 20:05 file-4_2.fastq.gz
-rw-r--r-- 1 roco roco 580M Aug 25 20:05 file-5.fastq.gz
-rwxr-xr-x 1 roco roco 1,7G Aug 25 20:05 file-5_1.fastq.gz
-rwxr-xr-x 1 roco roco 2,0G Aug 25 20:05 file-5_2.fastq.gz
-rw-r--r-- 1 roco roco 241M Aug 25 20:05 file-6.fastq.gz
-rwxr-xr-x 1 roco roco 772M Aug 25 20:05 file-6_1.fastq.gz
-rwxr-xr-x 1 roco roco 947M Aug 25 20:05 file-6_2.fastq.gz
-rwxr-xr-x 1 roco roco   99 Aug 26 16:49 list


L3337.fastq.gz  1       340164
L3338.fastq.gz  2       340165
L3339.fastq.gz  3       340166
347745  4
347786  5
347787  6

## Renaming files

### Add indiviual

```bash
python /mnt/data5/sarah/aDNA/resources/rename.py /mnt/data2/sarah/aDNA_snakemake/Phum/resources/rename_step1_runID_to_individual.csv /mnt/data2/sarah/sequencing_data/Phum/reads --test
```

## Move to pipeline folder

The target folder is `<species>/raw/reads`.

```bash
mv /mnt/data2/sarah/aDNA/sequencing_data/Phum/reads/Phum* /mnt/data2/sarah/aDNA_snakemake/Phum/raw/reads
```
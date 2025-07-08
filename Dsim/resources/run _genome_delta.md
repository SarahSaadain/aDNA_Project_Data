# Running GenomeDelta

https://github.com/rpianezza/GenomeDelta

```bash
conda activate GenomeDelta
```

## species and related ref genomes

Overview:

GCF_016746395.2_Prin_Dsim_3.1_genomic.fna

## With nohup for combined files

Below is the list of commands using `nohup` to run them in the background, redirecting output to a log file at `/mnt/data5/sarah/GenomeDelta_run<species>/run_<species>.log`

Each of these commands will:

* Run in the background (`&`)
* Continue running even if the terminal is closed (`nohup`)
* Save all standard output and errors to a log file (`run_<species>.log`)

Files:

- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim01.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim02.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim03.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim04.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim05.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim06.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim07.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim08.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim09.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim10.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim11.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim12.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim13.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim17.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim18.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim19.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim20.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim_combined.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam

## With nohup for individual files

Here's the full set of `nohup` commands for **all** the individual BAM files (excluding combined), with logs written to `/mnt/data5/sarah/GenomeDelta_<individual>/nohup.log`.

### Dbus (Reference: GCF\_011750605.1\_ASM1175060v1\_genomic.fna)

Files:

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim01.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim01 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim02.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim02 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim03.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim03 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim04.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim04 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim05.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim05 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim05/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim06.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim06 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim06/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim07.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim07 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim07/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim08.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim08 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim08/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim09.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim09 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim09/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim10.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim10 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim10/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim11.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim11 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim11/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim12.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim12 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim12/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim13.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim13 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim13/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim17.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim17 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim17/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim18.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim18 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim18/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim19.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim19 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim19/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dsim/processed/GCF_016746395.2_Prin_Dsim_3.1_genomic/mapped/Dsim20.fastq_GCF_016746395.2_Prin_Dsim_3.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dsim/raw/ref_genome/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dsim20 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dsim20/nohup.log 2>&1 &
```

mkdir -p \
/mnt/data5/sarah/GenomeDeltaResult/Dsim01 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim02 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim03 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim04 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim05 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim06 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim07 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim08 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim09 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim10 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim11 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim12 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim13 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim17 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim18 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim19 \
/mnt/data5/sarah/GenomeDeltaResult/Dsim20 

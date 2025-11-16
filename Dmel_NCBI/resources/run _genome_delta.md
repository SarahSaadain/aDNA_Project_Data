# Running GenomeDelta

https://github.com/rpianezza/GenomeDelta

```bash
conda activate GenomeDelta
```

## With nohup for combined files

Below is the list of commands using `nohup` to run them in the background, redirecting output to a log file at `/mnt/data5/sarah/GenomeDelta_run<species>/run_<species>.log`

Each of these commands will:

* Run in the background (`&`)
* Continue running even if the terminal is closed (`nohup`)
* Save all standard output and errors to a log file (`run_<species>.log`)

Files:

```
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel01.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel02.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel03.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel04.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel05.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel06.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel07.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel08.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel09.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel10.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel11.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel12.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel13.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel14.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel15.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel16.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel17.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel18.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel20.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel21.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel22.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel23.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel24.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped/Dmel25.fastq_GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_sorted.bam
```

Ref genome:

```
Dmel_NCBI/raw/ref/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna
```

### **Dbus**

```bash
#!/bin/bash

# Set reference genome
REFERENCE="/mnt/data5/sarah/aDNA/Dmel_NCBI/raw/ref/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"

# Input BAM directory
BAM_DIR="/mnt/data5/sarah/aDNA/Dmel_NCBI/processed/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic/mapped"

# Output base directory
OUTPUT_BASE="/mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Dmel"

# Thread count
THREADS=15

# Loop over BAM files
for BAM_FILE in ${BAM_DIR}/Dmel*.bam; do
    # Extract sample name (e.g., Dmel01)
    BASENAME=$(basename "$BAM_FILE")
    SAMPLE_NAME=$(echo "$BASENAME" | cut -d. -f1)

    # Define output and log paths
    OUTPUT_DIR="${OUTPUT_BASE}/${SAMPLE_NAME}"
    LOG_FILE="${OUTPUT_DIR}/nohup_${SAMPLE_NAME}.log"

    # Make sure output dir exists
    mkdir -p "$OUTPUT_DIR"

    # Run GenomeDelta with nohup
    echo "Launching GenomeDelta for ${SAMPLE_NAME}..."
    nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
        --bam "$BAM_FILE" \
        --fa "$REFERENCE" \
        --of "$OUTPUT_DIR" \
        --t "$THREADS" > "$LOG_FILE" 2>&1 &
done
```
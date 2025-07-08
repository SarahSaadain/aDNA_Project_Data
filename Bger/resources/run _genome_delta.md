# Running GenomeDelta

https://github.com/rpianezza/GenomeDelta

```bash
conda activate GenomeDelta
```

# BAMs mapped to Bgerracon1

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Bger/processed/Bgerracon1/mapped/Bger1.fastq_Bgerracon1_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/Bgerracon1.fasta \
  --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger1_Bgerracon1 \
  --t 15 > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger1_Bgerracon1/run.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Bger/processed/Bgerracon1/mapped/Bger2.fastq_Bgerracon1_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/Bgerracon1.fasta \
  --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger2_Bgerracon1 \
  --t 15 > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger2_Bgerracon1/run.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Bger/processed/Bgerracon1/mapped/Bger3.fastq_Bgerracon1_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/Bgerracon1.fasta \
  --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger3_Bgerracon1 \
  --t 15 > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger3_Bgerracon1/run.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Bger/processed/Bgerracon1/mapped/Bger_combined.fastq_Bgerracon1_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/Bgerracon1.fasta \
  --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger_combined_Bgerracon1 \
  --t 15 > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bgerracon1/Bger_combined_Bgerracon1.out 2>&1 &

```

# BAMs mapped to GCA_000762945.2_Bger_2.0_genomic
```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Bger/processed/GCA_000762945.2_Bger_2.0_genomic/mapped/Bger1.fastq_GCA_000762945.2_Bger_2.0_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_000762945.2_Bger_2.0_genomic.fna \
  --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger1_GCA000762945 \
  --t 15 > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger1_GCA000762945/run.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Bger/processed/GCA_000762945.2_Bger_2.0_genomic/mapped/Bger2.fastq_GCA_000762945.2_Bger_2.0_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_000762945.2_Bger_2.0_genomic.fna \
  --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger2_GCA000762945 \
  --t 15 > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger2_GCA000762945/run.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Bger/processed/GCA_000762945.2_Bger_2.0_genomic/mapped/Bger3.fastq_GCA_000762945.2_Bger_2.0_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_000762945.2_Bger_2.0_genomic.fna \
  --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger3_GCA000762945 \
  --t 15 > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger3_GCA000762945/run.out 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Bger/processed/GCA_000762945.2_Bger_2.0_genomic/mapped/Bger_combined.fastq_GCA_000762945.2_Bger_2.0_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_000762945.2_Bger_2.0_genomic.fna \
  --of /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger_combined_GCA000762945 \
  --t 15 > /mnt/data5/sarah/TE_Analysis/GenomeDeltaResult/Bger_combined_GCA000762945/run_BgerRef.out 2>&1 &

```

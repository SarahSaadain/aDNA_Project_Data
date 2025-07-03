# Running GenomeDelta

https://github.com/rpianezza/GenomeDelta

```bash
conda activate GenomeDelta
```

## species and related ref genomes

Overview:

| Species | Reference Genome                            |
| ------- | ------------------------------------------- |
| Dbus    | GCF\_011750605.1\_ASM1175060v1\_genomic.fna |
| Dfun    | GCA\_018901825.1\_ASM1890182v1\_genomic.fna |
| Drep    | GCA\_018903745.1\_ASM1890374v1\_genomic.fna |
| Dhis    | GCA\_958299025.2\_idDroHist2.2\_genomic.fna |
| Dimm    | GCA\_963583835.1\_idDroImmi1.1\_genomic.fna |

with paths:

| Species  | Reference Genome Filename                   | Reference Genome Path                                                                | BAM File Path                                                                                                                                          | BAM Filename                                                          |
| -------- | ------------------------------------------- | ------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------ | --------------------------------------------------------------------- |
| **Dbus** | GCF\_011750605.1\_ASM1175060v1\_genomic.fna | `/mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna` | `/mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus_combined.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam` | `Dbus_combined.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam` |
| **Dfun** | GCA\_018901825.1\_ASM1890182v1\_genomic.fna | `/mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna` | `/mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun_combined.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam` | `Dfun_combined.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam` |
| **Drep** | GCA\_018903745.1\_ASM1890374v1\_genomic.fna | `/mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna` | `/mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep_combined.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam` | `Drep_combined.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam` |
| **Dhis** | GCA\_958299025.2\_idDroHist2.2\_genomic.fna | `/mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna` | `/mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis_combined.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam` | `Dhis_combined.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam` |
| **Dimm** | GCA\_963583835.1\_idDroImmi1.1\_genomic.fna | `/mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna` | `/mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm_combined.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam` | `Dimm_combined.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam` |



## Commands without nohup:

Command template:

```bash
bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/<species>/processed/<ref_genome_file_wo_extension>/mapped/<species>_combined.fastq_<ref_genome_file_wo_extension>_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/<ref_genome_file> \
  --of /mnt/data5/sarah/GenomeDelta_run<species> \
  --t 15
```

Here's a list of the full `bash` command for each species, with the correct substitutions for `<species>` and `<ref_genome_file>`:

### **Dbus**

```bash
bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus_combined.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDbus \
  --t 15
```

### **Dfun**

```bash
bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun_combined.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDfun \
  --t 15
```

### **Drep**

```bash
bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep_combined.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDrep \
  --t 15
```

### **Dhis**

```bash
bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis_combined.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDhis \
  --t 15
```

### **Dimm**

```bash
bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm_combined.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDimm \
  --t 15
```

## With nohup

Here’s the updated list of commands using `nohup` to run them in the background, redirecting output to a log file at `/mnt/data5/sarah/GenomeDelta_run<species>/run_<species>.log`

Each of these commands will:

* Run in the background (`&`)
* Continue running even if the terminal is closed (`nohup`)
* Save all standard output and errors to a log file (`run_<species>.log`)

Let me know if you want this wrapped into a single script to launch all five jobs at once.

### **Dbus**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus_combined.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDbus \
  --t 15 > /mnt/data5/sarah/GenomeDelta_runDbus/run_Dbus.log 2>&1 &
```

### **Dfun**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun_combined.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDfun \
  --t 15 > /mnt/data5/sarah/GenomeDelta_runDfun/run_Dfun.log 2>&1 &
```

### **Drep**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep_combined.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDrep \
  --t 15 > /mnt/data5/sarah/GenomeDelta_runDrep/run_Drep.log 2>&1 &
```

### **Dhis**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis_combined.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDhis \
  --t 15 > /mnt/data5/sarah/GenomeDelta_runDhis/run_Dhis.log 2>&1 &
```

### **Dimm**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm_combined.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Bger/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDimm \
  --t 15 > /mnt/data5/sarah/GenomeDelta_runDimm/run_Dimm.log 2>&1 &
```
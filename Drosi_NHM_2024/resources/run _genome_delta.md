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

## With nohup for combined files

Below is the list of commands using `nohup` to run them in the background, redirecting output to a log file at `/mnt/data5/sarah/GenomeDelta_run<species>/run_<species>.log`

Each of these commands will:

* Run in the background (`&`)
* Continue running even if the terminal is closed (`nohup`)
* Save all standard output and errors to a log file (`run_<species>.log`)

Files:

- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus_combined.fastq_GCF_011750605.- 1_/mnt/data5/sarah/aDNA/ASM1175060v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun_combined.fastq_GCA_018901825.- 1_/mnt/data5/sarah/aDNA/ASM1890182v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis_combined.fastq_GCA_958299025.2_idDroHist2.- 2_/mnt/data5/sarah/aDNA/genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm_combined.fastq_GCA_963583835.1_idDroImmi1.- 1_/mnt/data5/sarah/aDNA/genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep_combined.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam

### **Dbus**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus_combined.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDbus \
  --t 15 > /mnt/data5/sarah/GenomeDelta_runDbus/run_Dbus.log 2>&1 &
```

### **Dfun**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun_combined.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDfun \
  --t 15 > /mnt/data5/sarah/GenomeDelta_runDfun/run_Dfun.log 2>&1 &
```

### **Drep**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep_combined.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDrep \
  --t 15 > /mnt/data5/sarah/GenomeDelta_runDrep/run_Drep.log 2>&1 &
```

### **Dhis**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis_combined.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDhis \
  --t 15 > /mnt/data5/sarah/GenomeDelta_runDhis/run_Dhis.log 2>&1 &
```

### **Dimm**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm_combined.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_runDimm \
  --t 15 > /mnt/data5/sarah/GenomeDelta_runDimm/run_Dimm.log 2>&1 &
```

## With nohup for individual files

Here's the full set of `nohup` commands for **all** the individual BAM files (excluding combined), with logs written to `/mnt/data5/sarah/GenomeDelta_<individual>/nohup.log`.

### Dbus (Reference: GCF\_011750605.1\_ASM1175060v1\_genomic.fna)

Files:


- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus01.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus02.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus03.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus04.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus05.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus06.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus07.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus08.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus09.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus10.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus01.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dbus01 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dbus01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus02.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dbus02 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dbus02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus03.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dbus03 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dbus03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus04.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dbus04 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dbus04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus05.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dbus05 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dbus05/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus06.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dbus06 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dbus06/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus07.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dbus07 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dbus07/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus08.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dbus08 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dbus08/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus09.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dbus09 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dbus09/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus10.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dbus10 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dbus10/nohup.log 2>&1 &
```

### Dfun (Reference: GCA\_018901825.1\_ASM1890182v1\_genomic.fna)

Files:

- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun01.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun02.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun03.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun04.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun05.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun06.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun07.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun08.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun09.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun10.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun01.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dfun01 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dfun01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun02.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dfun02 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dfun02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun03.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dfun03 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dfun03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun04.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dfun04 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dfun04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun05.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dfun05 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dfun05/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun06.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dfun06 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dfun06/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun07.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dfun07 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dfun07/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun08.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dfun08 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dfun08/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun09.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dfun09 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dfun09/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun10.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dfun10 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dfun10/nohup.log 2>&1 &
```

### Dhis

Files:

- /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis01.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis02.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis03.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis04.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis05.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis01.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dhis01 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dhis01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis02.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dhis02 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dhis02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis03.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dhis03 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dhis03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis04.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dhis04 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dhis04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis05.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dhis05 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dhis05/nohup.log 2>&1 &
```

### Dimm

Files:

- /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm01.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm02.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm03.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm04.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm05.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm01.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dimm01 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dimm01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm02.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dimm02 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dimm02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm03.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dimm03 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dimm03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm04.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dimm04 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dimm04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm05.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Dimm05 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Dimm05/nohup.log 2>&1 &
```

### Drep (Reference: GCA\_018903745.1\_ASM1890374v1\_genomic.fna)

Files:

- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep01.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep02.fastq_GCA_018903745.1_ASM1890374v1_genomic.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep02.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep03.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep04.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep05.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep06.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep07.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep08.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep09.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam
- /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep10.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep01.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Drep01 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Drep01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep02.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Drep02 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Drep02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep03.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Drep03 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Drep03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep04.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Drep04 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Drep04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep05.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Drep05 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Drep05/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep06.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Drep06 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Drep06/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep07.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Drep07 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Drep07/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep08.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Drep08 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Drep08/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep09.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Drep09 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Drep09/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep10.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDelta_Drep10 \
  --t 15 > /mnt/data5/sarah/GenomeDelta_Drep10/nohup.log 2>&1 &
```
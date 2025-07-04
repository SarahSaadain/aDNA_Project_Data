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
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/DbusCombined/nohup_Dbus.log 2>&1 &
```

### **Dfun**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun_combined.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/DfunCombined/nohup_Dfun.log 2>&1 &
```

### **Drep**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep_combined.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/DrepCombined/nohup_Drep.log 2>&1 &
```

### **Dhis**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis_combined.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dhis \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/DhisCombined/nohup_Dhis.log 2>&1 &
```

### **Dimm**

```bash
nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm_combined.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dimm \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/DimmCombined/nohup_Dimm.log 2>&1 &
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
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus01 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dbus01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus02.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus02 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dbus02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus03.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus03 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dbus03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus04.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus04 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dbus04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus05.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus05 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dbus05/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus06.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus06 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dbus06/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus07.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus07 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dbus07/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus08.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus08 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dbus08/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus09.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus09 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dbus09/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dbus/processed/GCF_011750605.1_ASM1175060v1_genomic/mapped/Dbus10.fastq_GCF_011750605.1_ASM1175060v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dbus/raw/ref_genome/GCF_011750605.1_ASM1175060v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dbus10 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dbus10/nohup.log 2>&1 &
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
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun01 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dfun01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun02.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun02 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dfun02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun03.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun03 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dfun03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun04.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun04 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dfun04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun05.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun05 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dfun05/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun06.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun06 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dfun06/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun07.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun07 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dfun07/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun08.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun08 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dfun08/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun09.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun09 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dfun09/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dfun/processed/GCA_018901825.1_ASM1890182v1_genomic/mapped/Dfun10.fastq_GCA_018901825.1_ASM1890182v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dfun/raw/ref_genome/GCA_018901825.1_ASM1890182v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dfun10 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dfun10/nohup.log 2>&1 &
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
  --of /mnt/data5/sarah/GenomeDeltaResult/Dhis01 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dhis01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis02.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dhis02 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dhis02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis03.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dhis03 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dhis03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis04.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dhis04 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dhis04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dhis/processed/GCA_958299025.2_idDroHist2.2_genomic/mapped/Dhis05.fastq_GCA_958299025.2_idDroHist2.2_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dhis/raw/ref_genome/GCA_958299025.2_idDroHist2.2_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dhis05 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dhis05/nohup.log 2>&1 &
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
  --of /mnt/data5/sarah/GenomeDeltaResult/Dimm01 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dimm01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm02.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dimm02 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dimm02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm03.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dimm03 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dimm03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm04.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dimm04 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dimm04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Dimm/processed/GCA_963583835.1_idDroImmi1.1_genomic/mapped/Dimm05.fastq_GCA_963583835.1_idDroImmi1.1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Dimm/raw/ref_genome/GCA_963583835.1_idDroImmi1.1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Dimm05 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Dimm05/nohup.log 2>&1 &
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
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep01 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Drep01/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep02.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep02 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Drep02/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep03.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep03 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Drep03/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep04.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep04 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Drep04/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep05.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep05 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Drep05/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep06.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep06 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Drep06/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep07.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep07 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Drep07/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep08.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep08 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Drep08/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep09.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep09 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Drep09/nohup.log 2>&1 &

nohup bash /mnt/data5/sarah/GenomeDelta/linux/main.sh \
  --bam /mnt/data5/sarah/aDNA/Drep/processed/GCA_018903745.1_ASM1890374v1_genomic/mapped/Drep10.fastq_GCA_018903745.1_ASM1890374v1_genomic_sorted.bam \
  --fa /mnt/data5/sarah/aDNA/Drep/raw/ref_genome/GCA_018903745.1_ASM1890374v1_genomic.fna \
  --of /mnt/data5/sarah/GenomeDeltaResult/Drep10 \
  --t 10 > /mnt/data5/sarah/GenomeDeltaResult/Drep10/nohup.log 2>&1 &
```

```bash
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dbus01
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dbus02
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dbus03
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dbus04
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dbus05
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dbus06
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dbus07
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dbus08
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dbus09
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dbus10
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dfun01
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dfun02
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dfun03
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dfun04
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dfun05
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dfun06
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dfun07
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dfun08
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dfun09
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dfun10
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dhis01
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dhis02
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dhis03
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dhis04
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dhis05
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dimm01
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dimm02
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dimm03
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dimm04
mkdir /mnt/data5/sarah/GenomeDeltaResult/Dimm05
mkdir /mnt/data5/sarah/GenomeDeltaResult/Drep01
mkdir /mnt/data5/sarah/GenomeDeltaResult/Drep02
mkdir /mnt/data5/sarah/GenomeDeltaResult/Drep03
mkdir /mnt/data5/sarah/GenomeDeltaResult/Drep04
mkdir /mnt/data5/sarah/GenomeDeltaResult/Drep05
mkdir /mnt/data5/sarah/GenomeDeltaResult/Drep06
mkdir /mnt/data5/sarah/GenomeDeltaResult/Drep07
mkdir /mnt/data5/sarah/GenomeDeltaResult/Drep08
mkdir /mnt/data5/sarah/GenomeDeltaResult/Drep09
mkdir /mnt/data5/sarah/GenomeDeltaResult/Drep10
```
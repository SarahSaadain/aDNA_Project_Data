On roco:

```bash
cd /mnt/data5/sarah
```

```bash
git clone https://github.com/SarahSaadain/aDNA_Pipeline_Snakemake.git demo_pipeline
```

```bash
mv demo_pipeline_dataset demo_pipeline/demo
```

```bash
cd demo_pipeline
```

```bash
conda activate snakemake
```

```bash
snakemake --cores 8 --use-conda --keep-going -n
```

```bash
snakemake --cores 20 --use-conda --keep-going
```

Download data to laptop:

```bash
bash ~/Documents/aDNA/fetch_results_demo.sh
```
# Working Directory

/mnt/data5/sarah/TE_Analysis

# Check candidates vs NCBI and DFAM Database/API

Use: 
* [run_ncbi_dfam_analysis.sh](run_ncbi_dfam_analysis.sh) to call
* [heck_fasta_against_ncbi.py](check_fasta_against_ncbi.py)
* [check_fasta_against_dfam.py](check_fasta_against_dfam.py)

```bash
nohup bash run_ncbi_dfam_analysis.sh > ncbi_dfam_analysis.log 2>&1 &
```

# TE Candidate Processing with DeviaTE

[un_gd_result_to_deviate.py](run_gd_result_to_deviate.py)

## Processing Workflow (Per Species)

For each species (e.g. `DSIM`, `DFUN`, ...), the script performs:

1. **Gather Input Files**

   * Scans `GenomeDeltaResult/*/file-GD-candidates.fasta`
   * Groups them by species (first 4 letters of folder name, uppercased)

2. **Combine Sequences**

   * All FASTA files for a species are merged into one
   * Each sequence header is prefixed with its folder name (indivisual):

     ```
     >Dsim01|original_header
     ```

3. **Deduplicate**

   * Runs `fastp` with deduplication enabled
   * Produces a single deduplicated `.fastq` file per species

4. **Run deviaTE**

   * Passes the deduplicated `.fastq` to `deviaTE`:

     ```bash
     deviaTE --input species_deduplicated_candidates.fastq
     ```

5. **Organize Output**

   * All results for a species are saved under:

     ```
     deviate/species/
     ├── species_combined_candidates.fasta
     ├── species_deduplicated_candidates.fastq
     └── deviaTE results
     ```

## Notes

* Skips steps if output files already exist
* Stops execution if any expected input file is missing
* Creates output folders automatically

## Example Output Structure

```
deviate/
├── Dsim/
│   ├── Dsim_combined_candidates.fasta
│   ├── Dsim_deduplicated_candidates.fastq
│   └── deviaTE output
├── Dfun/
│   └── ...
```

## Tools Used

| Tool      | Role                   |
| --------- | ---------------------- |
| `fastp`   | Sequence deduplication |
| `deviaTE` | TE insertion analysis  |


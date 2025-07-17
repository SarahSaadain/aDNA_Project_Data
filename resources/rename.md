
# Usage  

The script reads a CSV file containing old and new filename mappings and renames the files accordingly in the specified folder. 

**How It Works**
- The script reads the CSV file and stores the old-to-new name mappings.
- It scans the specified folder for filenames that contain any of the old names.
- If a match is found, the filename is updated accordingly.
- If --test is enabled, it prints the changes without renaming the files.

This ensures efficient and structured renaming of raw read files within the pipeline.

# CSV Format  

The CSV file should have two columns:  
- **Column 1:** The original filename or a substring to be replaced.  
- **Column 2:** The new filename or replacement substring.  

**Example (`rename_list.csv`):**  

```
344209,Dsim19_2trial_344209
344210,Dsim19_2trial_344210
```

# Running the Script  

Navigate to the project root and execute the script:  

```bash
python resources/rename.py rename_list.csv path/to/raw_reads/
```

This will rename the files in `path/to/raw_reads/` based on the mappings in `rename_list.csv`.

## Test Mode
To preview the changes without renaming files, use the `--test` flag:

```bash
python resources/rename.py rename_list.csv path/to/raw_reads/ --test
```

This will print the planned renaming actions without modifying any files.
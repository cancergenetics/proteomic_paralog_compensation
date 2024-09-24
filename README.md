# Proteomic paralog compensation analysis

This repository contains code to analyze proteomic data and associate gene losses with changes in the abundance of paralogous proteins.

## Input data

Data sources are listed in data_sources.md.

## Analysis pipeline

### 1. HAP1 proteomic analysis

- **Script**: `./src/run_HAP1_paralog_tests.py`
- **Description**: Processes proteomic data from 34 knockouts of paralogous proteins in HAP1 cell lines, associating KO with changes in paralog abundance.
- **Output**:
   - Processed HAP1 proteomic dataset (data in Supplementary table 1) 
  - All self-tests (associations of KOs with changes in the abundance of the knocked-out protein) (Data in Supplementary table 2)
  - All paralog tests (associations of KOs with changes in paralog abundance) (Data in Supplementary table 3)

### 2. Protein residuals generation

- **Script**: `./src/create_residuals.py`
- **Description**: Generates protein residuals (residuals of linear models predicting study/lineage-corrected protein abundance with study/lineage-corrected mRNA abundance).
- **Output**: Prot residual dataset and additional/intermediate residual datasets used for plotting.
- **Can speed up with parallel processing if setting nb_workers.**

### 3. CPTAC paralog compensation analysis

- **Script**: `./src/main.py`
- **Description**: Analyzes the CPTAC proteomic and transcriptomic datasets, associating hemizygous loss with changes in paralog abundance.
- **Arguments**:
  - `--runwith`: Specifies whether to run with proteomic, transcriptomic, or protein residual data.
  - `--nb_workers`: Set to 1 by default. Can be increased for parallel processing using pandarallel.
  - `--HAP1_overlap`: When set, only tests pairs testable with the HAP1 dataset.
- **Output**:
    - All CPTAC proteomic self-tests
    - All CPTAC proteomic paralog tests 
    - All CPTAC transcriptomic and residual paralog tests 
    - All CPTAC Fishers Exact Tests for biological enrichments 
    - All CPTAC annotations of biological enrichments
    - All t-test results for quantitative enrichments 
- **Default Output Directory**: `./output/output_CPTAC/<runwith>`
- **Notes**: 
  - Transcriptomic and residual analyses are only run on pairs testable with the proteomic dataset.
  - Proteomic results must be in their respective directory for these options to work.
  - HAP1 output needs to be present for the `--HAP1_overlap` option to work.

### 4. Visualization

- **Notebook**: `Visualize.ipynb` 
- **Description**: Generates all separate panels of the main and supplementary figures.
- **Additional Calculations**: 
  - P-values for ubiquitination enrichment
  - Enrichment of transcriptional hits in proteomic hits
- **Default Output Directory**: `./output/figures`

### 5. Generate supplementary tables 

- **Script**: `clean_and_format_supptables.py` 
- **Description**: Processes results CSVs output by scripts to generate supplementary tables of the paper.
- **Default Output Directory**: `./output/supp_tables`




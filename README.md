# ILS Ceramide Ring Trial Results

## Prerequisites

The evaluation code requires R 4.3 or later. Additionally, some packages may require system library dependencies. Please consult the package installation output for details for your environment.

### Ubuntu

```
sudo apt update && sudo apt install cmake libcurl4-openssl-dev build-essential libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

## Starting the Evaluation

The data for the evaluation is located in the `data` folder. Scripts used during the evaluation are located in the `R` folder, the `_targets.R` and `evaluation.R` files.

We use `renv` to manage the dependencies of the code. To install all packages, as defined in the `renv.lock` file, run:

```
renv::restore() 
```

and follow the instructions. Hit `CTRL+SHIFT+F10` to restart the R session.

In order to run the evaluation, execute the following script, for example within RStudio:

```
evaluation.R
```

This will run the code steps contained in the `_targets.R` file contained within the same folder.
Please see the respective files for further information.

## Project Structure

### Protocols and Templates Folder

The folder `protocol_and_templates` contains a reporting template for the recommended ring trial protocol in the `Report Template Standard Protocol Harmonization Ceramides.xlsx`, 
the file `Report Template Preferred Protocol Harmonization Ceramides.xlsx` that can be used as a reporting template for your own, lab-specific methods, and the file `Standard Protocol Harmonization Ceramides March2020.docx` containing the detailed standard analytical protocol for the ring trial. 

The EXCEL files contain the following sheets:

1. Sheet *Participant Information*
2. Sheet *Summary*
3. Sheet *Extraction - Overview*
4. Sheet *Extraction - Procedure*
5. Sheet *LC-MRM - Overview*
6. Sheet *LC-MRM - Sample Sequence*
7. Sheet *LC-MRM - Intra Assay QC Proc*
8. Sheet *LC-MRM - Intra Assay QC Res*
9. Sheet *LC-MRM - Raw Peak Areas*

Sheet no. 8, *LC-MRM - Intra Assay QC Res* is the one relevant for the concentration calculation and intra lab evaluation.

### Data Folder

The folder `data` contains the original report tables provided by each participating lab in the `original-reports` folder and the final ones used for the comparison in the `reports` folder. 
Some original reports needed to be edited to conform to the format required by the evaluation pipeline. An overview of dataset completeness, remarks by the submitters, as well as by curation is available in the `dataset-curation-summary.xlsx` file. 

The folder `definitions` contains information related to blank types and ranges in the different lab reports, theoretical ratios between labeled und unlabeled standards, labeled and unlabeled standard concentrations, expected calibration line sample concentrations, as well as a mapping table of MS instruments to the respective lab reports.

The folder `ring-trial-comparison` contains tabular results for the other ring trials that we compared our ceramide concentrations against.

### R folder

The `R` folder contains reusable portions of code used within the evaluation workflow, e.g. for plotting, calibration curve-related calculations etc.

## Contributors

- Nils Hoffmann - data integration workflow and plots
- Bo J. Burla - ring trial comparison and plots


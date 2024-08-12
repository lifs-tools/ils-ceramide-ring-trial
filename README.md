# ILS Ceramide Ring Trial Results
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10081970.svg)](https://doi.org/10.5281/zenodo.10081970)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12632988.svg)](https://zenodo.org/doi/10.5281/zenodo.12632988)

## Prerequisites

The evaluation code requires R 4.3 or later. Additionally, some packages may require system library dependencies. Please consult the package installation output for details for your environment.

This git project uses lfs. In order to clone the repository, you need to have git-lfs installed. Please see the git-lfs [installation instructions](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage). The file `.gitattributes` contains the file extensions that are tracked by git-lfs. To add additional file extensions, please use `git lfs track "*.ext"`, replacing `ext` with the file extension you want to track.

Clone the repository with:

```
git clone https://github.com/lifs-tools/ils-ceramide-ring-trial.git
```

or using ssh:

```
git clone git@github.com:lifs-tools/ils-ceramide-ring-trial.git
```

### Ubuntu

In order to run the evaluation in Ubuntu, you will need to install system packages before installing the R packages:

```
sudo apt update && sudo apt install cmake libcurl4-openssl-dev build-essential libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

## Opening the Project

If you use RStudio, the provided .Rprofile file should handle automatic activation of the renv environment and installation of missing dependencies. 
Should this not work, please install the `renv` package (either via RStudio's Tools->Install Packages or from an R console)

```
install.packages("renv")
```

Then run 

```
renv::restore()
```

to install any missing packages from the `renv.lock` file.

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

Please note that raw data provided by the labs is not included in this repository. The data was provided by the labs to the study coordination at the [SLING](https://sling.sg/) / [SingMass](https://singmass.sg/) at NUS Singapore. To maintain the privacy of the labs, the raw data was converted to mzML and uploaded to [Zenodo](https://zenodo.org/doi/10.5281/zenodo.12632988) for the labs that have provided their raw data. Please note that the lab IDs used within this repository differ from the ones used in the final manuscript. The [mapping table](https://zenodo.org/doi/10.5281/zenodo.12632988) maps the ids used in this repository (`LabId`) to the ids used in the final manuscript (`LabNum`).

### R folder

The `R` folder contains reusable portions of code used within the evaluation workflow, e.g. for plotting, calibration curve-related calculations etc.

The file `R/definitions.R` contains configurable options for data input, calibration curve calculation, NA value definitions, ceramide names and plotting defaults. 

### Manuscript folder

The `manuscript` folder contains additional data used by the `manuscript-figures-tables.Rmd` Rmarkdown notebook that was used as the last step in the targets pipeline to generate figures, tables and supplementary files for use in the final manuscript. The notebook uses the data created by the previous pipeline steps and combines them with detailed data from the other ring trials for comparison.

### Supplementary Data

To simplify review and comparison of the raw and quantified values for all labs, we provide the file `Suppl table concentration values.xlsx`.
This spreadsheet contains three result sheets: `Results table 1` with the full dataset of all labs, `Results table 2`, which reports mean and SD of injection replicates and `Results table 3`, which reports mean and SD of extraction replicates (injection replicates averaged). Please see the `Legend` sheet for an explanation of the column headers and remarks concerning missing values.

## Contributors

- Nils Hoffmann - data integration workflow and plots
- Bo J. Burla - ring trial comparison and plots
- Federico T. Torta - participating lab data, supplementary tables

# Help & Support

Please [file an issue](https://github.com/lifs-tools/ils-ceramide-ring-trial/issues/new/choose) within this repository if you have questions about the workflow or have identified any bugs.

# Acknowledgements

We would like to thank all members of the International Lipidomics Society's [interest group for Reference materials and biological reference ranges](https://lipidomicssociety.org/interest_groups/reference-materials-and-biological-reference-ranges/) for their support and valuable input and all participants of the ring trial for their feedback and patience.

# Inquiring about your lab id

Please note that the data analysis team had no access to lab identities due to the blinding of the analysis reports before transmission to us. We can therefor not provide you your lab id to re-identify your results. Please contact the study coordinator Federico Torta at the [SLING](https://sling.sg/) / [SingMass](https://singmass.sg/) at NUS Singapore for details. Also, please see the [mapping table](https://zenodo.org/doi/10.5281/zenodo.12632988), which maps the ids used in this repository (`LabId`) to the ids used in the final manuscript (`LabNum`). 

# References

Main publication

[1] Torta F, Hoffmann N, Burla B *et al.*, Concordant inter-laboratory derived concentrations of ceramides in human plasma reference materials via authentic standards. 2024

[1a] Evaluation code and data for the ceramide ring trial. https://doi.org/10.5281/zenodo10081970

[1b] mzML datasets for the ceramide ring trial. https://doi.org/10.5281/zenodo.12632988

---

[2] A. Gustavo González, M. Ángeles Herrador. A practical guide to analytical method validation, including measurement uncertainty and accuracy profiles. TrAC Trends in Analytical Chemistry 2007;26(3). doi: [10.1016/j.trac.2007.01.009](https://doi.org/10.1016/j.trac.2007.01.009)

[3] FDA, Bioanalytical Method Validation - Guidance for Industry. 2018 May. https://www.fda.gov/files/drugs/published/Bioanalytical-Method-Validation-Guidance-for-Industry.pdf 

[4] DA Armbruster, T Pry. Limit of blank, limit of detection and limit of quantitation. Clin Biochem Rev. 2008 Aug;29 Suppl 1(Suppl 1):S49-52. https://pubmed.ncbi.nlm.nih.gov/18852857

[5] V. Barwick (Ed), Eurachem/CITAC Guide: Guide to Quality in Analytical Chemistry: An Aid to Accreditation (3rd ed. 2016). ISBN 978-0-948926-32-7. Available from https://www.eurachem.org.

[6] B. Magnusson and U. Örnemark (eds.) Eurachem Guide: The Fitness for Purpose of Analytical Methods – A Laboratory Guide to Method Validation and Related Topics, (2nd ed. 2014). ISBN 978-91-87461-59-0. Available from www.eurachem.org.

[7] D. Kauhanen, M. Sysi-Ah, KM. Koistinen, R. Laaksonen, J. Sinisalo, K. Ekroos. Development and validation of a high-throughput LC-MS/MS assay for routine measurement of molecular ceramides. Anal Bioanal Chem. 2016 May;408(13):3475-83. doi: [10.1007/s00216-016-9425-z](https://doi.org/10.1007/s00216-016-9425-z)

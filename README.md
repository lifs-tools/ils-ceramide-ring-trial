# ILS Ceramide Ring Trial Results

## Starting the evaluation

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

## Contributors

- Nils Hoffmann - data integration workflow and plots
- Bo J. Burla - ring trial comparison and plots

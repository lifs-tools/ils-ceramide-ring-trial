library(targets)
library(future)
library(future.callr)
plan(callr)
options(tidyverse.quiet = TRUE)

tar_option_set(
  packages = c(
    "tidyverse",
    "knitr",
    "ggpubr",
    "broom",
    "scales",
    "ggbeeswarm",
    "ggthemes",
    "ggsci",
    "progress",
    "here",
    "patchwork",
    "glue",
    "ggh4x",
    "tarchetypes"
  )
)

# print debug output for specific function
#tar_option_set(debug = "readAssayTables")  

source("R/definitions.R")
source("R/functions_io.R")
source("R/functions_filterData.R")
source("R/functions_steps.R")
source("R/functions_methods.R")
source("R/functions_ggpubr-custom.R")
source("R/functions_plots.R")

list(
  ################################################################################
  # Load and combine reports
  ################################################################################  
  tar_target(standardReports, standardReportsToInclude),
  tar_target(preferredReports, preferredReportsToInclude),
  tar_target(outputDirectory, {
    dir.create(outputDir, showWarnings = F, recursive = T)
    outputDir
  }),
  tar_target(dataDirectory, dataDir),
  tar_target(
    datasetSummary,
    readDatasetSummary(
      file.path(dataDir, "dataset-curation-summary.xlsx"),
      sheet="Sheet1",
      range=datasetSummaryRange,
      columnNames=datasetSummaryColumnNames
    )
  ),
  tar_target(
    datasetSummaryPlot,
    plotDatasetSummary(
      datasetSummary,
      standardReports,
      preferredReports,
      theme = mythemeXRot,
      outputDirectory
    )
  ),
  tar_target(ceramideColumnNames, ceramideColNames),
  tar_target(reportsDirectory, reportsDir),
  tar_target(naValueDefaults, naValues),
  tar_target(blankTypes, blankTypesTable),
  tar_target(
    intraAssayTable,
    readAssayTables(
      standardReports,
      preferredReports,
      ceramideColumnNames,
      reportsDirectory,
      outputDirectory,
      naValueDefaults,
      blankTypes
    )
  ),
  tar_target(
    intraAssayTableOutput,
    readr::write_csv(intraAssayTable, file.path(outputDirectory, "intraAssayTableFiltered.csv"))
  ),
  ################################################################################
  # Stop execution if we have any unassigned instruments
  ################################################################################
  tar_target(
    checkIntraAssayTable,
    stopifnot(nrow(intraAssayTable |> filter(Instrument %in% ("Unknown")))==0)
  ),
  ################################################################################
  # Subset data and survey plots to visually spot potential issues
  ################################################################################
  tar_target(
    intraAssayCalibrationLines,
    intraAssayTable |> filter(SampleType %in% c("Calibration Line 1", "Calibration Line 2"))
  ),
  tar_target(
    calibrationLineSurveyPlotObj,
    calibrationLineSurveyPlot(
      intraAssayCalibrationLines,
      theme = mythemeXRot,
      colorscale = mycolorscale,
      outlierShape = outlierShape,
      outputDirectory
    )
  ),
  tar_target(
    intraAssayQC,
    intraAssayTable |>
      filter(grepl("QC", SampleType)) |>
      mutate(SampleType = forcats::as_factor(as.character(SampleType)))
  ),
  tar_target(
    qcSurveyPlotObj,
    qcSurveyPlot(
      intraAssayQC,
      theme = mytheme,
      colorscale = mycolorscale,
      outlierShape = outlierShape,
      outputDirectory
    )
  ),
  tar_target(intraAssayNIST,
             intraAssayTable |> filter(grepl("NIST", SampleType))),
  tar_target(
    nistSurveyPlotObj,
    nistSurveyPlot(
      intraAssayNIST,
      theme = mythemeXRot,
      colorscale = mycolorscale,
      outlierShape = outlierShape,
      outputDirectory
    )
  ),
  ################################################################################
  # Calibration Line: Combine measurement data with STD concentrations
  ################################################################################
  tar_target(
    calibLinesWithStandardsConcs,
    prepareCalibLinesWithStandardConcs(
      intraAssayTable,
      expectedStdsConcentrations,
      outputDirectory
    )
  ),
  tar_target(
    cunit,
    unique(calibLinesWithStandardsConcs$Unit)
  ),
  tar_target(
    ceramideIds,
    unique(calibLinesWithStandardsConcs$ceramideId)
  ),
  tar_target(
    linRegFormula,
    y ~ x
  ),
  tar_target(
    scatterRatioPlotsObj,
    # Plots the measured ratio of ceramide and matching internal standard 
    # against the expected concentration of the ceramide
    scatterRatioPlot(
      linRegFormula,
      calibLinesWithStandardsConcs,
      ceramideIds,
      outputDirectory
    ),
    pattern = map(ceramideIds)
  ),
  ################################################################################
  # Calibration Line: Regression Model Fitting
  ################################################################################
  # select calibration line samples as defined in definitions.
  tar_target(
    calibLineData,
    {
      print(
        paste(
          "Calculating calibration line models for LabIds",
          paste(sort(unique(calibLinesWithStandardsConcs$LabId)), collapse = ", "),
          "and calibration line samples",
          paste(calibrationLineSamplesToUse, collapse = ", ")
        )
      )
      calibLinesWithStandardsConcs |>
        filter(Sample %in% calibrationLineSamplesToUse) |>
        ungroup() |>
        nest_by(LabId, SampleType, ceramideId, ceramideName) 
    }
  ),
  tar_target(
    # calculate linear models
    calibLineDataLm,
    calibLineData |>
      dplyr::mutate(CalibrationLine = list(
        lm(
          RatioLipidToIS ~ Concentration,
          weights = 1 / Concentration ^ 2,
          data = data
        ))) |>
      dplyr::mutate(CalibrationLine_C16 = list(
          lm(
            RatioLipidToIS_C16 ~ Concentration,
            weights = 1 / Concentration ^ 2,
            data = data
          )
      )) |>
      dplyr::mutate(CalibrationLine_C18 = list(
        lm(
          RatioLipidToIS_C18 ~ Concentration,
          weights = 1 / Concentration ^ 2,
          data = data
        )
      )) |>
      dplyr::mutate(CalibrationLine_C24 = list(
        lm(
          RatioLipidToIS_C24 ~ Concentration,
          weights = 1 / Concentration ^ 2,
          data = data
        )
      )) |>
      dplyr::mutate(CalibrationLine_C241 = list(
        lm(
          RatioLipidToIS_C241 ~ Concentration,
          weights = 1 / Concentration ^ 2,
          data = data
        )))
  ),
  tar_target(
    calibLineDataLmCoeffs,
    calibLineDataLm |> summarise(broom::tidy(CalibrationLine))
  ),
  tar_target(
    calibLineDataLmSumm,
    calibLineDataLm |> summarise(broom::glance(CalibrationLine))
  ),
  tar_target(
    calibLineDataLmPred,
    calibLineDataLm |> summarise(
      broom::augment(CalibrationLine, se_fit = TRUE, interval = "confidence")
    )
  ),
  tar_target(
    calibLineDataLmCoeffs_C16,
    calibLineDataLm |> summarise(broom::tidy(CalibrationLine_C16))
  ),
  tar_target(
    calibLineDataLmSumm_C16,
    calibLineDataLm |> summarise(broom::glance(CalibrationLine_C16))
  ),
  tar_target(
    calibLineDataLmPred_C16,
    calibLineDataLm |> summarise(
      broom::augment(CalibrationLine_C16, se_fit = TRUE, interval = "confidence")
    )
  ),
  tar_target(
    calibLineDataLmCoeffs_C18,
    calibLineDataLm |> summarise(broom::tidy(CalibrationLine_C18))
  ),
  tar_target(
    calibLineDataLmSumm_C18,
    calibLineDataLm |> summarise(broom::glance(CalibrationLine_C18))
  ),
  tar_target(
    calibLineDataLmPred_C18,
    calibLineDataLm |> summarise(
      broom::augment(CalibrationLine_C18, se_fit = TRUE, interval = "confidence")
    )
  ),
  tar_target(
    calibLineDataLmCoeffs_C24,
    calibLineDataLm |> summarise(broom::tidy(CalibrationLine_C24))
  ),
  tar_target(
    calibLineDataLmSumm_C24,
    calibLineDataLm |> summarise(broom::glance(CalibrationLine_C24))
  ),
  tar_target(
    calibLineDataLmPred_C24,
    calibLineDataLm |> summarise(
      broom::augment(CalibrationLine_C24, se_fit = TRUE, interval = "confidence")
    )
  ),

  tar_target(
    calibLineDataLmCoeffs_C241,
    calibLineDataLm |> summarise(broom::tidy(CalibrationLine_C241))
  ),
  tar_target(
    calibLineDataLmSumm_C241,
    calibLineDataLm |> summarise(broom::glance(CalibrationLine_C241))
  ),
  tar_target(
    calibLineDataLmPred_C241,
    calibLineDataLm |> summarise(
      broom::augment(CalibrationLine_C241, se_fit = TRUE, interval = "confidence")
    )
  ),
  tar_target(
    # signed log scales
    slog,
    scales::trans_new(
      "signed_log",
      transform = function(x)
        sign(x) * log(abs(x)),
      inverse = function(x)
        sign(x) * exp(abs(x))
    )
  ),
  tar_target(
    d_theoratios,
    read_csv(here::here(
      file.path(dataDir, "definitions", "Ceramide_STD_ISTD_TheoreticalRatios.csv")
    )) |> mutate(ceramideId = as.character(.data$ceramideId))
  ),
  tar_target(
    meanSdCalibrationLinesPlotsObj,
    meanSdCalibrationLinesPlots(
      calibLineDataLmPred,
      mytheme,
      avgLipidToISRatios = NULL,
      suffix = "",
      cunit = cunit,
      d_theoratios = d_theoratios,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    calibLineDataLmCoeffsWide,
    { 
      calibLineDataLmCoeffs |>
        select(-statistic, -p.value) |>
        group_by(LabId, SampleType, ceramideId, ceramideName, estimate) |>
        pivot_wider(names_from = term, values_from = c(estimate, std.error)) |>
        rename(
          Intercept = `estimate_(Intercept)`,
          Intercept_SE = `std.error_(Intercept)`,
          SlopeX = `estimate_Concentration`,
          SlopeX_SE = `std.error_Concentration`
        ) |>
        left_join(
          calibLineDataLmSumm |> select(sigma),
          by = c("LabId", "SampleType", "ceramideId", "ceramideName")
        ) 
    }
  ),
  tar_target(
    combinedCalibLinesWithConcs,
    { 
      browser()
      calibLineDataLm |> unnest(cols = c(data)) |>
        left_join(
          calibLineDataLmCoeffsWide,
          by = c("LabId", "SampleType", "ceramideId", "ceramideName")
        )
    }
  ),
  tar_target(
    calibLineDataLmCoeffsWide_C16,
    calibLineDataLmCoeffs_C16 |>
      select(-statistic, -p.value) |>
      group_by(LabId, SampleType, ceramideId, ceramideName, estimate) |>
      pivot_wider(names_from = term, values_from = c(estimate, std.error)) |>
      rename(
        Intercept_C16 = `estimate_(Intercept)`,
        Intercept_SE_C16 = `std.error_(Intercept)`,
        SlopeX_C16 = `estimate_Concentration`,
        SlopeX_C16_SE = `std.error_Concentration`
      ) |>
      left_join(
        calibLineDataLmSumm_C16 |> select(sigma_C16 = sigma),
        by = c("LabId", "SampleType", "ceramideId", "ceramideName")
      )
  ),
  tar_target(
    combinedCalibLinesWithConcs_C16,
    calibLineDataLm |> unnest(cols = c(data)) |>
      left_join(
        calibLineDataLmCoeffsWide_C16,
        by = c("LabId", "SampleType", "ceramideId", "ceramideName")
      )
  ),
  tar_target(
    calibLineDataLmCoeffsWide_C18,
    calibLineDataLmCoeffs_C18 |>
      select(-statistic, -p.value) |>
      group_by(LabId, SampleType, ceramideId, ceramideName, estimate) |>
      pivot_wider(names_from = term, values_from = c(estimate, std.error)) |>
      rename(
        Intercept_C18 = `estimate_(Intercept)`,
        Intercept_SE_C18 = `std.error_(Intercept)`,
        SlopeX_C18 = `estimate_Concentration`,
        SlopeX_C18_SE = `std.error_Concentration`
      ) |>
      left_join(
        calibLineDataLmSumm_C18 |> select(sigma_C18 = sigma),
        by = c("LabId", "SampleType", "ceramideId", "ceramideName")
      )
  ),
  tar_target(
    combinedCalibLinesWithConcs_C18,
    calibLineDataLm |> unnest(cols = c(data)) |>
      left_join(
        calibLineDataLmCoeffsWide_C18,
        by = c("LabId", "SampleType", "ceramideId", "ceramideName")
      )
  ),
  tar_target(
    calibLineDataLmCoeffsWide_C24,
    calibLineDataLmCoeffs_C24 |>
      select(-statistic, -p.value) |>
      group_by(LabId, SampleType, ceramideId, ceramideName, estimate) |>
      pivot_wider(names_from = term, values_from = c(estimate, std.error)) |>
      rename(
        Intercept_C24 = `estimate_(Intercept)`,
        Intercept_SE_C24 = `std.error_(Intercept)`,
        SlopeX_C24 = `estimate_Concentration`,
        SlopeX_C24_SE = `std.error_Concentration`
      ) |>
      left_join(
        calibLineDataLmSumm_C24 |> select(sigma_C24 = sigma),
        by = c("LabId", "SampleType", "ceramideId", "ceramideName")
      )
  ),
  tar_target(
    combinedCalibLinesWithConcs_C24,
    calibLineDataLm |> unnest(cols = c(data)) |>
      left_join(
        calibLineDataLmCoeffsWide_C24,
        by = c("LabId", "SampleType", "ceramideId", "ceramideName")
      )
  ),
  tar_target(
    calibLineDataLmCoeffsWide_C241,
    calibLineDataLmCoeffs_C241 |>
      select(-statistic, -p.value) |>
      group_by(LabId, SampleType, ceramideId, ceramideName, estimate) |>
      pivot_wider(names_from = term, values_from = c(estimate, std.error)) |>
      rename(
        Intercept_C241 = `estimate_(Intercept)`,
        Intercept_SE_C241 = `std.error_(Intercept)`,
        SlopeX_C241 = `estimate_Concentration`,
        SlopeX_C241_SE = `std.error_Concentration`
      ) |>
      left_join(
        calibLineDataLmSumm_C241 |> select(sigma_C241 = sigma),
        by = c("LabId", "SampleType", "ceramideId", "ceramideName")
      )
  ),
  tar_target(
    combinedCalibLinesWithConcs_C241,
    calibLineDataLm |> unnest(cols = c(data)) |>
      left_join(
        calibLineDataLmCoeffsWide_C241,
        by = c("LabId", "SampleType", "ceramideId", "ceramideName")
      )
  ),
  ################################################################################
  # Calibration Curves: Calculate Analyte Concentrations using CC linear models from areas
  ################################################################################
  # calculate unknown concentration back from ratio of Non-labeled Cer/Labeled Cer, 
  # slope of model and intercept and single point using the corresponding
  # internal standard, calculate LoD and LoQ from Intercept standard error divided by Slope
  # this may be an overestimation in this case, since the calibration lines were fitted 
  # over a larger concentration range than typically recommended (>10x the minimum concentration),
  # however, matrix / reagent blanks were not available with enough replicates or values different from 0 
  # to calculate reliable estimates for LoD and LoQ based on the LoB. 
  tar_target(
    analyteConcentrationsFromCalibLines_Auth,
    combinedCalibLinesWithConcs |>
      mutate(
        ReagentBlank = 0,
        C_Adj = C_A_cal(
          S_A = RatioLipidToIS,
          a = SlopeX,
          b = Intercept - ReagentBlank
        ),
        C_LoD_C = (3.3 * Intercept_SE) / SlopeX,
        C_LoQ_C = (10 * Intercept_SE) / SlopeX,
        C_SinglePoint = (area_l / area_h * ISConcentration) - ReagentBlank,
        C_Adj_Rec = C_Adj/Concentration,
        C_Adj_Rec_Perc = C_Adj_Rec*100,
        C_SinglePoint_Rec = C_SinglePoint/Concentration,
        C_SinglePoint_Rec_Perc = C_SinglePoint_Rec*100
      )
  ),
  tar_target(
    analyteConcentrationsFromCalibLines_C16,
    combinedCalibLinesWithConcs_C16 |>
      mutate(
        ReagentBlank = 0,
        C_Adj_C16 = C_A_cal(
          S_A = RatioLipidToIS_C16,
          a = SlopeX_C16,
          b = Intercept_C16 - ReagentBlank
        )
      )
  ),
  tar_target(
    analyteConcentrationsFromCalibLines_C18,
    combinedCalibLinesWithConcs_C18 |>
      mutate(
        ReagentBlank = 0,
        C_Adj_C18 = C_A_cal(
          S_A = RatioLipidToIS_C18,
          a = SlopeX_C18,
          b = Intercept_C18 - ReagentBlank
        )
      )
  ),
  tar_target(
    analyteConcentrationsFromCalibLines_C24,
    combinedCalibLinesWithConcs_C24 |>
      mutate(
        ReagentBlank = 0,
        C_Adj_C24 = C_A_cal(
          S_A = RatioLipidToIS_C24,
          a = SlopeX_C24,
          b = Intercept_C24 - ReagentBlank
        )
      )
  ),
  tar_target(
    analyteConcentrationsFromCalibLines_C241,
    combinedCalibLinesWithConcs_C241 |>
      mutate(
        ReagentBlank = 0,
        C_Adj_C241 = C_A_cal(
          S_A = RatioLipidToIS_C241,
          a = SlopeX_C241,
          b = Intercept_C241 - ReagentBlank
        )
      )
  ),

  tar_target(
    analyteConcentrationsFromCalibLines,
    analyteConcentrationsFromCalibLines_function(
      analyteConcentrationsFromCalibLines_Auth,
      analyteConcentrationsFromCalibLines_C16,
      analyteConcentrationsFromCalibLines_C18,
      analyteConcentrationsFromCalibLines_C24,
      analyteConcentrationsFromCalibLines_C241
    )
  ),

  #stopifnot(manual_test_C_SinglePoint==analyteConcentrationsFromCalibLines[1,]$C_SinglePoint)
  tar_target(
    analyteConcentrationsFromCalibLinesFile,
    readr::write_csv(
      file = file.path(outputDirectory, "analyteConcentrationsFromCalibLines.csv"),
      analyteConcentrationsFromCalibLines
    )
  ),
  
  tar_target(
    analyteConcentrationsFromCalibLinesSummaryFile,
    readr::write_csv(
      file = file.path(outputDirectory, "LoD_LoQ_FromCalibLines.csv"),
      analyteConcentrationsFromCalibLines |> group_by(
        LabId, 
        SampleType, 
        ceramideId, 
        ceramideName,
        Intercept, 
        Intercept_SE, 
        SlopeX, 
        SlopeX_SE, 
        sigma, 
        C_LoD_C, 
        C_LoQ_C
      )
      |> summarise()
    )
  ),
  
  tar_target(
    linesLabelData,
    analyteConcentrationsFromCalibLines |>
      group_by(LabId, ceramideId, Sample) |>
      filter(Concentration == min(Concentration))
  ),
  tar_target(
    calibrationLineVsSinglePointPlotObjs,
    calibrationLineVsSinglePointPlot(
      analyteConcentrationsFromCalibLines,
      selectedCeramideId = ceramideIds,
      theme=mytheme,
      mycolorscale=mycolorscale,
      outputDirectory=outputDirectory
    ),
    pattern = map(ceramideIds)
  ),
  tar_target(
    calibrationLineVsSinglePointPlotObjsQQQ,
    calibrationLineVsSinglePointPlot(
      analyteConcentrationsFromCalibLines |> filter(MassAnalyzerType=="QQQ"),
      selectedCeramideId = ceramideIds,
      suffix="-QQQ",
      theme=mytheme,
      mycolorscale=mycolorscale,
      outputDirectory=outputDirectory
    ),
    pattern = map(ceramideIds)
  ),
  ################################################################################
  # Calibration Line: Plots for each lab
  ################################################################################
  tar_target(
    labIds,
    unique(calibLineDataLmPred$LabId)
  ),
  tar_target(
    calibLineLabPlotObj,
    calibrationLinePlot(
      linRegFormula = linRegFormula,
      lmPred = calibLineDataLmPred,
      adjAnalyteConcentrations = analyteConcentrationsFromCalibLines,
      theme = mytheme,
      labId = labIds,
      digits = 4,
      cunit,
      d_theoratios,
      outputDirectory
    ),
    pattern = map(labIds)
  ),
  tar_target(
    calibLineLabResidualsPlotObj,
    calibrationLineResidualsPlot(
      linRegFormula = linRegFormula,
      lmPred = calibLineDataLmPred,
      adjAnalyteConcentrations = analyteConcentrationsFromCalibLines,
      theme = mytheme,
      labId = labIds,
      digits = 4,
      cunit,
      d_theoratios,
      outputDirectory
    ),
    pattern = map(labIds)
  ),

  ################################################################################
  # Calibration Curves: Adjust concentrations as mean corrected values of Calibration Curves 1 and 2
  ################################################################################
  # Reference: Kauhanen et al., Development and validation of a high-throughput LC–MS/MS
  # assay for routine measurement of molecular ceramides, Anal Bioanal Chem (2016) 408:3475–3483
  # calculate averaged concentrations using the sum of the adjusted concentration using Calibration Line 1 and
  # the adjusted concentration using Calibration Line 2, divided by two
  # average adjusted concentrations -> (C1 + C2) / 2
  tar_target(
    averagedConcentrations,
    analyteConcentrationsFromCalibLines |>
      ungroup() |>
      select(
        LabId,
        SampleType,
        Sample,
        ceramideId,
        ceramideName,
        Unit,
        replicate,
        Concentration,
        ISConcentration,
        C_Adj,
        C_Adj_C16,
        C_Adj_C18,
        C_Adj_C24,
        C_Adj_C241,
        C_LoD_C,
        C_LoQ_C,
        C_SinglePoint
      ) |>
      group_by(LabId, Sample, ceramideId, ceramideName, Unit, replicate) |>
      pivot_wider(
        names_from = SampleType,
        values_from = c(C_Adj, C_Adj_C16, C_Adj_C18, C_Adj_C24, C_Adj_C241, C_LoD_C, C_LoQ_C, C_SinglePoint)
      ) |>
      # remove NAs that may have been introduced by grouping over replicates
      # some datasets only have two replicates at this point and pivot_wider may
      # introduce NAs for those to fill gaps
      filter(
        !is.na(`C_Adj_Calibration Line 1`) &
          !is.na(`C_Adj_Calibration Line 2`)
      ) |>
      mutate(
        Avg_C_Adj = (`C_Adj_Calibration Line 1` + `C_Adj_Calibration Line 2`) /
          2,
        Avg_C_Adj_C16 = (`C_Adj_C16_Calibration Line 1` + `C_Adj_C16_Calibration Line 2`) /
          2,
        Avg_C_Adj_C18 = (`C_Adj_C18_Calibration Line 1` + `C_Adj_C18_Calibration Line 2`) /
          2,
        Avg_C_Adj_C24 = (`C_Adj_C24_Calibration Line 1` + `C_Adj_C24_Calibration Line 2`) /
          2,
        Avg_C_Adj_C241 = (`C_Adj_C241_Calibration Line 1` + `C_Adj_C241_Calibration Line 2`) /
          2,
        Avg_C_LoD_C = (`C_LoD_C_Calibration Line 1` + `C_LoD_C_Calibration Line 2`) /
          2,
        Avg_C_LoQ_C = (`C_LoQ_C_Calibration Line 1` + `C_LoQ_C_Calibration Line 2`) /
          2,
        Avg_C_SinglePoint = (
          `C_SinglePoint_Calibration Line 1` + `C_SinglePoint_Calibration Line 2`
        ) / 2,
        Avg_C_Adj_Rec = Avg_C_Adj/Concentration,
        Avg_C_Adj_Rec_Perc = Avg_C_Adj_Rec*100,
        Avg_C_SinglePoint_Rec = Avg_C_SinglePoint/Concentration,
        Avg_C_SinglePoint_Rec_Perc = Avg_C_SinglePoint_Rec*100,
        label = paste(Sample, " ", Concentration, " ", Unit)
      )
  ),
  # boxplots of averaged adjusted concentrations
  tar_target(
    boxplotCer1Obj,
    stdsBoxplot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "1",
      theme = mythemeXRot,
      outlierShape=outlierShape,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    stdsVarHistoCer1Obj,
    stdsVarHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "1",
      theme = mythemeXRot,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    stdsRecoveryPercentHistoCer1Obj,
    stdsRecoveryPercentHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "1",
      theme = mythemeXRot,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    boxplotCer2Obj,
    stdsBoxplot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "2",
      theme = mythemeXRot,
      outlierShape=outlierShape,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    stdsVarHistoCer2Obj,
    stdsVarHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "2",
      theme = mythemeXRot,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    stdsRecoveryPercentHistoCer2Obj,
    stdsRecoveryPercentHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "2",
      theme = mythemeXRot,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    boxplotCer3Obj,
    stdsBoxplot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "3",
      theme = mythemeXRot,
      outlierShape=outlierShape,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    stdsVarHistoCer3Obj,
    stdsVarHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "3",
      theme = mythemeXRot,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    stdsRecoveryPercentHistoCer3Obj,
    stdsRecoveryPercentHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "3",
      fillscale = myfillscale,
      theme = mythemeXRot,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    boxplotCer4Obj,
    stdsBoxplot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "4",
      theme = mythemeXRot,
      outlierShape=outlierShape,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    stdsVarHistoCer4Obj,
    stdsVarHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "4",
      theme = mythemeXRot,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    stdsRecoveryPercentHistoCer4Obj,
    stdsRecoveryPercentHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "4",
      theme = mythemeXRot,
      fillscale = myfillscale,
      outputDirectory = outputDirectory
    )
  ),
  ################################################################################
  # QC Samples figures of merit calculation
  ################################################################################
  tar_target(
    qcAveragedConcentrations,
    createQcAveragedConcentrations(
      intraAssayQC = intraAssayQC,
      expectedStdsConcentrations = expectedStdsConcentrations,
      expectedCalibrationLineConcentrations = expectedCalibrationLineConcentrations,
      calibLineDataLmCoeffsWide = calibLineDataLmCoeffsWide
    )
  ),
  tar_target(
    qcAveragedConcentrationsFile,
    readr::write_csv(
      qcAveragedConcentrations,
      file = file.path(outputDirectory, "qcCal1Cal2AveragedConcentrations.csv")
    )
  ),
  tar_target(
    qcAveragedConcentrationsLong,
    createQcAveragedConcentrationsLong(qcAveragedConcentrations)
  ),
  tar_target(
    qcConcentrationsPlotObj,
    qcConcentrationsPlot(
      qcAveragedConcentrationsLong,
      theme = mytheme,
      fillscale = myfillscale,
      outlierShape = 19,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    qcConcentrationsByLabIdPlotObj,
    qcConcentrationsByLabIdPlot(
      qcAveragedConcentrations,
      theme = mytheme,
      fillscale = myfillscale,
      outlierShape = 19,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    qcConcentrationsByLabIdDensityPlotObj,
    qcConcentrationsByLabIdDensityPlot(
      qcAveragedConcentrations,
      theme = mytheme,
      fillscale = myfillscale,
      outlierShape = 19,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    qcSampleStats,
    createQcConcentrationsStats(qcAveragedConcentrations)
  ),
  tar_target(
    qcSampleStatsPlotObj,
    qcSampleStatsPlot(qcSampleStats,
                      theme = mytheme,
                      fillscale = myfillscale,
                      outputDirectory = outputDirectory)
  ),
  ################################################################################
  # NIST Samples concentration calculation
  ################################################################################
  tar_target(
    IntraAssayNISTwide,
    {intraAssayNIST_Part1 <- intraAssayNIST |>
      left_join(
        expectedStdsConcentrations |> filter(Sample == "STD 1") |> select(-Sample),
        by = c("ceramideId", "ceramideName")
      ) |>
      mutate(Unit = unique(
        expectedCalibrationLineConcentrations$Unit
      ))

    intraAssayNIST_Part2 <- intraAssayNIST_Part1 |>
      group_by(
        LabId,
        SampleType,
        Sample,
        ceramideId,
        ceramideName,
        Unit,
        replicate
      ) |>
      mutate(
        RatioLipidToIS = area[isotope == "l"] / area[isotope == "h"],
        PAR = 1 / RatioLipidToIS
      ) |>
      pivot_wider(
        names_from = isotope,
        values_from = area,
        names_prefix = "area_"
      ) |> ungroup() |>
      # here, Samples are independent analytical replicates, while replicates
      # are repeated injections of the same sample
      group_by(LabId, SampleType, ceramideId, ceramideName, Unit) |>
      filter(n() >= 3 &
              (!is.na(area_l) |
                 !is.na(area_h))) # remove all NA entries and groups with less than three entries

    intraAssayNIST_Part3 <- intraAssayNIST_Part1 |>
    group_by(
      LabId,
      SampleType,
      Sample,
      Unit,
      replicate
    ) |>
    mutate(
      RatioLipidToIS_C16 = area / area[isotope == "h" & ceramideName == "Cer 18:1;O2/16:0"],
      RatioLipidToIS_C18 = area / area[isotope == "h" & ceramideName == "Cer 18:1;O2/18:0"],
      RatioLipidToIS_C24 = area / area[isotope == "h" & ceramideName == "Cer 18:1;O2/24:0"],
      RatioLipidToIS_C241 = area / area[isotope == "h" & ceramideName == "Cer 18:1;O2/24:1"]
    ) |>
    filter(isotope == "l") |>
    pivot_wider(
        names_from = isotope,
        values_from = area,
        names_prefix = "area_"
      ) |> ungroup()
     intraAssayNIST_Part2 |>
     left_join(
         intraAssayNIST_Part3 |> dplyr::select(
           Sample,
           LabId,
           ceramideName,
           replicate,
           RatioLipidToIS_C16,
           RatioLipidToIS_C18,
           RatioLipidToIS_C24,
           RatioLipidToIS_C241
         ),
         by = c("Sample", "LabId", "ceramideName", "replicate")
     )}
    ## ------ Bo

  ),
  tar_target(
    IntraAssayNISTwideFile,
    readr::write_csv(file = file.path(outputDirectory, "IntraAssayNISTwideFile.csv"),
                     x = IntraAssayNISTwide)
  ),
  tar_target(
    combinedIntraAssayNISTWithLms,
    IntraAssayNISTwide |> ungroup() |> group_by(LabId, ceramideId, ceramideName) |>
      left_join(
        calibLineDataLmCoeffsWide |>
          ungroup() |>
          mutate(CalibrationLine = SampleType) |> select(-SampleType),
        by = c("LabId", "ceramideId", "ceramideName")
      ) |>
      left_join(
        calibLineDataLmCoeffsWide_C16 |>
          ungroup() |>
          mutate(CalibrationLine = SampleType) |> select(-SampleType),
        by = c("LabId", "ceramideId", "ceramideName", "CalibrationLine")
      ) |>
      left_join(
        calibLineDataLmCoeffsWide_C18 |>
          ungroup() |>
          mutate(CalibrationLine = SampleType) |> select(-SampleType),
        by = c("LabId", "ceramideId", "ceramideName", "CalibrationLine")
      ) |>
      left_join(
        calibLineDataLmCoeffsWide_C24 |>
          ungroup() |>
          mutate(CalibrationLine = SampleType) |> select(-SampleType),
        by = c("LabId", "ceramideId", "ceramideName", "CalibrationLine")
      ) |>
      left_join(
        calibLineDataLmCoeffsWide_C241 |>
          ungroup() |>
          mutate(CalibrationLine = SampleType) |> select(-SampleType),
        by = c("LabId", "ceramideId", "ceramideName", "CalibrationLine")
      )
  ),
  tar_target(
    manual_nist_test_C_SinglePoint,
    2 * combinedIntraAssayNISTWithLms[1, ]$area_l / combinedIntraAssayNISTWithLms[1, ]$area_h *
      combinedIntraAssayNISTWithLms[1, ]$ISConcentration
  ),
  tar_target(
    nistAnalyteConcentrationsFromCalibLines,
    combinedIntraAssayNISTWithLms |>
      mutate(
        C_Adj = C_A_cal(S_A = 2* RatioLipidToIS, a = SlopeX, b = Intercept),
        C_Adj_C16 = C_A_cal(S_A = 2*RatioLipidToIS_C16, a = SlopeX_C16, b = Intercept_C16),
        C_Adj_C18 = C_A_cal(S_A = 2*RatioLipidToIS_C18, a = SlopeX_C18, b = Intercept_C18),
        C_Adj_C24 = C_A_cal(S_A = 2*RatioLipidToIS_C24, a = SlopeX_C24, b = Intercept_C24),
        C_Adj_C241 = C_A_cal(S_A = 2*RatioLipidToIS_C241, a = SlopeX_C241, b = Intercept_C241),
        #ISConcentration is only half as large in the QC and NIST samples, compensate by multiplying by 2
        C_SinglePoint = 2 * area_l / area_h * ISConcentration
      )
  ),
  tar_target(
    nistAnalyteConcentrationsFromCalibLinesFile,
    readr::write_csv(
      file = file.path(
        outputDirectory,
        "nistAnalyteConcentrationsFromCalibLines.csv"
      ),
      x = nistAnalyteConcentrationsFromCalibLines
    )
  ),
  tar_target(
    stopAtSinglePoint,
    stopifnot(
      manual_nist_test_C_SinglePoint == nistAnalyteConcentrationsFromCalibLines[1, ]$C_SinglePoint
    )
  ),
  tar_target(
    nistSampleTypes,
    unique(nistAnalyteConcentrationsFromCalibLines$SampleType)
  ),
  tar_target(
    nistAreaRatioPlotsObj,
    nistAreaRatioPlot(
      data = nistAnalyteConcentrationsFromCalibLines,
      selectedSampleType = nistSampleTypes,
      theme = mythemeXRot,
      fillscale = myfillscale,
      outlierShape = outlierShape,
      outputDirectory = outputDirectory
    ),
    pattern = map(nistSampleTypes)
  ),
  tar_target(
    nistSrm1950MeanAreaRatios,
    nistAnalyteConcentrationsFromCalibLines |>
      filter(SampleType == "NIST SRM") |>
      ungroup() |>
      group_by(ceramideId, ceramideName) |>
      summarise(
        meanAreaRatios = mean(RatioLipidToIS),
        sdAreaRatios = sd(RatioLipidToIS)
      )
  ),
  tar_target(
    meanSdCalibrationLinesPlotsAvgRatioObj,
    meanSdCalibrationLinesPlots(
      calibLineDataLmPred = calibLineDataLmPred,
      theme = mytheme,
      avgLipidToISRatios = nistSrm1950MeanAreaRatios,
      suffix = "-with-NISTSRM-AvgRatio",
      cunit = cunit,
      d_theoratios,
      outputDirectory = outputDirectory
    )
  ),
  tar_target(
    nistAveragedConcentrations,
    get_nistAveragedConcentrations(nistAnalyteConcentrationsFromCalibLines)
  ),
  ################################################################################
  # NIST Sample Details
  ################################################################################
  tar_target(
    massAnalyzerSummary,
    nistAveragedConcentrations |>
      ungroup() |>
      group_by(LabId) |>
      summarise(MassAnalyzerType = first(MassAnalyzerType)) |>
      group_by(MassAnalyzerType) |>
      summarise(n_massAnalyzer = n()) |>
      mutate(facetLabel = paste0(MassAnalyzerType, " n=", n_massAnalyzer))
  ),
  tar_target(
    nistConcentrationsStats,
    nistAveragedConcentrations |>
      group_by(
        LabId,
        SampleType,
        ceramideId,
        ceramideName,
        Unit,
        Protocol,
        Instrument,
        MassAnalyzerType
      ) |>
      summarise(
        AvgAvg_C_Adj = mean(Avg_C_Adj),
        SdAvg_C_Adj = sd(Avg_C_Adj),
        CV_C_Adj = SdAvg_C_Adj / AvgAvg_C_Adj,
        CV_Perc = CV_C_Adj * 100,
        ZScore_C_Adj = (Avg_C_Adj - AvgAvg_C_Adj) / SdAvg_C_Adj,
        Avg_C_SinglePoint = mean(C_SinglePoint),
        SdAvg_C_SinglePoint = sd(C_SinglePoint),
        CV_C_SinglePoint = SdAvg_C_SinglePoint / Avg_C_SinglePoint,
        CV_Perc_C_SinglePoint = CV_C_SinglePoint * 100,
        ZScore_C_SinglePoint = (C_SinglePoint - Avg_C_SinglePoint) / SdAvg_C_SinglePoint,
      ) |> left_join(massAnalyzerSummary, by = c("MassAnalyzerType"))
  ),
  tar_target(
    nistConcentrationsStatsLong,
    nistConcentrationsStats |>
      pivot_longer(
        cols = c("CV_C_Adj", "CV_C_SinglePoint"),
        names_to = "Calibration",
        values_to = "CV"
      )
  ),
  tar_target(
    nistConcentrationsZScoreStatsLong,
    nistConcentrationsStats |>
      pivot_longer(
        cols = c("ZScore_C_Adj", "ZScore_C_SinglePoint"),
        names_to = "Calibration",
        values_to = "ZScore"
      )
  ),
  tar_target(
    nistAveragedConcentrationsLong,
    nistAveragedConcentrations |>
      mutate(
        C_CalibCurve_Avg = Avg_C_Adj,
        C_CalibCurve_Avg_C16 = Avg_C_Adj_C16,
        C_CalibCurve_Avg_C18 = Avg_C_Adj_C18,
        C_CalibCurve_Avg_C24 = Avg_C_Adj_C24,
        C_CalibCurve_Avg_C241 = Avg_C_Adj_C241
      ) |>
      pivot_longer(
        cols = c(
          "C_CalibCurve_Avg",
          "C_CalibCurve_Avg_C16",
          "C_CalibCurve_Avg_C18",
          "C_CalibCurve_Avg_C24",
          "C_CalibCurve_Avg_C241",
          "C_SinglePoint"
        ),
        names_to = "Calibration",
        values_to = "Adj_Conc"
      ) |> left_join(massAnalyzerSummary, by = c("MassAnalyzerType"))
  ),
  tar_target(
    nistAveragedConcentrationsLongFile,
    write_csv(
      nistAveragedConcentrationsLong,
      file = file.path(
        outputDirectory,
        "nistCal1Cal2AveragedConcentrationsLong.csv"
      )
    )
  ),
  tar_target(
    nistAveragedConcentrationsFile,
    write_csv(
      nistAveragedConcentrations,
      file = file.path(outputDirectory, "nistCal1Cal2AveragedConcentrations.csv")
    )
  ),
  ################################################################################
  # NIST Sample and Ceramide Concentrations Plots
  ################################################################################
  tar_target(
    nistConcentrationsPlotObj,
    nistConcentrationsPlot(
      nistAveragedConcentrationsLong,
      fillscale = myfillscale,
      theme = mythemeXRot,
      outlierShape = outlierShape,
      outputDirectory = outputDirectory
    )
  ),
  ################################################################################
  # NIST Sample and Ceramide Cv Plots
  ################################################################################
  tar_target(
    nistCvPlotObj,
    nistCvPlot(
      nistConcentrationsStatsLong,
      fillscale = myfillscale,
      theme = mythemeXRot,
      outlierShape = outlierShape,
      outputDirectory = outputDirectory
    )
  ),
  ################################################################################
  # NIST Concentrations by Sample type, LabId and Mass Analyzer Plots
  ################################################################################
  tar_target(
    nistAveragedConcentrationsLongPlotObjs,
    nistAveragedConcentrationsLong |>
      #filter(LabId=="15") |>
      ungroup() |>
      group_by(SampleType, Unit) |>
      group_walk(
        ~ nistConcentrationsPlotBySampleType(
          data = .x,
          namesuffix = unique(.y$SampleType),
          unit = unique(nistAveragedConcentrations$Unit),
          theme = mythemeXRot,
          xfillscale = myfillscale,
          outlierShape=outlierShape,
          outputDirectory = outputDirectory
        ),
        .keep = TRUE
      ) |>
      group_walk(
        ~ nistConcentrationsPlotByMa(
          data = .x,
          namesuffix = unique(.y$SampleType),
          unit = unique(nistAveragedConcentrations$Unit),
          theme = mythemeXRot,
          xfillscale = myfillscale,
          outlierShape=outlierShape,
          outputDirectory = outputDirectory
        ),
        .keep = TRUE
      )
  ),
  ################################################################################
  # NIST Concentrations by LabId Plots
  ################################################################################
  tar_target(
    nistConcentrationsByLabIdPlotObj,
    nistConcentrationsByLabIdPlot(
      nistAveragedConcentrationsLong,
      theme = mythemeXRot,
      xfillscale = myfillscale,
      outlierShape=outlierShape,
      outputDirectory = outputDirectory
    )
  ),
  ################################################################################
  # NIST CV by LabId Plots
  ################################################################################
  tar_target(
    nistCvByLabIdPlotObj,
    nistCvByLabIdPlot(
      nistConcentrationsStatsLong,
      theme = mythemeXRot,
      xfillscale = myfillscale,
      outlierShape=outlierShape,
      outputDirectory = outputDirectory
    )
  ),
  ################################################################################
  # NIST Zscores by LabId Plots
  ################################################################################
  tar_target(
    nistZscoreByLabIdPlotObj,
    nistZscoreByLabIdPlot(
      nistConcentrationsZScoreStatsLong,
      theme = mythemeXRot,
      xfillscale = myfillscale,
      outlierShape=outlierShape,
      outputDirectory = outputDirectory
    )
  ),
  ################################################################################
  # Render Rmarkdown report to generate manuscript tables and plots
  ################################################################################  
  tarchetypes::tar_render(
    name = manuscriptRmarkdown,
    path = "manuscript/manuscript-figures-tables.Rmd",
    params = list(run_in_targets = TRUE)
  )
)
  

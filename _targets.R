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
    "here"
  )
)

tar_option_set(debug = "calibLineDataLmAnova")

reportsDir <- "./data"
Sys.setenv("REPORTS_DIR"=reportsDir)

outputDir <- "./output"
Sys.setenv("OUTPUT_DIR"=outputDir)

source("R/definitions.R")
source("R/functions_io.R")
source("R/functions_filterData.R")
source("R/functions_steps.R")
source("R/functions_methods.R")
source("R/functions_ggpubr-custom.R")
source("R/functions_plots.R")

list(
  tar_target(createOutputDir, dir.create(outputDir, recursive = T, showWarnings = F)),
  tar_target(standardReports, standardReportsToInclude),
  tar_target(preferredReports, preferredReportsToInclude),
  tar_target(
    datasetSummary, 
    readDatasetSummary(
      file.path(reportsDir, "dataset-curation-summary.xlsx"),
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
      theme = mythemeXRot
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
      file.path(reportsDirectory, "reports"),
      naValueDefaults,
      blankTypes
    )
  ),
  ################################################################################
  # Subset data and survey plots
  ################################################################################
  tar_target(
    intraAssayCalibrationLines,
    intraAssayTable |> filter(SampleType %in% c("Calibration Line 1", "Calibration Line 2"))
  ),
  tar_target(
    calibrationLineSurveyPlotObj,
    calibrationLineSurveyPlot(intraAssayCalibrationLines, theme=mythemeXRot, colorscale=mycolorscale, outlierShape=outlierShape)
  ),
  tar_target(
    intraAssayQC,
    intraAssayTable |>
      filter(grepl("QC", SampleType)) |>
      mutate(SampleType = forcats::as_factor(as.character(SampleType)))
  ),
  tar_target(
    qcSurveyPlotObj,
    qcSurveyPlot(intraAssayQC, theme=mytheme, colorscale=mycolorscale, outlierShape=outlierShape)
  ),
  tar_target(intraAssayNIST,
             intraAssayTable |> filter(grepl("NIST", SampleType))),
  tar_target(
    nistSurveyPlotObj,
    nistSurveyPlot(intraAssayNIST, theme=mythemeXRot, colorscale=mycolorscale, outlierShape=outlierShape)
  ),
  ################################################################################
  # Calibration Line: Combine measurement data with STD concentrations
  ################################################################################
  tar_target(
    calibLinesWithStandardsConcs,
    prepareCalibLinesWithStandardConcs(intraAssayTable,
                                       expectedStdsConcentrations)
  ),
  tar_target(cunit,
             unique(calibLinesWithStandardsConcs$Unit)),
  # tar_assert_expr(
  #   length(cunit)==1,
  #   "Expected only one concentration unit!"
  # ),
  tar_target(
    ceramideIds,
    unique(calibLinesWithStandardsConcs$ceramideId)
  ),
  tar_target(linRegFormula,
             y ~ x),
  tar_target(
    scatterRatioPlotsObj,
    # Plots the measured ratio of ceramide and matching internal standard against the expected concentration of the ceramide
    scatterRatioPlots(linRegFormula, calibLinesWithStandardsConcs, ceramideIds)
  ),
  ################################################################################
  # Calibration Line: Regression Model Fitting
  ################################################################################
  # report cases with too few replicates -> still go through model calculation, but may have "bad" models
  tar_target(
    calibLinesNotEnoughData,
    intraAssayTable |>
      filter(grepl("Calibration Line", SampleType)) |>
      nest_by(LabId, SampleType, ceramideId, ceramideName)
  ),
  tar_target(
    labResultsWithoutCalibrationLines,
    nrow(calibLinesNotEnoughData)
  ),
  tar_target(
    calibLineData,
    calibLinesWithStandardsConcs |>
      filter(Sample %in% calibrationLineSamplesToUse) |>
      ungroup() |>
      nest_by(LabId, SampleType, ceramideId, ceramideName)
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
        )
      ))
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
    meanSdCalibrationLinesPlotsObj,
    meanSdCalibrationLinesPlots(
      calibLineDataLmPred,
      mytheme,
      avgLipidToISRatios = NULL,
      suffix = "",
      cunit = cunit
    )
  ),
  tar_target(
    calibLineDataLmCoeffsWide,
    calibLineDataLmCoeffs |>
      select(-std.error, -statistic, -p.value) |>
      group_by(LabId, SampleType, ceramideId, ceramideName, estimate) |>
      pivot_wider(names_from = term, values_from = estimate) |>
      rename(Intercept = `(Intercept)`, SlopeX = Concentration) |>
      left_join(
        calibLineDataLmSumm |> select(sigma),
        by = c("LabId", "SampleType", "ceramideId", "ceramideName")
      )
  ),
  tar_target(
    combinedCalibLinesWithConcs,
    calibLineDataLm |> unnest(cols = c(data)) |>
      left_join(
        calibLineDataLmCoeffsWide,
        by = c("LabId", "SampleType", "ceramideId", "ceramideName")
      )
  ),
  ################################################################################
  # Calibration Curves: Calculate Analyte Concentrations using CC linear models from areas
  ################################################################################
  # calculate unknown concentration back from ratio of Light Cer/Heavy Cer, slope of model and intercept
  # and single point using the corresponding internal standard
  tar_target(
    analyteConcentrationsFromCalibLines,
    combinedCalibLinesWithConcs |>
      mutate(
        ReagentBlank = 0,
        C_Adj = C_A_cal(
          S_A = area_l / area_h,
          a = SlopeX,
          b = Intercept - ReagentBlank
        ),
        C_LoD_C = (3.3 * sdev_area_ratio) / SlopeX,
        C_LoQ_C = (10 * sdev_area_ratio) / SlopeX,
        C_SinglePoint = (area_l / area_h * ISConcentration) - ReagentBlank,
        C_Adj_Rec = C_Adj/Concentration,
        C_Adj_Rec_Perc = C_Adj_Rec*100,
        C_SinglePoint_Rec = C_SinglePoint/Concentration,
        C_SinglePoint_Rec_Perc = C_SinglePoint_Rec*100
      )
  ),
  tar_target(
    analyteConcentrationsFromCalibLinesFile,
    readr::write_csv(file = "output/analyteConcentrationsFromCalibLines.csv", analyteConcentrationsFromCalibLines)
  ),
  tar_target(
    linesLabelData,
    analyteConcentrationsFromCalibLines |> group_by(LabId, ceramideId, Sample) |> filter(Concentration ==
                                                                                             min(Concentration))
  ),
  tar_target(
    calibrationLineVsSinglePointPlotObj,
    calibrationLineVsSinglePointPlot(analyteConcentrationsFromCalibLines, mytheme, mycolorscale)
  ),
  tar_target(
    calibrationAnalyteConcentrationsPlotObj,
    calibrationAnalyteConcentrationsPlot(
      analyteConcentrationsFromCalibLines,
      calibLineDataLmPred,
      mytheme,
      cunit,
      linRegFormula
    )
  ),
  ################################################################################
  # Calibration Line: Plots for each lab
  ################################################################################
  tar_target(labIds,
             unique(calibLineDataLmPred$LabId)),
  tar_target(
    calibLineLabPlotObj,
    labIds |>
      map_dfr(
        ~ calibrationLinePlot(
          linRegFormula = linRegFormula,
          lmPred = calibLineDataLmPred,
          adjAnalyteConcentrations = analyteConcentrationsFromCalibLines,
          theme = mytheme,
          labId = .x,
          digits = 4,
          cunit
        )
      )
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
        C_LoD_C,
        C_LoQ_C,
        C_SinglePoint
      ) |>
      group_by(LabId, Sample, ceramideId, ceramideName, Unit, replicate) |>
      pivot_wider(
        names_from = SampleType,
        values_from = c(C_Adj, C_LoD_C, C_LoQ_C, C_SinglePoint)
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
        Avg_C_LoD_C = (`C_LoD_C_Calibration Line 1` + `C_LoD_C_Calibration Line 2`) /
          2,
        Avg_C_LoQ_C = (`C_LoQ_C_Calibration Line 1` + `C_LoQ_C_Calibration Line 2`) /
          2,
        # Avg_C_Adj = (`C_Adj_Calibration Line 1`) /
        #   1,
        # Avg_C_LoD_C = (`C_LoD_C_Calibration Line 1`) /
        #   1,
        # Avg_C_LoQ_C = (`C_LoQ_C_Calibration Line 1`) /
        #   1,
        Avg_C_SinglePoint = (
          `C_SinglePoint_Calibration Line 1` + `C_SinglePoint_Calibration Line 2`
        ) / 2,
        # Avg_C_SinglePoint = (
        #   `C_SinglePoint_Calibration Line 1`
        # ) / 1,
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
      outlierShape=outlierShape
    ) + myfillscale
  ),
  tar_target(
    stdsVarHistoCer1Obj,
    stdsVarHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "1",
      theme = mythemeXRot
    ) + myfillscale
  ),
  tar_target(
    stdsRecoveryPercentHistoCer1Obj,
    stdsRecoveryPercentHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "1",
      theme = mythemeXRot
    ) + myfillscale
  ),
  tar_target(
    boxplotCer2Obj,
    stdsBoxplot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "2",
      theme = mythemeXRot,
      outlierShape=outlierShape
    ) + myfillscale
  ),
  tar_target(
    stdsVarHistoCer2Obj,
    stdsVarHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "2",
      theme = mythemeXRot
    ) + myfillscale
  ),
  tar_target(
    stdsRecoveryPercentHistoCer2Obj,
    stdsRecoveryPercentHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "2",
      theme = mythemeXRot
    ) + myfillscale
  ),
  tar_target(
    boxplotCer3Obj,
    stdsBoxplot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "3",
      theme = mythemeXRot,
      outlierShape=outlierShape
    ) + myfillscale
  ),
  tar_target(
    stdsVarHistoCer3Obj,
    stdsVarHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "3",
      theme = mythemeXRot
    ) + myfillscale
  ),
  tar_target(
    stdsRecoveryPercentHistoCer3Obj,
    stdsRecoveryPercentHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "3",
      theme = mythemeXRot
    ) + myfillscale
  ),
  tar_target(
    boxplotCer4Obj,
    stdsBoxplot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "4",
      theme = mythemeXRot,
      outlierShape=outlierShape
    ) + myfillscale
  ),
  tar_target(
    stdsVarHistoCer4Obj,
    stdsVarHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "4",
      theme = mythemeXRot
    ) + myfillscale
  ),
  tar_target(
    stdsRecoveryPercentHistoCer4Obj,
    stdsRecoveryPercentHistoPlot(
      averagedConcentrations = averagedConcentrations,
      selectCeramideId = "4",
      theme = mythemeXRot
    ) + myfillscale
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
    readr::write_csv(qcAveragedConcentrations, file="output/qcCal1Cal2AveragedConcentrations.csv")
  ),
  tar_target(
    qcAveragedConcentrationsLong,
    createQcAveragedConcentrationsLong(qcAveragedConcentrations)
  ),
  tar_target(
    qcConcentrationsPlotObj,
    qcConcentrationsPlot(qcAveragedConcentrationsLong, theme=mytheme, fillscale=myfillscale, outlierShape=19)
  ),
  tar_target(
    qcConcentrationsByLabIdPlotObj,
    qcConcentrationsByLabIdPlot(qcAveragedConcentrations, theme=mytheme, fillscale=myfillscale, outlierShape=19)
  ),
  tar_target(
    qcConcentrationsByLabIdDensityPlotObj,
    qcConcentrationsByLabIdDensityPlot(qcAveragedConcentrations, theme=mytheme, fillscale=myfillscale, outlierShape=19)
  ),
  tar_target(
    qcSampleStats,
    createQcConcentrationsStats(qcAveragedConcentrations)
  ),
  tar_target(
    qcSampleStatsPlotObj,
    qcSampleStatsPlot(qcSampleStats, theme=mytheme, fillscale=myfillscale)
  ),
  ################################################################################
  # NIST Samples concentration calculation
  ################################################################################
  tar_target(
    IntraAssayNISTwide,
    intraAssayNIST |>
      left_join(
        expectedStdsConcentrations |> filter(Sample == "STD 1") |> select(-Sample),
        by = c("ceramideId", "ceramideName")
      ) |>
      mutate(Unit = unique(
        expectedCalibrationLineConcentrations$Unit
      )) |> # add unit
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
      group_by(LabId, SampleType, ceramideId, ceramideName, Unit) |> # here, Samples are independent analytical replicates, while replicates are repeated injections of the same sample
      filter(n() >= 3 &
               (!is.na(area_l) |
                  !is.na(area_h))) # remove all NA entries and groups with less than three entries
  ),
  tar_target(
    combinedIntraAssayNISTWithLms,
    IntraAssayNISTwide |> ungroup() |> group_by(LabId, ceramideId, ceramideName) |>
      left_join(
        calibLineDataLmCoeffsWide |>
          ungroup() |>
          mutate(CalibrationLine = SampleType) |>
          select(-SampleType) |>
          group_by(LabId, ceramideId, ceramideName),
        by = c("LabId", "ceramideId", "ceramideName")
      )
  ),
  tar_target(
    manual_nist_test_C_SinglePoint,
    combinedIntraAssayNISTWithLms[1, ]$area_l / combinedIntraAssayNISTWithLms[1, ]$area_h *
      combinedIntraAssayNISTWithLms[1, ]$ISConcentration
  ),
  tar_target(
    nistAnalyteConcentrationsFromCalibLines,
    combinedIntraAssayNISTWithLms |>
      mutate(
        C_Adj = C_A_cal(S_A = area_l / area_h, a = SlopeX, b = Intercept),
        C_SinglePoint = area_l / area_h * ISConcentration
      )
  ),
  tar_target(
    stopAtSinglePoint,
    stopifnot(
      manual_nist_test_C_SinglePoint == nistAnalyteConcentrationsFromCalibLines[1, ]$C_SinglePoint
    )
  ),
  tar_target(
    nistAreaRatioPlotsObj,
    nistAreaRatioPlots(
      data = nistAnalyteConcentrationsFromCalibLines,
      sampleTypes = unique(nistAnalyteConcentrationsFromCalibLines$SampleType),
      theme = mythemeXRot,
      fillscale = myfillscale,
      outlierShape=outlierShape
    )
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
      cunit = cunit
    )
  ),
  tar_target(
    nistAveragedConcentrations,
    nistAnalyteConcentrationsFromCalibLines |>
      ungroup() |>
      select(
        LabId,
        SampleType,
        Sample,
        ceramideId,
        ceramideName,
        Unit,
        replicate,
        Protocol,
        Instrument,
        MassAnalyzerType,
        CalibrationLine,
        C_Adj,
        C_SinglePoint
      ) |>
      group_by(
        LabId,
        SampleType,
        Sample,
        ceramideId,
        ceramideName,
        Unit,
        replicate
      ) |>
      pivot_wider(names_from = CalibrationLine, values_from = C_Adj) |>
      mutate(Avg_C_Adj = (`Calibration Line 1` + `Calibration Line 2`) /
               2,
             # Avg_C_SinglePoint=(`C_SinglePoint_Calibration Line 1` + `C_SinglePoint_Calibration Line 2`)/2,)
      )
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
        mutate(C_CalibCurve_Avg = Avg_C_Adj) |>
        pivot_longer(
          cols = c("C_CalibCurve_Avg", "C_SinglePoint"),
          names_to = "Calibration",
          values_to = "Adj_Conc"
        ) |> left_join(massAnalyzerSummary, by = c("MassAnalyzerType"))
    ),
    tar_target(
      nistAveragedConcentrationsLongFile,
      write_csv(nistAveragedConcentrationsLong, file = "output/nistCal1Cal2AveragedConcentrationsLong.csv")
    ),
    tar_target(
      nistAveragedConcentrationsFile,
      write_csv(nistAveragedConcentrations, file = "output/nistCal1Cal2AveragedConcentrations.csv")
    ),
    ################################################################################
    # Comparison to other ring trial results, Bowden, Quehenberger, Giera, Biocrates
    ################################################################################
    tar_target(
      ringTrialsComparisonTbl,
      loadAndCombineRingTrialData(nistAveragedConcentrations)
    ),
    tar_target(
      ringTrialsComparisonPlotObj,
      ringTrialsComparisonPlot(
        ringTrialsComparisonTbl,
        fillscale=myfillscale,
        theme=mythemeXRot,
        outlierShape=outlierShape
      )
    ),
    ################################################################################
    # NIST Sample and Ceramide Concentrations Plots
    ################################################################################
    tar_target(
      nistConcentrationsPlotObj,
      nistConcentrationsPlot(nistAveragedConcentrationsLong, 
                             fillscale=myfillscale,
                             theme=mythemeXRot,
                             outlierShape=outlierShape)
    ),
    ################################################################################
    # NIST Sample and Ceramide Cv Plots
    ################################################################################
    tar_target(
      nistCvPlotObj,
      nistCvPlot(nistConcentrationsStatsLong, 
                 fillscale=myfillscale,
                 theme=mythemeXRot,
                 outlierShape=outlierShape)
    ),
    ################################################################################
    # NIST Concentrations by Sample type, LabId and Mass Analyzer Plots
    ################################################################################
    tar_target(
      nistAveragedConcentrationsLongPlotObjs,
      nistAveragedConcentrationsLong |>
        ungroup() |>
        group_by(SampleType, Unit) |>
        group_walk(
          ~ nistConcentrationsPlotBySampleType(
            data = .x,
            namesuffix = unique(.y$SampleType),
            unit = unique(nistAveragedConcentrations$Unit),
            theme = mythemeXRot,
            xfillscale = myfillscale,
            outlierShape=outlierShape
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
            outlierShape=outlierShape
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
        outlierShape=outlierShape
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
        outlierShape=outlierShape
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
        outlierShape=outlierShape
      )
    )  
)
  
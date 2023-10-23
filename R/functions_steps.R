createAssayTable <- function(reportsToInclude, ceramideColNames, filePrefix="Report Template Standard #", protocol = "Standard", reportsDir = ".", na, blankTypes, outputDirectory = "output") {
  assayTable <-
    loadLabDataAndBindRows(
      labIds = reportsToInclude, 
      ceramideColNames = ceramideColNames,
      filePrefix = filePrefix,
      protocol = protocol,
      reportsDir = reportsDir,
      na = na,
      blankTypes = blankTypes
    ) |> 
    left_join(instruments, by="LabId") |> 
    mutate(Instrument=replace_na(Instrument, "Unknown"), MassAnalyzerType=replace_na(MassAnalyzerType, "Unknown"))
  readr::write_csv(assayTable, file.path(outputDirectory, paste0("all-labs-", tolower(protocol), "-unfiltered.csv")))
  assayTable
}

prepareCalibLinesWithStandardConcs <- function(assayTable, expectedStdsConcentrations, outputDirectory = "output") {
  print(paste("Calculation Calibration Lines with Expected STD concentrations"))
  d_part1 <- assayTable |>
    filter(grepl("Calibration Line", SampleType)) |>
    left_join(expectedStdsConcentrations, by=c("ceramideId", "ceramideName", "Sample")) |>
    group_by(LabId, SampleType, Sample, ceramideId, ceramideName, Unit, replicate) |>
    mutate(RatioLipidToIS=area[isotope=="l"]/area[isotope=="h"],
           PAR=1/RatioLipidToIS,
           RatioIsConcToStdConc=ISConcentration/Concentration
    ) |>
    pivot_wider(names_from=isotope, values_from=area, names_prefix="area_") |>
    ungroup() |>
    group_by(LabId, SampleType, Sample, ceramideId, ceramideName, Unit) |>
    filter(n()>=1 & !is.na("area_l") & !is.na("area_h")) |> 
    mutate(
      sdev_area_ratio=sd(area_l/area_h)
    ) |>
    ungroup()

  #readr::write_csv(d_part1, file.path(outputDirectory, "debug_1.csv"))
  d_part2 <- assayTable |>
    filter(grepl("Calibration Line", SampleType)) |>
    left_join(expectedStdsConcentrations, by=c("ceramideId", "ceramideName", "Sample")) |>
    group_by(
      LabId,
      SampleType,
      Sample,
      Unit,
      replicate
    ) |>
    mutate(
      RatioLipidToIS_C16 = area / if_else(length(ceramideName[isotope == "h" & ceramideName == "Cer 18:1;O2/16:0"]) > 0, (area[isotope == "h" & ceramideName == "Cer 18:1;O2/16:0"])[1], NA_real_),
      RatioLipidToIS_C24 = area / if_else(length(ceramideName[isotope == "h" & ceramideName == "Cer 18:1;O2/24:0"]) > 0, (area[isotope == "h" & ceramideName == "Cer 18:1;O2/24:0"])[1], NA_real_)
    ) |> 
    filter(isotope == "l") |>  
     pivot_wider(names_from=isotope, values_from=area, names_prefix="area_") |> 
     ungroup()
  
  #readr::write_csv(d_part2, file.path(outputDirectory, "debug_2.csv"))
      
  d_final <- d_part1 |> 
    inner_join(d_part2 |> dplyr::select(Sample, SampleType, LabId, ceramideName, replicate, RatioLipidToIS_C16, RatioLipidToIS_C24), by = c("Sample", "SampleType", "LabId","ceramideName", "replicate" ))
  #readr::write_csv(d_final, file.path(outputDirectory, "debug.csv"))
  d_final
}

createQcAveragedConcentrations <- function(intraAssayQC, expectedStdsConcentrations, expectedCalibrationLineConcentrations, calibLineDataLmCoeffsWide) {
  # LOD from calibration lines
  # R_LoD_C Ratio limit of detection (y-axis)
  
  # LOD from total blanks and matrix blanks
  # LoD_TB, LoD_MB
  
  # dilution steps, MQC has an unknown concentration of unlabeled ceramides
  # but a known concentration of labeled ceramides -> use ratio to calculate the unlabeled concentration
  # MQC: 1, LQC: 0.5, LLQC: 0.33
  # HQC: MQC + 10muL of STD3 (unlabeled and labeled ceramides)
  # HLQC: MQC + 20muL of STD3 (unlabeled and labeled ceramides)
  # 570 muL EtAc:isoprop + 10muL pooled QC (NIST) + 20muL labelled IS mixture 
  # -> 600 muL total volume
  # 550 muL calibr. line samples + 10 muL 5% BSA + 20muL labelled IS mixture + 20muL non-labelled mix dilutions (STD1-6) 
  # -> 600 muL total volume
  
  # Accuracy: Percentage of nominal concentrations -> C_adj/C_nom*100
  
  # TODO: Blank handling
  # IntraAssayBlankQCwide <- intraAssayQC |> 
  #   filter(SampleType=="Blank QC") |>
  #   mutate(Unit=unique(expectedCalibrationLineConcentrations$Unit)) |> # add unit
  #   group_by(LabId, SampleType, Sample, ceramideId, ceramideName, Unit, replicate) |>
  #   mutate(RatioLipidToIS=area[isotope=="l"]/area[isotope=="h"],
  #          PAR=1/RatioLipidToIS
  #   ) |>
  #   pivot_wider(names_from=isotope, values_from=area, names_prefix="area_") |> ungroup() |>
  #   group_by(LabId, SampleType, ceramideId, ceramideName, Unit) |> # here, Samples are independent analytical replicates, while replicates are repeated injections of the same sample
  #   filter(n()>=3 & (!is.na(area_l) | !is.na(area_h))) # remove all NA entries and groups with less than three entries
  
  IntraAssayQCwide <- intraAssayQC |> 
    filter(SampleType!="Blank QC") |>
    left_join(expectedStdsConcentrations |> filter(Sample=="STD 1") |> select(-Sample), by=c("ceramideId","ceramideName")) |>
    mutate(Unit=unique(expectedCalibrationLineConcentrations$Unit)) |> # add unit
    group_by(LabId, SampleType, Sample, ceramideId, ceramideName, Unit, replicate) |>
    mutate(RatioLipidToIS=area[isotope=="l"]/area[isotope=="h"],
           PAR=1/RatioLipidToIS
    ) |>
    pivot_wider(names_from=isotope, values_from=area, names_prefix="area_") |> ungroup() |>
    group_by(LabId, SampleType, ceramideId, ceramideName, Unit) |> # here, Samples are independent analytical replicates, while replicates are repeated injections of the same sample
    filter(n()>=3 & (!is.na(area_l) | !is.na(area_h))) # remove all NA entries and groups with less than three entries
  
  combinedIntraAssayQCWithLms <-
    IntraAssayQCwide |> ungroup() |> group_by(LabId, ceramideId, ceramideName) |>
    left_join(
      calibLineDataLmCoeffsWide |> 
        ungroup() |> 
        mutate(CalibrationLine=SampleType) |> 
        select(-SampleType) |> 
        group_by(LabId, ceramideId, ceramideName),
      by = c("LabId", "ceramideId", "ceramideName"),
      relationship = "many-to-many"
    ) 
  
  qcAnalyteConcentrationsFromCalibLines <- combinedIntraAssayQCWithLms |>
    mutate(
      C_Adj=C_A_cal(S_A = 2 * area_l/area_h, a = SlopeX, b = Intercept),
      C_SinglePoint=2 * area_l/area_h*ISConcentration,
      SampleType = forcats::as_factor(as.character(SampleType))
    )
  
  qcAveragedConcentrations <- qcAnalyteConcentrationsFromCalibLines |> 
    ungroup() |> 
    select(LabId, SampleType, Sample, ceramideId, ceramideName, Unit, Protocol, Instrument, MassAnalyzerType, replicate, CalibrationLine, C_Adj, C_SinglePoint) |> 
    group_by(LabId, Sample, ceramideId, ceramideName, Unit, replicate) |> 
    pivot_wider(names_from=CalibrationLine, values_from=C_Adj) |> 
    mutate(
      Avg_C_Adj=(`Calibration Line 1` + `Calibration Line 2`)/2,
      C_SinglePoint=C_SinglePoint,
    )
  qcAveragedConcentrations
}

createQcAveragedConcentrationsLong <- function(qcAveragedConcentrations) {
  massAnalyzerSummary <- qcAveragedConcentrations |> 
    ungroup() |> 
    group_by(LabId) |> 
    summarise(MassAnalyzerType=first(MassAnalyzerType)) |>
    group_by(MassAnalyzerType) |> 
    summarise(n_massAnalyzer=n()) |> 
    mutate(facetLabel=paste0(MassAnalyzerType, " n=", n_massAnalyzer))
  
  qcAveragedConcentrationsLong <- qcAveragedConcentrations |> 
    mutate(C_CalibCurve_Avg=Avg_C_Adj) |> 
    pivot_longer(
      cols=c("Avg_C_Adj", "C_SinglePoint"),
      names_to = "Calibration",
      values_to = "Adj_Conc"
    ) |> left_join(massAnalyzerSummary, by = c("MassAnalyzerType"))
  qcAveragedConcentrationsLong
}

createQcConcentrationsStats <- function(qcAveragedConcentrations, outputDirectory = "output") {
  qcConcentrationsStats <- qcAveragedConcentrations |> 
    group_by(LabId, SampleType, ceramideId, ceramideName, Unit, Protocol, Instrument, MassAnalyzerType) |>
    summarise(
      AvgAvg_C_Adj = mean(Avg_C_Adj), 
      SdAvg_C_Adj = sd(Avg_C_Adj),
      CV_C_Adj = SdAvg_C_Adj/AvgAvg_C_Adj,
      Avg_C_SinglePoint = mean(C_SinglePoint),
      SdAvg_C_SinglePoint = sd(C_SinglePoint),
      CV_C_SinglePoint = SdAvg_C_SinglePoint/Avg_C_SinglePoint
    )
  readr::write_csv(file=file.path(outputDirectory, "qcConcentrationsStats.csv"), qcConcentrationsStats)
  qcConcentrationsStats
}


analyteConcentrationsFromCalibLines_function <- function(Auth, C16, C24) {
  Auth |>
    left_join(C16 |> select(LabId, ceramideName, replicate, Sample, SampleType, C_Adj_C16)) |> 
    left_join(C24 |> select(LabId, ceramideName, replicate, Sample, SampleType, C_Adj_C24))
}


get_nistAveragedConcentrations <- function(data){
  data |>
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
      C_Adj_C16,
      C_Adj_C24,
      C_SinglePoint,
      RatioLipidToIS
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
    pivot_wider(
      names_from = CalibrationLine,
      values_from = c(C_Adj, C_Adj_C16, C_Adj_C24, C_SinglePoint, RatioLipidToIS)
    ) |>
    mutate(
      Avg_C_Adj = rowMeans(cbind(`C_Adj_Calibration Line 1`, `C_Adj_Calibration Line 2`), na.rm=TRUE),
      Avg_C_Adj_C16 = rowMeans(cbind(`C_Adj_C16_Calibration Line 1`, `C_Adj_C16_Calibration Line 2`), na.rm=TRUE),
      Avg_C_Adj_C24 =  rowMeans(cbind(`C_Adj_C24_Calibration Line 1`, `C_Adj_C24_Calibration Line 2`), na.rm=TRUE),
      C_SinglePoint =  rowMeans(cbind(`C_SinglePoint_Calibration Line 1`, `C_SinglePoint_Calibration Line 2`), na.rm=TRUE),
      RatioLipidToIS = rowMeans(cbind(`RatioLipidToIS_Calibration Line 1`, `RatioLipidToIS_Calibration Line 2`), na.rm=TRUE))
    
}

readLCMSDetailsSpecs <- function(labId, filename, na = c("","N/F","NA","N/A", "NF", "<LOD", "not inj", "N.D.", "-", "ND", "failed inj", "n/a")) {
  readxl::read_excel(
    filename, 
    sheet = "LC-MRM - Overview", 
    range = "B11:C20",
    col_names = c("name", "value"),
    col_types = c("text", "text"),
    na = na
  ) |> mutate(LabId=labId) |> pivot_wider()
}

readSubTable <- function(labId, filename, sheet, range, columnNames, sampleType, ceramideIdLookup=cerIdLookup, na = c("","N/F","NA","N/A", "NF", "<LOD", "not inj", "N.D.", "-", "ND", "failed inj", "n/a")) {
  readxl::read_excel(
    filename, 
    sheet = sheet, 
    range = range,
    col_names = columnNames,
    col_types = c("text", rep("numeric", 4*3*2)),
    na = na
  ) |> pivot_longer(
    cols = starts_with("Cer"), 
    names_to = c("ceramideId", "isotope", "replicate"),
    names_pattern = "Cer(.*)-(.*):R(.*)",
    values_to = "area"
  ) |> mutate(
    ceramideId = as.character(ceramideId),
    SampleType = sampleType, 
    LabId = labId
  ) |> left_join(
    ceramideIdLookup, by = c("ceramideId" = "ceramideId")
  )
}

readDatasetSummary <- function(filename, sheet, range, columnNames, na=c("")) {
  readxl::read_excel(
    filename, 
    sheet = sheet, 
    range = range,
    col_names = columnNames,
    col_types = "text",
    na = na
  ) |> mutate(LabId=as.character(LabId))
}

loadLabDataLcMRMIntraAssayQc <- function(labId, filename, ceramideColNames, na, blankTypesColumnRange) {
  intraAssayCalLine1 <- readSubTable(
    labId = labId,
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B7:Z12",
    columnNames = ceramideColNames,
    sampleType = "Calibration Line 1",
    na = na
  )
  
  intraAssayCalLine2 <- readSubTable(
    labId = labId,
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B18:Z23",
    columnNames = ceramideColNames,
    sampleType = "Calibration Line 2",
    na = na
  )
  
  intraAssayNIST_SRM <- readSubTable(
    labId = labId,  
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B29:Z34",
    columnNames = ceramideColNames,
    sampleType = "NIST SRM",
    na = na
  )
  
  intraAssayNIST_hTAG <- readSubTable(
    labId = labId,  
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B38:Z43",
    columnNames = ceramideColNames,
    sampleType = "NIST hTAG",
    na = na
  )
  
  intraAssayNIST_T1D <- readSubTable(
    labId = labId,  
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B47:Z52",
    columnNames = ceramideColNames,
    sampleType = "NIST T1D",
    na = na
  )
  
  intraAssayNIST_YAA <- readSubTable(
    labId = labId,  
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B56:Z61",
    columnNames = ceramideColNames,
    sampleType = "NIST YAA",
    na = na
  )
  
  intraAssayLowestLevelQC <- readSubTable(
    labId = labId,  
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B67:Z72",
    columnNames = ceramideColNames,
    sampleType = "Lowest Level QC",
    na = na
  )
  
  intraAssayLowQC <- readSubTable(
    labId = labId,  
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B76:Z81",
    columnNames = ceramideColNames,
    sampleType = "Low QC",
    na = na
  )
  
  intraAssayMiddleQC <- readSubTable(
    labId = labId,  
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B85:Z90",
    columnNames = ceramideColNames,
    sampleType = "Middle QC",
    na = na
  )
  
  intraAssayHighQC <- readSubTable(
    labId = labId,  
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B94:Z99",
    columnNames = ceramideColNames,
    sampleType = "High QC",
    na = na
  )
  
  intraAssayHighestLevelQC <- readSubTable(
    labId = labId,  
    filename = filename,
    sheet = "LC-MRM - Intra Assay QC Res",
    range = "B103:Z108",
    columnNames = ceramideColNames,
    sampleType = "Highest Level QC",
    na = na
  )
  
  intraAssayBlankQC <- data.frame()
  # select only matrix QCs here
  if (length(blankTypesColumnRange)>0) {
    intraAssayBlankQC <- readSubTable(
      labId = labId,  
      filename = filename,
      sheet = "LC-MRM - Intra Assay QC Res",
      range = blankTypesColumnRange,
      columnNames = ceramideColNames,
      sampleType = "Blank QC",
      na = na
    )
  }
  
  intraAssayTable <- bind_rows(
    intraAssayCalLine1,
    intraAssayCalLine2,
    intraAssayNIST_SRM,
    intraAssayNIST_hTAG,
    intraAssayNIST_T1D,
    intraAssayNIST_YAA,
    intraAssayLowestLevelQC,
    intraAssayLowQC,
    intraAssayMiddleQC,
    intraAssayHighQC,
    intraAssayHighestLevelQC,
    intraAssayBlankQC
  )
  intraAssayTable
}

loadLCMSDetails <- function(labIds, filePrefix="Report Template Standard #", fileSuffix=".xlsx", reportsDir=".", na) {
  labIds |> 
    map_dfr(
      ~ readLCMSDetailsSpecs(
        .x, 
        file.path(reportsDir, paste0(filePrefix, .x, fileSuffix)),
        na
      )
    )
}

loadLabDataAndBindRows <- function(labIds, ceramideColNames, filePrefix="Report Template Standard #", fileSuffix=".xlsx", protocol="Standard", reportsDir=".", na, blankTypes) {
  blankTypesColumnRange <- blankTypes[blankTypes$LabId %in% labIds, ]
  df <- labIds |> 
    map_dfr(
      ~ loadLabDataLcMRMIntraAssayQc(
        .x, 
        file.path(reportsDir, paste0(filePrefix, .x, fileSuffix)), 
        ceramideColNames,
        na,
        blankTypesColumnRange[blankTypesColumnRange$LabId==.x, ]$Range
      )
    ) %>% mutate(
      LabId=factor(.data$LabId, levels = labIds), 
      Protocol=protocol
    )
}

loadAndCombineRingTrialData <- function(nistAveragedConcentrations, reportsDir="data") {
  # d_ils_orig <- read_csv(here:here("data/Cer-ILS-Ringtrial_FullDataset_SLING_V1.csv"))
  d_bowden_orig <- read_csv(file.path(reportsDir, "ring-trial-comparison","Bowden_InterlabStudy_Ceramides_Conc.csv"))
  d_quehenberger_orig <- read_csv(file.path(reportsDir, "ring-trial-comparison","Quehenberger_Ceramides_Conc.csv"))
  d_biocrates <- read_csv(file.path(reportsDir, "ring-trial-comparison","Biocrates_Ringtrial_Ceramides.csv"))
  d_sum_giera_orig <- read_csv(file.path(reportsDir, "ring-trial-comparison","Giera_MEDM_COD_V1.csv"))
  
  d_ils <- nistAveragedConcentrations |> 
    dplyr::select(LabID = LabId, SampleName = SampleType, Replicate = replicate, 
                  Compound = ceramideName, Protocol, MassAnalyzer = MassAnalyzerType, 
                  Conc_SP = C_SinglePoint, Conc_MP = Avg_C_Adj)
  
  d_ils_labs <- d_ils |> 
    group_by(LabID, SampleName,Compound) |> 
    summarise(
      Conc_SP_mean = mean(Conc_SP, na.rm = TRUE), 
      Conc_MP_mean = mean(Conc_MP, na.rm = TRUE),
      CV_intra_SP = sd(Conc_SP, na.rm = TRUE)/mean(Conc_SP, na.rm = TRUE) *100,
      CV_intra_MP = sd(Conc_MP, na.rm = TRUE)/mean(Conc_MP, na.rm = TRUE) * 100
    ) |> 
    mutate(Study = "ILS Ceramides (MP)", .before = 1) |> 
    ungroup()
  
  d_biocrates_mean <- d_biocrates |>
    mutate(id = row_number(), .before =1) |> 
    pivot_longer(-id:-LabId, names_to = "Compound", values_to = "Conc") |> 
    group_by(LabId, Compound) |>
    summarise(
      Conc_mean = mean(Conc, na.rm = TRUE),
      CV_intra = sd(Conc, na.rm = TRUE)/mean(Conc, na.rm = TRUE)*100) |> 
    ungroup() |> 
    mutate(LabId = as.character(LabId),
           Study = "Moseley 2019 (SP)", .before = 1) |> 
    dplyr::select(Study, Compound, LabID=LabId, Conc_mean, CV_intra)
  
  d_biocrates_mean[d_biocrates_mean == "NaN"] <- NA
  
  d_bowden <- d_bowden_orig |> 
    select(LabID = LabId, Compound = Species, Conc_mean = Conc) |> 
    mutate(LabID = as.character(LabID)) |> 
    mutate(Study = "Bowden 2017 (MP)")
  
  d_quehenberger <- d_quehenberger_orig |> 
    select(LabID = LabId, Compound = Species, Conc_mean = Conc) |> 
    mutate(LabID = as.character(LabID)) |> 
    mutate(Study = "Quehenberger 2010")
  
  d_sum_giera <- d_sum_giera_orig |>
    select(Study, Compound = Species, everything()) |>
    mutate(Study = "Giera 2021 (SP)")
  
  d_plot_1 <- d_ils_labs |> 
    filter(SampleName == "NIST SRM") |> 
    select(Study, LabID, Compound, Conc_mean = Conc_SP_mean, CV_intra = CV_intra_SP) |> 
    bind_rows(d_biocrates_mean) |> 
    bind_rows(d_bowden) |> 
    bind_rows(d_quehenberger) |> 
    bind_rows(d_sum_giera) |>
    mutate(Study = factor(Study, 
                          levels = c("Quehenberger 2010", "Bowden 2017 (MP)", 
                                     "Moseley 2019 (SP)" ,"Giera 2021 (SP)", 
                                     "ILS Ceramides (MP)"))) |> 
    arrange(Study)
  
  d_plot_1
}

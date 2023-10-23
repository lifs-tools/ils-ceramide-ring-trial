readAssayTables <- function(standardReportsToInclude,
                            preferredReportsToInclude, 
                            ceramideColNames,
                            reportsDir,
                            na,
                            blankTypes) {
  print(paste("Creating Preferred Method Reports Table"))
  preferredIntraAssayTable <- createAssayTable(
    preferredReportsToInclude, 
    ceramideColNames,
    "Report Template Preferred #",
    protocol = "Preferred",
    reportsDir = reportsDir,
    na = naValues,
    blankTypes = blankTypes
  )
  
  print(paste("Creating Standard Method Reports Table"))
  intraAssayTable <- createAssayTable(
    standardReportsToInclude,
    ceramideColNames,
    filePrefix = "Report Template Standard #",
    protocol = "Standard",
    reportsDir = reportsDir,
    na = naValues,
    blankTypes = blankTypes
  )
  print(paste("Combining Method Reports Tables"))
  # combine standard and preferred protocols
  intraAssayTable <- bind_rows(intraAssayTable, preferredIntraAssayTable) |> as_tibble()
  
  filteredIntraAssayTable <- intraAssayTable |>
    filter(LabId != "14" |
             Sample != "STD 5" |
             replicate != 1) |> # filter outlier for LabId 14, STD 5, replicate 1 inc calibration lines 1 and 2
    filter(LabId != "14" |
             replicate != 3 |
             ceramideId != 3 |
             SampleType != "Calibration Line 1") |> # filter outlier for LabId 14, replicate 3, ceramide 3
    filter(LabId != "16" |
             Sample != "STD 6" |
             replicate != 1 |
             ceramideId != 1 |
             SampleType != "Calibration Line 1") |> # filter N/F / NAs for LabId 16, STD 6, repl. 1, Calibr. Line 1 and matching heavy
    filter(LabId != "16" |
             Sample != "STD 6" |
             replicate != 2 |
             ceramideId != 2 |
             SampleType != "Calibration Line 2") |> # filter N/F / NAs for LabId 16, STD 6, repl. 2, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 4" |
             #replicate != 1 |
             ceramideId != 1 |
             SampleType != "Calibration Line 1") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 4" |
             #replicate != 1 |
             ceramideId != 1 |
             SampleType != "Calibration Line 2") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 4" |
             #replicate != 1 |
             ceramideId != 2 |
             SampleType != "Calibration Line 1") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 4" |
             #replicate != 1 |
             ceramideId != 2 |
             SampleType != "Calibration Line 2") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 5" |
             #replicate != 1 |
             ceramideId != 1 |
             SampleType != "Calibration Line 1") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 5" |
             #replicate != 1 |
             ceramideId != 1 |
             SampleType != "Calibration Line 2") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 5" |
             #replicate != 1 |
             ceramideId != 2 |
             SampleType != "Calibration Line 1") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 5" |
             #replicate != 1 |
             ceramideId != 2 |
             SampleType != "Calibration Line 2") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 6" |
             #replicate != 1 |
             ceramideId != 1 |
             SampleType != "Calibration Line 1") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 6" |
             #replicate != 1 |
             ceramideId != 1 |
             SampleType != "Calibration Line 2") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 6" |
             #replicate != 1 |
             ceramideId != 2 |
             SampleType != "Calibration Line 1") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(LabId != "32" |
             Sample != "STD 6" |
             #replicate != 1 |
             ceramideId != 2 |
             SampleType != "Calibration Line 2") %>% #filter ND for LabId 32, STD 6, repl. 1, Calibr. Line 2 and matching heavy
    filter(!str_detect(SampleType, "^Calibration Line") | # filter NAs for calibration line data
             !is.na(area))
    return(filteredIntraAssayTable)
}

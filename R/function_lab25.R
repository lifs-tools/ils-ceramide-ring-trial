library(tidyverse)
filename <- "/home/nilshoffmann/Downloads/Report Template Standard #25-calibrationcurves.xlsx"
#xl <- readxl::read_xlsx(filename)

na <- c("","N/F","NA","N/A", "NF", "<LOD", "not inj", "N.D.", "-", "ND", "failed inj", "n/a")
labId1 <- "25a"
labId2 <- "25b"
ceramideColNames <- c("Sample",
                      "Cer1-l:R1", "Cer1-l:R2", "Cer1-l:R3", 
                      "Cer2-l:R1", "Cer2-l:R2", "Cer2-l:R3",
                      "Cer3-l:R1", "Cer3-l:R2", "Cer3-l:R3",
                      "Cer4-l:R1", "Cer4-l:R2", "Cer4-l:R3",
                      "Cer1-h:R1", "Cer1-h:R2", "Cer1-h:R3", 
                      "Cer2-h:R1", "Cer2-h:R2", "Cer2-h:R3",
                      "Cer3-h:R1", "Cer3-h:R2", "Cer3-h:R3",
                      "Cer4-h:R1", "Cer4-h:R2", "Cer4-h:R3"
)

intraAssayCalLine1a <- readSubTable(
  labId = labId1,
  filename = filename,
  sheet = "Tabelle1",
  range = "A5:Y10",
  columnNames = ceramideColNames,
  sampleType = "Calibration Line 1",
  na = na
)

intraAssayCalLine2a <- readSubTable(
  labId = labId1,
  filename = filename,
  sheet = "Tabelle1",
  range = "A16:Y21",
  columnNames = ceramideColNames,
  sampleType = "Calibration Line 2",
  na = na
)

intraAssayCalLine1b <- readSubTable(
  labId = labId2,
  filename = filename,
  sheet = "Tabelle1",
  range = "A27:Y32",
  columnNames = ceramideColNames,
  sampleType = "Calibration Line 1",
  na = na
)

intraAssayCalLine2b <- readSubTable(
  labId = labId2,
  filename = filename,
  sheet = "Tabelle1",
  range = "A38:Y43",
  columnNames = ceramideColNames,
  sampleType = "Calibration Line 2",
  na = na
)

calibrationLines <- bind_rows(list(intraAssayCalLine1a, intraAssayCalLine1b, intraAssayCalLine2a, intraAssayCalLine2b))

library(ggplot2)
ggplot(data=calibrationLines, mapping=aes(x=Sample, y=area, shape=isotope, color=LabId)) + geom_point() + facet_grid(ceramideName+isotope~SampleType, scales = "free_y")

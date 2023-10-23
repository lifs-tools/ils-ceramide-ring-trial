require("tidyverse")
require("ggsci")

if(Sys.getenv("DATA_DIR")!="") {
  dataDir <- Sys.getenv("DATA_DIR")
}else{
  dataDir <- "data"
}
reportsDir <- file.path(dataDir, "reports")
outputDir <- "output"


datasetSummaryColumnNames <- c(
  "LabId",
  "Type",
  "Processing Status",
  "Calibration Line 1",
  "Calibration Line 2",
  "NIST SRM",
  "NIST hTAG",
  "NIST T1D",
  "NIST Young AA",
  "LLQC",
  "LQC",
  "MQC",
  "HQC",
  "HLQC",
  "Blank"
)

datasetSummaryRange <- "A2:O45"

nPreferredReports <- 10
nStandardReports <- 21

naValues <- c("","N/F","NA","N/A", "NF", "<LOD", "not inj", "N.D.", "-", "ND", "failed inj", "n/a")

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

outlierShape <- NA

cerIdLookup <- data.frame(
  ceramideId = c("1","2","3","4"),
  ceramideName = c("Cer 18:1;O2/16:0","Cer 18:1;O2/18:0","Cer 18:1;O2/24:0","Cer 18:1;O2/24:1")
) |> mutate(ceramideId = as.character(ceramideId))

calibrationLineSamplesToUse <- 
 c("STD 1", "STD 2", "STD 3", "STD 4", "STD 5", "STD 6")

standardReportsToInclude <- c(
  #"02a", # removed for final evaluation, replaced by #02b which has now triplicate injections
  "02b",
  "04",
  #"06", # replaced by 18a, which is more complete 
  #"07",
  #"08", duplicate of 14
  "10a",
  #"11", duplicate of 14
  "12",
  "14",
  # massive outlier in STD5 replicate 1 for all Ceramides -> remove
  "15",
  "16",
  # 1 matrix blank, 2 total blanks, N/F
  "17a",
  # missing replicates 2 and 3 for NIST samples, matrix blank, total blank, no QC samples
  "18a",
  # matrix blank, total blank, NF
  "18b",
  # matrix blank, total blank
  "19",
  # missing replicates 2 and 3 for NIST samples (NA), 3 x total blank, 3 x matrix blank, NA
  "20a", 
  # no ISTD signals for sample yA-A 4, set all to empty
  "21", ### REMOVED BEFORE  only three qc samples instead of 6 with 3 replicates
  # missing replicates 2 and 3 for NIST samples, no blanks, empty for missing values
  "22",
  # renamed sheet "Main LC-MRM-Intra Assay QC Res" to "LC-MRM - Intra Assay QC Res", three blank rows, N/A
  "23",
  # missing replicates 2 and 3 for NIST samples, three blank rows, missing replicates 2 and 3
  "24", 
  # matrix blank, total blank, not inj / <LOD for missing values, missing replicates 2 and 3, two-times calibration lines 1 and 2
  #"25", 
   ### REMOVED BEFORE edited by switching STD 4 and 5 rows for calibration line 2
  "27",
  "28",
  "29a",
  "30",
  "31",
  "32",
  "34", # replicates 2 and 3 are missing
  "36",
  "38"
)

preferredReportsToInclude <- c(
  # "01", sheet with quantities is missing
  "03",
  "05", # seems to have an issue with calibration line data being reversed
  "07",
  "09", # FIA-MRM, table layout is different, HQC and HLQC seem to be too high
  "10b",
  "13", ### REMOVED BEFORE only raw areas
  "16",
  "17a",
  # missing replicates 2 and 3 for NIST samples, matrix blank, total blank, no QC samples
  #"17b", table layout is totally different, sheets missing, not all ceramides included, onl 18:1/16:0, 18:1/18:0 and 18:1/24:0 as labeled
  "20b", 
  ### REMOVED BEFORE has 4 replicates instead of three, removed smallest one, only 1 STD curve, issue with table? Zero value, high var of replicates. Bo: FTT communicated with lab, purple sheet is a differernt non-SOP extraction. #20a is the SOP
  #"26a", ### QQQ -  REMOVED BEFORE  FIA and QQresults, table layout is different, after communination lab provided peak areas intesities, added to template sheet by Bo 
  #"26b", ### FIA -  REMOVED BEFORE  FIA and QQresults, table layout is different, after communination lab provided peak areas intesities. added to template sheet by Bo 
  "29b", # UHPSFC different sheet name,
  "32",
  "33", # additional STD dilutions, using only defaults
  "35",
  "37"
)

singlePointCalibrationSamples <- c()

CalibrationLineSamples <- c()

expectedHeavyStdsConcentrations <-
  read_csv(file.path(dataDir, "definitions", "labeledISConcentrations.csv"),
    col_types = cols(
      ceramideName = col_character(),
      ceramideId = col_character(),
      Unit = col_character()
    )
  ) |> rename(ISConcentration = Concentration)

expectedLightStdsConcentrations <-
  read_csv(file.path(dataDir, "definitions", "nonlabeledStandardsConcentrations.csv"),
    col_types = cols(
      ceramideName = col_character(),
      ceramideId = col_character(),
      Unit = col_character()
    )
  )

expectedCalibrationLineConcentrations <-
  read_csv(file.path(dataDir, "definitions", "calibrationLineConcentrations.csv"),
    col_types = cols(
      ceramideName = col_character(),
      ceramideId = col_character(),
      Unit = col_character()
    )
  ) |> pivot_longer(
    cols = c("STD 1", "STD 2", "STD 3", "STD 4", "STD 5", "STD 6"),
    names_to = "Sample",
    values_to = "Concentration"
  )

instruments <- 
  read_csv(file.path(dataDir, "definitions", "instruments.csv"),
    col_types = cols(
      LabId = col_character(),
      Instrument = col_character(),
      MassAnalyzerType = col_character()
    )
  )

blankTypesTable <- 
  read_csv(file.path(dataDir, "definitions", "blanktypes.csv"),
    col_types = cols(
      LabId = col_character(),
      BlankReportType = col_character()
    )
  )
thetheme <- theme_bw()

mythemeXRot <- thetheme + theme(axis.text.x = element_text(
  angle = 90,
  hjust = 1,
  vjust = 0.5
), legend.position="bottom")

mythemeXRot45 <- thetheme + theme(axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 1
), legend.position="bottom")

mytheme <- thetheme + theme(legend.position="bottom")

mycolorscale <- scale_colour_nejm()
myfillscale <- scale_fill_nejm()

################################################################################
# Nominal concentrations of labeled standards and ratios of light/heavy
################################################################################
expectedStdsConcentrations <- expectedCalibrationLineConcentrations |>
  left_join(
    expectedHeavyStdsConcentrations |> 
      select(-isotope)
  ) |>
  mutate(TheorRatioAnalyteIS=Concentration/ISConcentration)

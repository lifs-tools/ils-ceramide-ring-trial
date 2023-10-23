# Convert STD data in raw data table from dataset #13 to wide format


library(tidyverse)
library(readxl)
library(here)

d_orig <- read_xlsx(path = here("data/reports/Report Template Preferred #13.xlsx"),sheet = "Raw Peak Areas")

d <- d_orig |> 
  filter(str_detect(Sample_Name, "STD")) |> 
  arrange(Sample_Name) |> 
  separate(col = Sample_Name, into = c("SampleID","Replicate"), sep = "_") |> 
  pivot_wider(names_from = c("Replicate"), values_from = starts_with("Cer")) 

write_csv(d, here("data/reports/temp_dataset13_STD.csv"))


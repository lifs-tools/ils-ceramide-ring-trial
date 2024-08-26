library(targets)
library(gt)
library(parallel)
library(here)

# please see evaluation_debug.R for details on how to reset / debug a targets pipeline
# enable the following line to print details on warnings or on errors 
targets::tar_meta(fields = warnings, complete_only = TRUE)
#targets::tar_meta(fields = error, complete_only = TRUE)
#use a faster implementation for storing of intermediate tibbles
tar_option_set(
  memory = "transient", 
  garbage_collection = TRUE, 
  format = "fst_tbl", 
  storage = "worker", 
  retrieval = "worker"
)
# use the following command to run independent parts of the workflow in parallel
tar_make_future(workers = max(1, parallel::detectCores() - 1), garbage_collection = TRUE)
# run the following line to print details on warnings or on errors 
#targets::tar_meta(fields = warnings, complete_only = TRUE)
#targets::tar_meta(fields = error, complete_only = TRUE)

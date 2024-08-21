
rstudioapi::restartSession()

library(targets)
library(gt)
library(parallel)
library(here)


# use the following target to delete all intermediate results
#tar_destroy(ask = FALSE)

# use the following make target to allow browser() and other statements
#tar_make(callr_function = NULL)
# 
# enable the following line to print details on warnings or on errors 
#targets::tar_meta(fields = warnings, complete_only = TRUE)
#targets::tar_meta(fields = error, complete_only = TRUE)
tar_option_set(
  memory = "transient", 
  garbage_collection = TRUE, 
  format = "fst_tbl", #use a faster implementation for storing of intermediate tibbles
  storage = "worker", 
  retrieval = "worker"
)
# use the following command to run independent parts of the workflow in parallel
tar_make_future(workers = max(1, parallel::detectCores() - 1), garbage_collection = TRUE)
#tar_glimpse()
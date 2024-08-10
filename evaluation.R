
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
tar_option_set(memory = "transient", garbage_collection = TRUE)
# use the following command to run independent parts of the workflow in parallel
tar_make_future(workers = min(4, parallel::detectCores() - 1), garbage_collection = TRUE)
# 1#tar_make_future(workers = 1)
#tar_glimpse()
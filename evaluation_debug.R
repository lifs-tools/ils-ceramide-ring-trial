rstudioapi::restartSession()

library(targets)
library(gt)
library(parallel)
library(here)

# use the following target to delete all intermediate results
#tar_destroy(ask = FALSE)
# 
# use the following make target to allow browser() and other statements instead of using tar_make_future
# make sure to start it from a clean session and DO NOT USE the tar_option_set command above!
tar_make(callr_function = NULL)

# run the following line to print details on warnings or on errors 
#targets::tar_meta(fields = warnings, complete_only = TRUE)
#targets::tar_meta(fields = error, complete_only = TRUE)

library(targets)
library(gt)
library(parallel)
# To halt on errors with browser() use the following line to execute
# tar_make(callr_function = NULL)

# Run the workflow
tar_make_future(workers = min(1, parallel::detectCores() - 1))
tar_glimpse()

 2018-04-30 millard@insa-toulouse.fr
#
# Calculation of growth rate and extracellular (production and consumption) fluxes, under metabolic steady-states.
# See physio_fit.R for details.
#
# Copyright 2018, INRA, France
# License: GPL v2

# Reinitialize the session
rm(list=ls(all=TRUE))

# Go to PhysioFit directory
setwd("D:/GIT/PhysioFit")

# Load the main script and its dependencies
source("physio_fit.R")

# Go to the working directory
setwd("D:/GIT/PhysioFit/Example")

# Run calculation: a simple example here, please refers to PhysioFit documentation for details

# estimate first-order degradation constants
res_deg <- estimate_k("example_blk")

# estimate fluxes
res <- physio_fit(datfile="example", mc=TRUE, it=100, vini=0.04, weight = c("DHA"=1, "X"=0.001), deg=res_deg$deg_cst["DHA"])


# # Run in batch mode
# res_batch <- list()
# for (i in seq(4)){
#     res_batch[[i]] <- physio_fit(datfile=paste("example", i, sep=""), mc=TRUE, it=100)
# }
# lapply(res_batch, FUN = `[[`, "sens")









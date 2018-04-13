# 2016-08-18 millard@insa-toulouse.fr
#
# Calculation of growth rate and extracellular (production and consumption) fluxes, under metabolic steady-states.
# See physio_fit.R for details.
#
# Copyright 2016, INRA, France
# License: GPL v2

# Reinitialize the session
rm(list=ls(all=TRUE))

# Go to PhysioFit directory
setwd("D:/Users/millard/Desktop/Old_desktop/Drafts en cours/PhysioFit/PhysioFit")

# Load the main script and its dependencies
source("physio_fit.R")

# Go to the working directory
setwd("D:/Users/millard/Desktop/Old_desktop/Drafts en cours/FucRha")

# Run calculation: a simple example here, please refers to PhysioFit documentation for details
res <- physio_fit(datfile="3b", mc=TRUE, it=100, vini=0.2,weight=c("X"=0.002,"Fuc"=0.2,"Pdiol"=0.05,"Ace"=0.02,"Lac"=1e-3,"For"=0.01))

# Run in batch mode
res_batch <- list()
for (i in seq(4)){
    res_batch[[i]] <- physio_fit(datfile=paste("example", i, sep=""), mc=TRUE, it=100)
}
lapply(res_batch, FUN = `[[`, "sens")









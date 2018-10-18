####################
### PhysioFit    ###
####################
#
# version 1.0.1
#
# This R script allows to i) quantify exchange (production and consumption) fluxes and ii) cell growth
# rate during (batch) cultivations from time-course measurements of extracellular metabolites and biomass
# concentration, assuming that cells are in metabolic (pseudo) steady-state.
#
# Initial concentrations (i.e. concentrations of biomass - X - and metabolites - M - at t=0)  and fluxes
# (exchange fluxes - qM - and growth rate - mu -) are estimated by fitting to experimental data.
#
# Simulations, i.e. calculation of X(t) and M(t), are performed using analytical functions (see the code
# for details). Fluxes can be corrected for non-enzymatic degradation of some metabolites (assuming a
# first-order kinetics), and lag time in the cultivations can be considered.
#
# The sign of qM can be used to automatically identify products and substrates in high throughput workflows
# for functional analysis metabolic systems.
#
# Optimization is performed using the nlsic algorithm (see Sokol et al., Bioinformatics, 2012 for details) by
# minimizing the following cost function:
#
#     residuum = ((sim-meas)/weight)**2
#
# where 'sim' is the simulated data, 'meas' denotes measurements, and 'weight' is a weighting factor (e.g.,
# standard deviation on measurements).
#
# A global sensitivity analysis (Monte-Carlo approach) is implemented to evaluate the precision of the
# estimated parameters, and plots are generated to facilitate the evaluation of the fitting quality.
#
# Author: pierre.millard@insa-toulouse.fr
#         MetaSys team, LISBP, Toulouse, France
#
# Copyright 2016-2018, INRA, France
# License: GPL v2
#
#
####################
### Dependencies ###
####################
#
# The following R packages must be installed: nnls, numDeriv, RColorBrewer.
# They can be installed by running the following commands in R:
#
#    install.packages("nnls")
#    install.packages("numDeriv")
#    install.packages("RColorBrewer")
#
# The nlsic.R script is required for optimization and must be included in the physio_fit.R directory.
# 
#
####################
### Input file   ###
####################
#
# Experimental data must be provided in a tabulated text file ('.txt' extension), with the following columns:
#
#       time    experimental times
#       X       biomass concentration
#       meta_y  metabolite y concentration (1 column per metabolite)
#
# Units of estimated parameters depend on the units of experimental data, e.g. fluxes will be given in
# mmol/gDW/h if metabolite concentrations are in mM, times in h, and biomass concentration in gDW/L.
#
# Missing data must be given as NAs in the input file.
#
# Weighting factors (e.g. experimental standard deviation on measurements) used to calculate cost during
# optimization can be provided in a file with the same format (and named 'xxx_sd.txt', where xxx is the
# name of data file), except for the 'time' column which must be removed. Weighting factors may also be
# provided as an argument (see the Usage section).
#
# The files 'example.txt' and 'example_sd.txt' are provided as example.
#
#
####################
### Usage        ###
####################
# 
# Fit experimental data using the function physio_fit(file, ...), which takes the following arguments:
#
#       file    name of input file (without '.txt' extension)
#       mc      boolean, Monte-Carlo sensitivity analysis is performed if True (default = True)
#       it      number of Monte-Carlo iterations (default = 50)
#       pos     negative concentrations of noisy datasets generated during Monte-Carlo iterations are
#               set to 0 if True (default = True)
#       vini    initial value for fluxes and concentrations (default = 1)
#       upcf    upper constraints on fluxes (default = 50)
#       locf    lower constraints on fluxes (default = -50)
#       upcc    upper constraints on initial concentrations (default = 50)
#       locc    lower constraints on initial concentrations (default = 1e-6)
#       lag     estimate lag time (boolean, default = False)
#       deg     first order degradation constant for metabolites(s) which are non-enzymatically degraded
#               (named vector, default = NULL)
#       weight  weight matrix used for residuum calculation (as (sim-meas)/weight), can be:
#               - a matrix with the same format as the input file (without the 'time' column) and
#                 containing weights (e.g. SDs on measurements).
#               - a named vector containing weights for *all* the metabolites provided in the input file
#               - NULL (by default), in this case the matrix is automatically loaded from the file
#                 'xxx_sd.txt' - where xxx is the data file name - if this file exists, otherwise
#                 'weight' is constructed from 'sd_X' and 'sd_M' arguments.
#       sd_X    standard deviation on biomass concentration (default = 0.002),
#               used only if 'weight=NULL'
#       sd_M    standard deviation on metabolite concentrations (default = 0.5),
#               used only if 'weight=NULL'
#       save    boolean, save results (in txt & pdf files) if True (default = True)
#       summary boolean, results of the khi-2 test are displayed if True (default = False)
#
#  Returns a list with the following items (type):
#
#       $sys    $to_est         names of parameters to estimate (vector)
#               $params         parameters (vector)
#               $nconc          names of all concentrations (vector)
#               $nflux          names of fluxes (vector)
#               $metab          names of metabolites (vector)
#               $weight         standard deviation on measurements (matrix)
#               $te_upc         upper bounds of free parameters (vector)
#               $te_loc         lower bounds of free parameters (vector)
#               $u              linear constraint matrix, u*par>=co (matrix)
#               $co             constraints vector (vector)
#               $times          experimental times (vector)
#               $data_meas      experimental data (matrix)
#               $deg            degradation constants (vector)
#               $lag            lag time is estimated (boolean)
#       $opt                    optimization results returned by NLSIC (list)
#       $sens                   mean, median, standard deviation and 95% CI on the free parameters,
#                               estimated using a Monte-Carlo analysis (matrix)
#
#
####################
### Output files ###
####################
#
# Results are stored in the following files (where xxx is the name of the input file) created in the
# 'xxx_res' folder of the working directory:
#
#       xxx.pdf         plot of simulated and experimental data
#       xxx_res.txt     optimal solution and goodness of fit (khi2 test)
#                       results of the sensitivity analysis (mean, median, SD, 95% CI) if mc=TRUE
#       xxx_log.txt     detailed information on the system (list of metabolites, initial parameters,
#                       constraints, etc) and nlsic results (parameters, jacobian, residuals, etc)
#
#
####################
### Changes log  ###
####################
#
#   2018-10-18 version 1.0.1
#       - check names of degradation constants before calculations
#
#   2018-04-30 version 1.0
#       - new features:
#          . non-enzymatic degradation of some compounds can be taken into account (assuming a 1st order
#            kinetics) using the 'deg'argument
#          . initial lag time before start of growth and substrate consumption can be considered and
#            estimated using the 'lag' argument
#       - progress of monte-carlo iterations are displayed as a progress bar
#
#   2017-02-28 version 0.9.1
#       - raise an error if the 'weight' matrix contains 0.
#
#   2016-09-10 version 0.9
#       - evaluate the goodness of fit (khi2 test)
#       - minor changes in the pdf report and result txt file
#
#   2016-08-18 version 0.8
#       - the argument 'mode' has been removed, the two types of plot are generated in the pdf report
#       - additional plot of the differences between simulated and experimental data
#
#   2016-04-25 version 0.7
#       - weights can also be provided as a named vector
#       - handling missing values (given as NAs in the input file) which are not used during fitting
#
#   2016-02-23 version 0.6
#       - code optimized (simulations are performed using analytical functions instead of ODEs, no
#         need to indicate if metabolites are products or substrates)
#
#   2016-02-23 version 0.5
#       - experimental SDs can be provided to physio_fit() using the 'weight' argument
#       - results are saved in a folder *_res of the working directory
#       - metabolite concentrations can now be plotted on several y axes (one per metabolite)
#         to facilitate visualization when ranges of concentrations between metabolites widely
#         differ (e.g. µM vs mM)
#       - when Monte-Carlo analysis is performed, simulations using the mean of estimated parameters
#         are also plotted (dotted line)
#
#   2016-02-18 version 0.4
#       - biomass concentration is plotted on a second y axis (log scale)
#       - legend panel replaced by colored axes titles
#       - bug fix: in some situations, parameters returned by monte carlo analysis were mixed up
#
#   2016-02-17 version 0.3
#       - code optimized
#       - now returns an error if a required package or the input file is missing
#
#   2016-02-16 version 0.2
#       - monte-carlo sensitivity analysis is now included in the physio_fit() function
#       - keys of the list returned by physio_fit() have been changed
#       - changes in output files content (more details are included in the log file)
#       - changes in the plot of exp vs sim data
#       - bug fix: noisy datasets for MC analyses where generated from the measured data,
#                  now they are generated from the dataset simulated from the best solution
#
#   2016-02-15 version 0.1
#       - initial release
#
####################
####################


################################################
### Load libraries                           ###
################################################

load_pack <- function(packages){
    # load an R package, raise an error if not installed
    for (i in packages){
        if (!suppressWarnings(require(i, character.only=TRUE, quietly=TRUE))){
            stop("The package '", i, "' is missing. To install this package, please run the following command:\n   install.packages('", i, "')\n", sep="")
        }
    }
    return(invisible(NULL))
}

# load R packages
load_pack(c("nnls", "numDeriv", "RColorBrewer"))

# load NLSIC
result = tryCatch({source("nlsic.R")},
                  warning = function(w) {stop("NLSIC cannot be loaded, please check that 'nlsic.R' is in the PhysioFit directory", sep="")},
                  error = function(e) {stop("NLSIC cannot be loaded, please check that 'nlsic.R' is in the PhysioFit directory", sep="")},
                  finally = {cat('', sep="")})

# create the default color palette function
fun_col <- colorRampPalette(rev(brewer.pal(9, "Set1")))

# 'not in' function
"%ni%" <- Negate("%in%")


################################################
### Load functions                           ###
################################################


createSys <- function(data, vini=0.1, sd_X=0.002, sd_M=0.5, upcf=50, locf=-50, upcc=50, locc=1e-6, weight=NULL, sd_file=NULL, deg=NULL, lag=FALSE) {
    # Check input data and initialize parameters, constraints, and some variables
    if (!("time" %in% colnames(data))){
        stop("column 'time' not found, please check the input data file\n", sep="")
    }
    # get variables names
    nconc <- colnames(data)[colnames(data) != "time"]
    # extract metabolites (substrates and products) and fluxes
    metab <- nconc[nconc != "X"]
    nflux <- c("mu", paste("q", metab, sep=""))
    # check names of 'deg'
    if (!is.null(deg)){
        ndeg <- names(deg)
        # if no names are provided for degradation constants
        if (is.null(ndeg)){
            stop("degradation constants should be provided as a *named* vector\n", sep="")
        }
        # if some degradation constants are not named
        unnamed <- ndeg[ndeg == ""]
        if (length(unnamed)){
            stop("names should be provided for *all* degradation constants\n", sep="")
        }
        # if some degradation constants are not in metab
        not_in <- setdiff(names(deg), metab)
        if (length(not_in)){
            warning("a degradation constant is provided for metabolite(s) [", paste(not_in, collapse=","), "], but no measurements found in the datafile\n", sep="")
        }
    }
    # parameters
    to_est <- c(nconc, nflux)
    params <- rep(vini, length(to_est))
    # add lag
    if (lag){
        to_est <- c(to_est, "tlag")
        params <- c(params, max(data[,"time"])*0.2)
    }
    names(params) <- to_est
    # construct weight matrix
    if (is.null(weight)){
        # try to load weights (if not provided) from file 'datfile_sd.txt'
        if (file.exists(sd_file)){
            weight <- as.matrix(read.table(sd_file, sep="\t", header=TRUE))
            if (ncol(weight) != length(nconc)){
                stop("matrix 'weight' should have ", length(nconc), " columns, check the file '", sd_file, "'\n", sep="")
            }
            if (nrow(weight) != nrow(data)){
                stop("matrix 'weight' should have ", nrow(data), " rows, check the file '", sd_file, "'\n", sep="")
            }
        }else{
        # set weights to their default value
            weight <- matrix(1, nrow=nrow(data), ncol=length(nconc), dimnames=list(row=NULL, col=nconc))
            weight[,"X"] <- sd_X
            weight[,metab] <- sd_M
        }
        # or create the weight matrix from the vector provided
    }else if (is.vector(weight)){
        if (!all(nconc %in% names(weight))){
            stop("weight of all the metabolite(s) quantified must be provided in the 'weight' vector\n", sep="")
        }
        weight <- matrix(rep(weight[nconc], nrow(data)), nrow=nrow(data), ncol=length(nconc), byrow=TRUE, dimnames=list(row=NULL, col=nconc))
    }
    # raise an error if 'weight' contains 0
    if (any(weight <= 0)){
        stop("all weight values must be strictly positive\n", sep="")
    }
    # upper and lower constraints
    te_upc <- c(rep(upcc, length(nconc)), rep(upcf, length(nflux)))
    if (lag){
        te_upc <- c(te_upc, max(data[,"time"])-1e-6)
    }
    names(te_upc) <- to_est
    te_loc <- c(rep(locc, length(nconc)), 0, rep(locf, length(nflux)-1))
    if (lag){
        te_loc <- c(te_loc, 1e-6)
    }
    names(te_loc) <- to_est
    # construct matrix u and vector co such that u%*%param-co>=0 (luco$u, luco$co)
    luco <- uplo2uco(params, te_upc, te_loc)
    # gather items
    res <- list("params"    = params,
                "to_est"    = to_est,
                "nconc"     = nconc,
                "nflux"     = nflux,
                "metab"     = metab,
                "weight"    = weight[,nconc],
                "te_upc"    = te_upc,
                "te_loc"    = te_loc,
                "u"         = luco$u,
                "co"        = luco$co,
                "times"     = data[,"time"],
                "data_meas" = data[,nconc],
                "nb_par"    = length(to_est),
                "nb_conc"   = length(nconc),
                "deg"       = deg,
                "lag"       = lag)
    return(res)
}

sim_deg <- function(sys, times, params){
    # Calculate M(t) (where M is the concentration of the metabolite M that is degraded assuming 
    # a first order kinetics) using the following analytical functions:
    #     M(t) = M_0 * exp(-k*t)
    sim <- matrix(NA, nrow=length(times), ncol=sys$nb_conc, dimnames=list(row=NULL, col=sys$nconc))
    for (s in sys$metab){
        sim[,s] <- params[s] * exp(-params[paste("k", s, sep="")] * times)
    }
    return(sim)
}

sim_physio <- function(sys, times, params){
    # Calculate X(t) and M(t) (where X is the biomass concentration and M is the concentration of
    # the metabolite M) using the following analytical functions:
    #     X(t) = X_0 * exp(mu * t)
    #     M(t) = (qM * X_0 / mu) * (exp(mu * t) - 1) + M_0
    # note: qM is positive (negative) when M is produced (consummed)
    sim <- matrix(NA, nrow=length(times), ncol=sys$nb_conc, dimnames=list(row=NULL, col=sys$nconc))
    exp_mu_t <- exp(params["mu"] * times)
    sim[,"X"] <- params["X"] * exp_mu_t
    for (s in sys$metab){
        sim[,s] <- (params[paste("q", s, sep="")] * params["X"] / params["mu"]) * (exp_mu_t - 1) + params[s]
    }
    return(sim)
}

sim_physio_lag <- function(sys, times, params){
    # Calculate X(t) and M(t) (where X is the biomass concentration and M is the concentration of
    # the metabolite M) using the following analytical functions:
    #     if t<lag:  X(t) = X_0
    #                M(t) = M_0 * exp(-k * t)
    #     if t>=lag: X(t) = X_0 * exp(mu * (t-lag))
    #                M(t) = (qM * X_0 / mu) * (exp(mu * (t-lag)) - 1) + M_0
    # note: qM is positive (negative) when M is produced (consummed)
    sim <- matrix(NA, nrow=length(times), ncol=sys$nb_conc, dimnames=list(row=NULL, col=sys$nconc))
    exp_mu_t <- exp(params["mu"] * (times-params["tlag"]))
    tiu <- (times>=params["tlag"])
    til <- !tiu
    sim[til,"X"] <- params["X"]
    sim[tiu,"X"] <- params["X"] * exp_mu_t[tiu]
    for (s in sys$metab){
        sim[tiu,s] <- (params[paste("q", s, sep="")] * params["X"] / params["mu"]) * (exp_mu_t[tiu] - 1) + params[s]
        sim[til,s] <- params[s]
    }
    return(sim)
}

sim_physio_deg <- function(sys, times, params){
    # Calculate X(t) and M(t) (where X is the biomass concentration and M is the concentration of
    # the metabolite M) using the following analytical functions:
    #     X(t) = X_0 * exp(mu * t)
    #     for M which are degraded:     M(t) = M_0 * exp(-k * t) + (qM * X_0 * exp(k * t - mu * t) / (mu + k))
    #     for M which are not degraded: M(t) = (qM * X_0 / mu) * (exp(mu * (t-lag)) - 1) + M_0
    # note: qM is positive (negative) when M is produced (consummed)
    sim <- matrix(NA, nrow=length(times), ncol=sys$nb_conc, dimnames=list(row=NULL, col=sys$nconc))
    exp_mu_t <- exp(params["mu"] * times)
    sim[,"X"] <- params["X"] * exp_mu_t
    for (s in sys$metab){
        if (s %in% names(sys$deg)){
            exp_k_t <- exp(-sys$deg[s] * times)
            sim[,s] <- params[s]*exp_k_t + (params[paste("q", s, sep="")] * params["X"] * (exp_mu_t-exp_k_t) / (params["mu"] + sys$deg[s]))
        }else{
            sim[,s] <- (params[paste("q", s, sep="")] * params["X"] / params["mu"]) * (exp_mu_t - 1) + params[s]
        }
    }
    return(sim)
}

sim_physio_deg_lag <- function(sys, times, params){
    # Calculate X(t) and M(t) (where X is the biomass concentration and M is the concentration of
    # the metabolite M) using the following analytical functions:
    #     if t<lag:  X(t) = X_0
    #                for M which are degraded:     M(t) = M_0 * exp(-k * t)
    #                for M which are not degraded: M(t) = M_0
    #     if t>=lag: X(t) = X_0 * exp(mu * (t-lag))
    #                for M which are degraded:     M(t) = M_0 * exp(-k * lag) + (qM * X_0 * exp(k * (t-lag) - mu * (t-lag)) / (mu + k))
    #                for M which are not degraded: M(t) = (qM * X_0 / mu) * (exp(mu * (t-lag)) - 1) + M_0
    # note: qM is positive (negative) when M is produced (consummed)
    sim <- matrix(NA, nrow=length(times), ncol=sys$nb_conc, dimnames=list(row=NULL, col=sys$nconc))
    exp_mu_t <- exp(params["mu"] * (times-params["tlag"]))
    tiu <- (times>=params["tlag"])
    til <- !tiu
    sim[til,"X"] <- params["X"]
    sim[tiu,"X"] <- params["X"] * exp_mu_t[tiu]
    for (s in sys$metab){
        if (s %in% names(sys$deg)){
            sim[til,s] <- params[s]*exp(-sys$deg[s] * times[til])
            s0 <- params[s]*exp(-sys$deg[s] * params["tlag"])
            exp_k_t <- exp(-sys$deg[s] * (times[tiu]-params["tlag"]))
            sim[tiu,s] <- s0*exp_k_t + (params[paste("q", s, sep="")] * params["X"] * (exp_mu_t[tiu]-exp_k_t) / (params["mu"] + sys$deg[s]))
        }else{
            sim[tiu,s] <- (params[paste("q", s, sep="")] * params["X"] / params["mu"]) * (exp_mu_t[tiu] - 1) + params[s]
            sim[til,s] <- params[s]
        }
    }
    return(sim)
}

Cost <- function(fp, cjac=F, data, sys, sim_fun) {
    # Cost function for optimization
    sim <- sim_fun(sys, sys$times, fp)
    res <- as.numeric((sim - data)/sys$weight)
    if (any(is.na(res))){
        res <- res[!is.na(res)]
    }
    return(list(res=res, sim=sim))
}

test_khi2 <- function(nb_points, k_val, nb_par){
    # Perform Khi2 statistical test
    df <- nb_points - nb_par
    p_val <- pchisq(k_val, df=df)
    khi2test <- list("khi2 value"                  = k_val,
                     "data points"                 = nb_points,
                     "fitted parameters"           = nb_par,
                     "degrees of freedom"          = df,
                     "khi2 reduced value"          = k_val/df,
                     "p-value, i.e. P(X^2<=value)" = p_val)
    if (p_val > 0.95){
        khi2test$conclusion <- "At level of 95% confidence, the model does not fit the data good enough with respect to the provided measurement SD."
    }else{
        khi2test$conclusion <- "At level of 95% confidence, the model fits the data good enough with respect to the provided measurement SD."
    }
    return(khi2test)
}

physio_fit <- function(datfile, mc=TRUE, it=50, save=TRUE, vini=1, sd_X=0.002, sd_M=0.5, upcf=50, locf=-50, upcc=50, locc=1e-6, weight=NULL, pos=TRUE, alpha="33", col=NULL, deg=NULL, lag=FALSE, summary=FALSE){
    
    # Load experimental data from file 'datfile.txt'
    datfile <- tools::file_path_sans_ext(datfile)
    fin <- paste(datfile, ".txt", sep="")
    if (file.exists(fin)){
        data_meas <- as.matrix(read.table(fin, sep="\t", header=TRUE))
    }else{
        stop("file '", fin, "' does not exist.\n", sep="")
    }
    
    # Initialize parameters, constraints, etc
    cat("Initialize the system...\n")
    sys <- createSys(data_meas, vini=vini, sd_X=sd_X, sd_M=sd_M, upcf=upcf, locf=locf, upcc=upcc, locc=locc, weight=weight, sd_file=paste(datfile, "_sd.txt", sep=""), deg=deg, lag=lag)
    
    # Analytical function used for simulations
    if (lag == FALSE){
        if (is.null(deg)){
            sim_fun = sim_physio
        }else{
            sim_fun = sim_physio_deg
        }
    }else{
        if (is.null(deg)){
            sim_fun = sim_physio_lag
        }else{
            sim_fun = sim_physio_deg_lag
        }
    }
    
    # Fit experimental data
    cat("Estimate fluxes...\n")
    result <- nlsic(sys$params[sys$to_est], # parameters to estimate
                    Cost, # cost function
                    u=sys$u, co=sys$co, # constraints
                    control=list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
                    e=NULL, eco=NULL, flsi=lsi_ln,
                    sys$data_meas,
                    sys,
                    sim_fun)
    
    # Check if nlsic returned an error
    if (result$error!=0){
        stop("nlsic returned the following error:\n", result$mes, "\n", sep="")
    }

    # goodness of fit (khi2 test)
    khi2test <- test_khi2(nb_points=sum(!is.na(data_meas[,-1])), k_val=sum(result$retres$res**2), nb_par=sys$nb_par)
    
    # Calculate difference between measured and simulated data
    result$resi <- sys$data_meas - result$retres$sim
    
    # Run Monte-Carlo analysis
    if (mc){
        cat("Run sensitivity analysis...\n")
        # create matrices to store the results, vectors to generate noisy datasets, etc
        mc_res <- matrix(NA, nrow=it, ncol=sys$nb_par+1, dimnames=list(row=NULL, col=c(sys$to_est,"res")))
        sim_times <- sort(unique(c(sys$times, seq(from=min(sys$times), to=max(sys$times), length.out=100))))
        mc_sim <- array(NA, dim=c(it, sys$nb_conc, length(sim_times)), dimnames=list("it"= seq(it), "var"=sys$nconc, "t"=sim_times))
        data_ori <- as.numeric(result$retres$sim[,sys$nconc])
        sdd <- as.numeric(sys$weight)
        lori <- length(data_ori)
        # create the progress bar
        pb <- txtProgressBar(min=0, max=it, style=3)
        # run iterations
        for (i in seq(it)){
            # update the progress bar
            setTxtProgressBar(pb, i)
            # create noisy datasets
            noise_dat <- rnorm(lori, mean=data_ori, sd=sdd)
            m_noise <- matrix(noise_dat, ncol=sys$nb_conc, dimnames=list(row=NULL, col=sys$nconc))
            # set negative concentrations to 0
            if (pos){
                m_noise[m_noise < 0] <- 0
            }
            # run optimization on noisy data
            res <- nlsic(result$par,
                         Cost,
                         u=sys$u, co=sys$co,
                         control=list(errx=1.e-6, maxit=50, trace=0, ci=list(report=T), history=F, adaptbt=T),
                         e=NULL, eco=NULL, flsi=lsi_ln,
                         m_noise,
                         sys,
                         sim_fun)
            # save results
            if (res$error==0){
                # parameters
                mc_res[i,] <- c(res$par, sum(res$res**2))
                # simulations
                mc_sim[i,,] <- sim_fun(sys, sim_times, res$par)
            }else{
                cat("... an unexpected error has happened during optimization\n")
                cat("    ***nlsic results***\n")
                print(res)
            }
        }
        # close progress bar
        close(pb)
        # calculate mean, median, sd and 95% CI of estimated parameters
        out_mc <- matrix(c(apply(mc_res, 2, mean),
                           apply(mc_res, 2, median),
                           apply(mc_res, 2, sd),
                           t(apply(mc_res, 2, FUN=function(x) quantile(x, c(0.025, 0.975))))),
                         ncol=5,
                         dimnames=list(row=c(sys$to_est,"res"), col=c("mean","median","sd","CI_2.5","CI_97.5")))
        # calculate upper and lower 95% CI of the fitted curves
        out_mc_sim_min <- matrix(apply(mc_sim, c(2,3), FUN=function(x) quantile(x, probs=0.025)),
                                 ncol=sys$nb_conc,
                                 dimnames=list(row=sim_times, col=sys$nconc))
        out_mc_sim_max <- matrix(apply(mc_sim, c(2,3), FUN=function(x) quantile(x, probs=0.975)),
                                 ncol=sys$nb_conc,
                                 dimnames=list(row=sim_times, col=sys$nconc))
    }else{
        out_mc <- NULL
        out_mc_sim_min <- NULL
        out_mc_sim_max <- NULL
    }
    
    # Gather results
    res_opt <- list(sys=sys, opt=result, khi2test=khi2test, sens=out_mc, sens_sim=list(min=out_mc_sim_min, max=out_mc_sim_max))
    
    # Save results and create plot of exp. vs sim.
    if (save){
        
        cat("Save results and plot the data...\n")
        
        # create results directory
        mainDir <- getwd()
        subDir  <- paste(datfile, "_res", sep="")
        if (!file.exists(subDir)){
            dir.create(file.path(mainDir, subDir))
        }
        setwd(file.path(mainDir, subDir))
        
        # save results
        fout <- paste(datfile, "_res.txt", sep="")
        list2file(list(sys=sys, result=result), file=paste(datfile, "_log.txt", sep=""))
        
        if (mc){
            write(paste("Results (", it, " mc iterations)\n", sep=""), file=fout)
            out_rep <- cbind(matrix(c(result$par, sum(result$res**2)), ncol=1, dimnames=list(row=c(sys$to_est, "res"), col="optimal")), out_mc)
            suppressWarnings(write.table(out_rep, fout, sep="\t", append=TRUE, quote=FALSE))
        }else{
            write("Best solution\n", file=fout)
            suppressWarnings(write.table(result$par, file=fout, append=TRUE, quote=FALSE, col.names=FALSE, sep="\t"))
        }
        write("\n\nGoodness of fit (khi2 test)\n", fout, append=TRUE)
        suppressWarnings(write.table(as.matrix(khi2test), fout, sep="\t", append=TRUE, quote=FALSE, col.names=FALSE))
        
        # create pdf
        plot_pdf(res_opt, datfile=datfile, alpha=alpha, cust_col=col, sim_fun=sim_fun)
        
        # go back to the original working directory
        setwd(mainDir)
    }
    
    # Print khi2 test results
    if (summary){
        cat("khi2 test results:\n")
        cat("   khi2 value: ", khi2test$`khi2 value`, "\n", sep="")
        cat("   data points: ", khi2test$`data points`, "\n", sep="")
        cat("   fitted parameters: ", khi2test$`fitted parameters`, "\n", sep="")
        cat("   degrees of freedom: ", khi2test$`degrees of freedom`, "\n", sep="")
        cat("   khi2 reduced value: ", khi2test$`khi2 reduced value`, "\n", sep="")
        cat("   p-value: ", khi2test$`p-value, i.e. P(X^2<=value)`, "\n", sep="")
        cat("   conclusion: ", khi2test$`conclusion`, "\n", sep="")
    }
    cat("Done.\n")
    
    return(res_opt)
}

plot_pdf <- function(res_opt, datfile="fit", alpha="33", cust_col=NULL, sim_fun){
    
    # number of metabolites
    ls <- res_opt$sys$nb_conc
    
    # times for simulated curves (with at least 100 points for a nice plot)
    tsim <- sort(unique(c(res_opt$sys$times, seq(from=min(res_opt$sys$times), to=max(res_opt$sys$times), length.out=100))))
    
    # run simulations using parameters from best solution
    sim_opt <- sim_fun(res_opt$sys, tsim, res_opt$opt$par)
    
    # run simulations using parameters from "mean" solution
    if (!is.null(res_opt$sens)){
        sim_mc <- sim_fun(res_opt$sys, tsim, res_opt$sens[res_opt$sys$to_est, "mean"])
    }

    # custom color palette
    if (is.null(cust_col)){
        cust_col <- fun_col(ls)
    }
    
    # x axis limit
    xlim <- c(0, tsim[length(tsim)]*1.1)
    
    # create pdf
    pdf(file=paste(datfile, ".pdf", sep=""), width=6.5+ls*0.75, height=6)
    
    ### PAGE 0
    par(mfrow=c(2,3))
    for (k in colnames(res_opt$sys$data_meas)){
        ylim <- c(min(c(res_opt$sys$data_meas[,k], sim_opt[,k]), na.rm=TRUE)*0.8, max(c(res_opt$sys$data_meas[,k], sim_opt[,k]), na.rm=TRUE)*1.1)
        # plot measurements & simulations, with biomass on the right y axis (log scale)
        cci <- cust_col[which(colnames(res_opt$sys$data_meas) == k)]
        plot(x = res_opt$sys$times,
             y = res_opt$sys$data_meas[,k],
             type = "p",
             xlim = xlim,
             ylim = ylim,
             pch = 21,
             cex = 1.4,
             bg = cci,
             lwd = 2,
             col = "black",
             xlab = "time",
             ylab = k,
             las=1)
        lines(x = tsim,
              y = sim_opt[,k],
              col = cci,
              lwd = 1.4)
        if (!is.null(res_opt$sens)){
            polygon(x=c(tsim, rev(tsim)), y=c(res_opt$sens_sim$max[,k], rev(res_opt$sens_sim$min[,k])), col=paste(cci, alpha, sep=""), border=NA)
            lines(x = tsim,
                  y = sim_mc[,k],
                  col = cci,
                  lty = 2,
                  lwd = 1.4)
        }
    }
    
    ### PAGE 1
        par(mar=c(5.5, ls*3-2, 2, 4) + 0.1, mfrow=c(1,1))
        ylim <- c(min(c(res_opt$sys$data_meas[,"X"], sim_opt[,"X"]), na.rm=TRUE)*0.8, max(c(res_opt$sys$data_meas[,"X"], sim_opt[,"X"]), na.rm=TRUE)*1.1)
        # plot measurements & simulations, with biomass on the right y axis (log scale)
        plot(x = res_opt$sys$times,
             y = res_opt$sys$data_meas[,"X"],
             type = "p",
             xlim = xlim,
             ylim = ylim,
             log = "y",
             pch = 21,
             cex = 1.4,
             lwd = 2,
             bg = cust_col[1],
             col = "black", axes = FALSE, bty = "n", xlab = "", ylab = "")
        lines(x = tsim,
              y = sim_opt[,"X"],
              col = cust_col[1],
              lwd = 1.4)
        axis(4, las=1, cex=0.9, line=-2)
        mtext(side=4, text="biomass", line=1.2, cex=1, col=cust_col[1])
        if (!is.null(res_opt$sens)){
            polygon(x=c(tsim, rev(tsim)), y=c(res_opt$sens_sim$max[,"X"], rev(res_opt$sens_sim$min[,"X"])), col=paste(cust_col[1], alpha, sep=""), border=NA)
            lines(x = tsim,
                  y = sim_mc[,"X"],
                  lwd = 1.4,
                  lty = 2,
                  col = cust_col[1])
        }
        # plot metabolite concentrations
        for (ic in seq(ls-1)){
            val <- res_opt$sys$metab[ic]
            ylim <- c(0, max(c(res_opt$sys$data_meas[,val], sim_opt[,val]), na.rm=TRUE)*(1+(ic-1)*0.15))
            par(new=TRUE)
            plot(x = tsim,
                 y = sim_opt[,val],
                 xlim = xlim,
                 ylim = ylim,
                 lty = 1,
                 las = 1,
                 type = "l",
                 lwd = 1.4,
                 axes = FALSE,
                 col = cust_col[ic+1],
                 ylab = "",
                 xlab = "",
                 main = "")
            if (!is.null(res_opt$sens)){
                polygon(x=c(tsim, rev(tsim)), y=c(res_opt$sens_sim$max[,val], rev(res_opt$sens_sim$min[,val])), col=paste(cust_col[ic+1], alpha, sep=""), border=NA)
                lines(x = tsim,
                      y = sim_mc[,val],
                      lty = 2,
                      lwd = 1.4,
                      col = cust_col[ic+1])
            }
            # plot measurements
            points(x = res_opt$sys$times,
                   y = res_opt$sys$data_meas[,val],
                   pch = 21,
                   cex = 1.4,
                   bg = cust_col[ic+1],
                   col = "black",
                   lwd = 2)
            axis(2, ylim=ylim, line=(ic-1)*3, las=1, cex=0.9, col=cust_col[ic+1])
            mtext(2, text=val, line=ic*3-2.2, adj=-0.2, cex=0.97, col=cust_col[ic+1])
        }
        axis(1, xlim=xlim, las=1)
        mtext(1, text="time", line=2.5, adj=0.5, cex=1.2)

    ### PAGE 2
        par(mar=c(5.5, 4, 2, 5.5) + 0.1)
        ylim <- c(min(c(res_opt$sys$data_meas[,"X"], sim_opt[,"X"]), na.rm=TRUE)*0.8, max(c(res_opt$sys$data_meas[,"X"], sim_opt[,"X"]), na.rm=TRUE)*1.1)
        # plot measurements & simulations, with biomass on the right y axis (log scale)
        plot(x = res_opt$sys$times,
             y = res_opt$sys$data_meas[,"X"],
             type = "p",
             xlim = xlim,
             ylim = ylim,
             log = "y",
             pch = 21,
             cex = 1.4,
             bg = cust_col[1], lwd = 2,
             col = "black", axes = FALSE, bty = "n", xlab = "", ylab = "")
        lines(x = tsim,
              y = sim_opt[,"X"],
              col = cust_col[1], lwd = 1.4)
        axis(4, las=1, cex = 1.4)
        mtext(side=4, text="biomass", line=3.2, cex=1.2, col=cust_col[1])
        if (!is.null(res_opt$sens)){
            polygon(x=c(tsim, rev(tsim)), y=c(res_opt$sens_sim$max[,"X"], rev(res_opt$sens_sim$min[,"X"])), col=paste(cust_col[1], alpha, sep=""), border=NA)
            lines(x = tsim,
                  y = sim_mc[,"X"],
                  col = cust_col[1],
                  lty = 2,
                  lwd = 1.4)
        }
        # plot metabolite concentrations
        par(new=TRUE)
        ylim <- c(0, max(c(as.numeric(res_opt$sys$data_meas[,res_opt$sys$metab]), as.numeric(sim_opt[,res_opt$sys$metab])), na.rm=TRUE)*1.1)
        matplot(x = tsim,
                y = sim_opt[,res_opt$sys$metab],
                xlim = xlim,
                ylim = ylim,
                lty = 1,
                las = 1,
                type = "l",
                lwd = 1.4,
                cex.lab = 1.2,
                cex.main = 1.2,
                col = cust_col[-1],
                ylab = "",
                xlab = "time",
                main = "")
        if (!is.null(res_opt$sens)){
            for (ic in seq(ls-1)){
                val <- res_opt$sys$metab[ic]
                polygon(x=c(tsim, rev(tsim)), y=c(res_opt$sens_sim$max[,val], rev(res_opt$sens_sim$min[,val])), col=paste(cust_col[ic+1], alpha, sep=""), border=NA)
            }
            matlines(x = tsim,
                     y = sim_mc[,res_opt$sys$metab],
                     xlim = xlim,
                     ylim = ylim,
                     lty = 2,
                     lwd = 1.4,
                     col = cust_col[-1])
        }
        # plot measurements
        matpoints(x = res_opt$sys$times,
                  y = res_opt$sys$data_meas[,res_opt$sys$metab],
                  pch = 21,
                  cex = 1.4, lwd = 2,
                  bg = cust_col[-1],
                  col = "black")
        # add legend as colored axis label
        for (j in seq(ls-1)){
            if (j>1){
                ph1 = res_opt$sys$metab[seq(j-1)]
                ins1 = ", "
            }else{
                ph1 = c()
                ins1 = ""
            }
            if (j<(ls-1)){
                ph2 = res_opt$sys$metab[seq(from=j+1, to=ls-1)]
                ins2 = ", "
            }else{
                ph2 = c()
                ins2 = ""
            }
            te <- paste("expression(phantom('", paste(ph1, collapse=", "), ins1, "')*'", res_opt$sys$metab[j],"'*phantom('",ins2, paste(ph2, collapse=", "), " concentration'))", sep="")
            mtext(side=2, text=eval(parse(text=te)), line=2.5, cex=1.1, col=cust_col[j+1])
        }
        mtext(side = 2,
              text = eval(parse(text=paste("expression(", paste(paste(paste("phantom('", res_opt$sys$metab, "')", sep=""), collapse="*', '*"), "*' concentration'", sep=""), ")", sep=""))),
              line = 2.5, cex = 1.1, col = "black")

    ### PAGE 3
        par(mar=c(5.5, 4, 2, 5.5) + 0.1)
        ymin <- min(res_opt$opt$resi[,"X"], na.rm=TRUE)
        ymax <- max(res_opt$opt$resi[,"X"], na.rm=TRUE)
        ylim <- c(ymin*(1-sign(ymin)*0.2), ymax*(1+sign(ymax)*0.2))
        # plot measurements & simulations, with biomass on the right y axis (log scale)
        plot(x = res_opt$sys$times,
             y = res_opt$opt$resi[,"X"],
             type = "l",
             xlim = xlim,
             ylim = ylim,
             pch = 21,
             cex = 1.4,
             lwd = 1.4,
             col = cust_col[1], axes = FALSE, bty = "n", xlab = "", ylab = "", main="measured - simulated")
        abline(h=0)
        points(x = res_opt$sys$times,
               y = res_opt$opt$resi[,"X"],
               col = "black",
               bg = cust_col[1], lwd = 2,
               pch = 21,
               cex = 1.4)
        axis(4, las=1, cex = 1.4)
        mtext(side=4, text="biomass", line=3.2, cex=1.2, col=cust_col[1])
        # plot metabolite concentrations
        par(new=TRUE)
        ymin <- min(res_opt$opt$resi[,res_opt$sys$metab], na.rm=TRUE)
        ymax <- max(res_opt$opt$resi[,res_opt$sys$metab], na.rm=TRUE)
        ylim <- c(ymin*(1-sign(ymin)*0.2), ymax*(1+sign(ymax)*0.2))
        matplot(x = res_opt$sys$times,
                y = res_opt$opt$resi[,res_opt$sys$metab],
                xlim = xlim,
                ylim = ylim,
                lty = 1,
                las = 1,
                type = "l",
                lwd = 1.4,
                cex.lab = 1.2,
                cex.main = 1.2,
                col = cust_col[-1],
                ylab = "",
                xlab = "time",
                main = "")
        matpoints(x = res_opt$sys$times,
                  y = res_opt$opt$resi[,res_opt$sys$metab],
                  pch = 21,
                  cex = 1.4, lwd = 2,
                  bg = cust_col[-1],
                  col = "black")
        # add legend as colored axis label
        for (j in seq(ls-1)){
            if (j>1){
                ph1 = res_opt$sys$metab[seq(j-1)]
                ins1 = ", "
            }else{
                ph1 = c()
                ins1 = ""
            }
            if (j<(ls-1)){
                ph2 = res_opt$sys$metab[seq(from=j+1, to=ls-1)]
                ins2 = ", "
            }else{
                ph2 = c()
                ins2 = ""
            }
            te <- paste("expression(phantom('", paste(ph1, collapse=", "), ins1, "')*'", res_opt$sys$metab[j],"'*phantom('",ins2, paste(ph2, collapse=", "), "'))", sep="")
            mtext(side=2, text=eval(parse(text=te)), line=2.5, cex=1.1, col=cust_col[j+1])
        }
        mtext(side = 2,
              text = eval(parse(text=paste("expression(", paste(paste(paste("phantom('", res_opt$sys$metab, "')", sep=""), collapse="*', '*"), "*''", sep=""), ")", sep=""))),
              line = 2.5, cex = 1.1, col = "black")
        
    dev.off()
        
    return(invisible(NULL))
}

uplo2uco <- function(param, upper=NULL, lower=NULL) {
    # Return a list with a matrix u and a vector co such that u%*%param-co>=0
    # translates the inequalities param <= upper and param >= lower
    # Names in upper and lower must be present in names of param
    u=matrix(0., nrow=length(upper)+length(lower), ncol=length(param))
    co=numeric(nrow(u))
    colnames(u)=names(param)
    rownames(u)=c(paste(names(upper), " <= ", upper, sep=""),
                  paste(names(lower), " >= ", lower, sep=""))
    names(co)=rownames(u)
    # fill u and co
    u[iseq(length(upper)), names(upper)]=diag(-1, length(upper))
    co[iseq(length(upper))]=-upper
    u[length(upper)+iseq(length(lower)), names(lower)]=diag(1, length(lower))
    co[length(upper)+iseq(length(lower))]=lower
    return(list(u=u, co=co))
}

iseq <- function(n) {
    # positive sequence of integer numbers 1:n
    # if n=0 then integer(0) is returned
    seq.int(from=1, to=n, length=n)
}

list2file <- function(l, file = paste(deparse(substitute(l)), ".txt", sep = "")) {
    # save the list 'l' in a file
    tmp <- getOption("width")
    options(width=10000)
    sink(file)
    print(l)
    sink()
    options(width=tmp)
    return(invisible(NULL))
}

estimate_k <- function(datfile, save=TRUE, vini=0.01, sd_M=0.5, upcf=50, locf=0, upcc=50, locc=1e-6, weight=NULL, pos=TRUE, alpha="33", col=NULL) {
    # Estimate first-order degradation constant for a list of compounds
    
    # Load experimental data from file 'datfile.txt'
    datfile <- tools::file_path_sans_ext(datfile)
    fin <- paste(datfile, ".txt", sep="")
    if (file.exists(fin)){
        data_meas <- as.matrix(read.table(fin, sep="\t", header=TRUE))
    }else{
        stop("file '", fin, "' does not exist.\n", sep="")
    }
    
    # Initialize parameters, constraints, etc
    cat("Initialize the system...\n")
    sys <- createDeg(data_meas, vini=vini, sd_M=sd_M, upcf=upcf, locf=locf, upcc=upcc, locc=locc, weight=weight, sd_file=paste(datfile, "_sd.txt", sep=""))

    # Fit experimental data
    cat("Estimate fluxes...\n")
    result <- nlsic(sys$params[sys$to_est], # parameters to estimate
                    Cost, # cost function
                    u=sys$u, co=sys$co, # constraints
                    control=list(errx=1.e-6, maxit=50, trace=1, ci=list(report=T), history=F, adaptbt=T),
                    e=NULL, eco=NULL, flsi=lsi_ln,
                    sys$data_meas,
                    sys,
                    sim_deg)
    
    # Check if nlsic returned an error
    if (result$error!=0){
        stop("nlsic returned the following error:\n", result$mes, "\n", sep="")
    }
    
    # Goodness of fit (khi2 test)
    khi2test <- test_khi2(nb_points=sum(!is.na(data_meas[,-1])), k_val=sum(result$retres$res**2), nb_par=sys$nb_par)
    
    # Calculate difference between measured and simulated data
    result$resi <- sys$data_meas - result$retres$sim
    
    # Gather results
    deg_cst <- result$par[sys$nflux]
    tmp <- regmatches(names(deg_cst), regexec('(^.)(.*)', names(deg_cst)))
    names(deg_cst) <- sapply(tmp, "[", 3)
    res_opt <- list(sys=sys, opt=result, khi2test=khi2test, deg_cst=deg_cst)
    
    # Save results and create plot of exp. vs sim.
    if (save){
        
        cat("Save results...\n")
        
        # create results directory
        mainDir <- getwd()
        subDir <- paste(datfile, "_res", sep="")
        if (!file.exists(subDir)){
            dir.create(file.path(mainDir, subDir))
        }
        setwd(file.path(mainDir, subDir))
        
        # save results
        fout <- paste(datfile, "_res.txt", sep="")
        list2file(list(sys=sys, result=result), file=paste(datfile, "_log.txt", sep=""))
        
        write("Best solution\n", file=fout)
        suppressWarnings(write.table(result$par, file=fout, append=TRUE, quote=FALSE, col.names=FALSE, sep="\t"))
        
        write("\n\nGoodness of fit (khi2 test)\n", fout, append=TRUE)
        suppressWarnings(write.table(as.matrix(khi2test), fout, sep="\t", append=TRUE, quote=FALSE, col.names=FALSE))
        
        # create pdf
        #plot_pdf(res_opt, datfile=datfile, alpha=alpha, cust_col=col, sim_fun=sim_deg)
        
        # go back to the original working directory
        setwd(mainDir)
    }
    
    cat("Done.\n")
    
    return(res_opt)
}


createDeg <- function(data, vini=0.1, sd_M=0.5, upcf=50, locf=0, upcc=50, locc=1e-6, weight=NULL, sd_file=NULL) {
    # Initialize parameters, constraints, and some variables
    # check the 'time' column is present
    if (!("time" %in% colnames(data))){
        stop("column 'time' not found, please check the input data file\n", sep="")
    }
    # get variables names
    nconc <- colnames(data)[colnames(data) != "time"]
    # extract metabolites (substrates and products) and degradation fluxes
    metab <- nconc
    nflux <- paste("k", metab, sep="")
    # parameters
    to_est <- c(nconc, nflux)
    params <- rep(vini, length(to_est))
    names(params) <- to_est
    # construct weight matrix
    if (is.null(weight)){
        # try to load weights (if not provided) from file 'datfile_sd.txt'
        if (file.exists(sd_file)){
            weight <- as.matrix(read.table(sd_file, sep="\t", header=TRUE))
            if (ncol(weight) != length(nconc)){
                stop("matrix 'weight' should have ", length(nconc), " columns, check the file '", sd_file, "'\n", sep="")
            }
            if (nrow(weight) != nrow(data)){
                stop("matrix 'weight' should have ", nrow(data), " rows, check the file '", sd_file, "'\n", sep="")
            }
        }else{
            # set weights to their default value
            weight <- matrix(1, nrow=nrow(data), ncol=length(nconc), dimnames=list(row=NULL, col=nconc))
            weight[,metab] <- sd_M
        }
        # or create the weight matrix from the vector provided
    }else if (is.vector(weight)){
        if (!all(nconc %in% names(weight))){
            stop("weight of all the metabolite(s) quantified must be provided in the 'weight' vector\n", sep="")
        }
        weight <- matrix(rep(weight[nconc], nrow(data)), nrow=nrow(data), ncol=length(nconc), byrow=TRUE, dimnames=list(row=NULL, col=nconc))
    }
    # raise an error if 'weight' contains 0
    if (any(weight <= 0)){
        stop("all weight values must be strictly positive\n", sep="")
    }
    # constraints
    te_upc <- c(rep(upcc, length(nconc)), rep(upcf, length(nflux)))
    names(te_upc) <- to_est
    te_loc <- c(rep(locc, length(nconc)), 0, rep(locf, length(nflux)-1))
    names(te_loc) <- to_est
    # construct matrix u and vector co such that u%*%param-co>=0 (luco$u, luco$co)
    luco <- uplo2uco(params, te_upc, te_loc)
    # gather items
    res <- list("params"    = params,
                "to_est"    = to_est,
                "nconc"     = nconc,
                "nflux"     = nflux,
                "metab"     = metab,
                "weight"    = weight[,nconc],
                "te_upc"    = te_upc,
                "te_loc"    = te_loc,
                "u"         = luco$u,
                "co"        = luco$co,
                "times"     = data[,"time"],
                "data_meas" = data[,nconc],
                "nb_par"    = length(to_est),
                "nb_conc"   = length(nconc))
    return(res)
}


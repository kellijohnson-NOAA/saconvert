#' Convert ICES input files to Stock Synthesis (SS) input files
#'
#' A series of files are used as input to ICES models.
#' These files are standard input files for many types of models.
#' \code{ICES2SS} offers a standardized way to convert these ICES files
#' to files that can be used in a Stock Synthesis (SS) model.
#'
#' @param user.wd A file path to the directory containing the ICES files.
#' @param user.od A file path to a directory where the resulting files will
#' be saved.
#' @param ices.id A character value to set the ices.id because, sometimes,
#' there is no stock id at the beginning of the file names.
#' @param slx A numerical value specifying the age-based selectivity type that
#' will be used in Stock Synthesis.
#' @param tvslx A logical value specifying whether or not to implement
#' time-varying selectivity.
#' @param ages Two values, a min and max age, used for reporting fishing
#' mortality estimates for each year. These values will be placed in the
#' Stock Synthesis starter file for F std report ages.
#' @param forN An integer value specifying the number of forecast years for
#' the projections.
#' @param Q_sd A logical value indicating if an added standard deviation
#' should be estimated.
#' The parameters will be fixed at zero and will not be estimated if
#' \code{FALSE}. An additional standard deviation for each survey included
#' in the model will be estimated if \code{TRUE}, which is the default.
#'
#' @export
#' @author Kelli Faye Johnson
#'
ICES2SS <- function(user.wd, user.od, ices.id = "",
  slx = 26, tvslx = TRUE, ages = c(1, 2),
  forN = 1, Q_sd = TRUE,
  div2mt = 1) {

  #### todo items
  # todo: options for how to implement selectivity in the forecasts
  # todo: utilize discard data

  if (!lastcharacter(user.wd, .Platform$file.sep)) {
    user.wd <- file.path(user.wd, .Platform$file.sep)
  }
  if (!lastcharacter(user.od, .Platform$file.sep)) {
    user.od <- file.path(user.od, .Platform$file.sep)
  }
  dir.create(user.od, showWarnings = FALSE, recursive = TRUE)

  # Read in ICES data
  # (01) Catch in numbers
  # (02) Catch mean weight
  # (03) Discard mean weight
  # (04) Spawning mean weight
  # (05) Survey
  # (06) Natural mortality
  # (07) Proportion mature
  # (08) Maturity ogive
  # (09) Landed fraction
  # (10) Landed mean weight
  # (11) pf
  cn <- read.ices(paste(user.wd,ices.id,"cn.dat",sep=""))
  attr(cn, "time") <- 0.5
  cw <- read.ices(paste(user.wd,ices.id,"cw.dat",sep=""))
  dw <- read.ices(paste(user.wd,ices.id,"dw.dat",sep=""))
  sw <- read.ices(paste(user.wd,ices.id,"sw.dat",sep=""))
  surveys <- read.ices(paste(user.wd,ices.id,"survey.dat",sep=""))
  nm <- read.ices(paste(user.wd,ices.id,"nm.dat",sep=""))
  pm <- read.ices(paste(user.wd,ices.id,"pm.dat",sep=""))
  mo <- read.ices(paste(user.wd,ices.id,"mo.dat",sep=""))
  lf <- read.ices(paste(user.wd,ices.id,"lf.dat",sep=""))
  lw <- read.ices(paste(user.wd,ices.id,"lw.dat",sep=""))
  #propf <- read.ices(paste(user.wd,ices.id,"pf.dat",sep=""))
  t.spawn <- pm[1,1] #assuming time/age invariant spawning time

  # Catch of a single fleet with catches each modelled year
  # setting Freport as (catch.nages):(catch.nages-1); unweighted F
  nfleets <- 1
  catch.yy <- range(as.numeric(rownames(cn)))
  catch.yrs <- catch.yy[1]:catch.yy[2]
  catch.nyrs <- length(catch.yrs)
  catch.nages <- dim(cn)[2]
  catch.ages <- min(as.numeric(colnames(cn))):max(as.numeric(colnames(cn)))
  tot.catch <- apply(cn*cw,1,sum)

  waa.array <- array(c(
    cw[1:catch.nyrs,],  dw[1:catch.nyrs,], sw[1:catch.nyrs,]),
    dim = c(catch.nyrs, catch.nages, 3))
  # ASAP uses a weight-at-age pointer vector to assign three weight-at-age
  # matrices to the following matrices needed for the model
  # catch, discards, SSB, Jan-1 B, and indices

  n.surveys <- length(surveys)
  recr.CV <- rep(0.5, catch.nyrs)
  catch.CV <- rep(0.1, catch.nyrs)
  disc.CV <- rep(0, catch.nyrs)
  Neff.catch <- rep(100, catch.nyrs)
  Neff.disc <- rep(0, catch.nyrs)

  fleet.names <- paste0("fleet", seq(nfleets))
  survey.names <- gsub("\\s", "", names(surveys))

  #### Starter
  starter <- r4ss::SS_readstarter(verbose = FALSE,
    file = dir(utils::tail(dir(system.file("extdata", package = "r4ss"),
    pattern = "simple", full.names = TRUE), 1),
    pattern = "starter", full.names = TRUE))
  starter$sourcefile <- paste0(user.od, "starter.ss")
  starter$datfile <- "data.ss"
  starter$ctlfile <- "control.ss"
  starter$parmtrace <- 1
  starter$F_age_range <- ages
  starter$F_report_basis <- 0
  r4ss::SS_writestarter(starter, dir = user.od,
    overwrite = TRUE, warn = FALSE, verbose = FALSE)

  #### Forecast
  forecast <- r4ss::SS_readforecast(verbose = FALSE,
    file = dir(utils::tail(dir(system.file("extdata", package = "r4ss"),
    pattern = "simple", full.names = TRUE), 1),
    pattern = "forecast", full.names = TRUE))
  forecast$sourcefile <- paste0(user.od, "forecast.ss")
  forecast$benchmarks <- 1
  forecast$MSY <- 2
  # forecast$SPRtarget
  # forecast$Btarget
  forecast$Bmark_years <- rep(c(-999, 0), 5)
  forecast$Bmark_relF_Basis <- 2
  forecast$Forecast <- 4
  forecast$Nforecastyrs <- forN
  # forecast$F_scalar
  forecast$Fcast_years <- rep(c(-999, 0), 3)
  forecast$Fcast_selex <- 0
  forecast$ControlRuleMethod <- 1
  # forecast$BforconstantF
  # forecast$BfornoF
  forecast$Flimitfraction <- 1
  forecast$FirstYear_for_caps_and_allocations <- max(catch.yrs) + 1
  r4ss::SS_writeforecast(forecast, dir = user.od,
    overwrite = TRUE, verbose = FALSE)

  #### Data
  data <- r4ss::SS_readdat(verbose = FALSE,
    file = dir(utils::tail(dir(system.file("extdata", package = "r4ss"),
    pattern = "simple", full.names = TRUE), 1),
    pattern = "data", full.names = TRUE))
  ctl <- r4ss::SS_readctl(verbose = FALSE,
    file = dir(utils::tail(dir(system.file("extdata", package = "r4ss"),
      pattern = "simple", full.names = TRUE), 1),
      pattern = "control", full.names = TRUE),
    use_datlist = TRUE, datlist = data)
  data$sourcefile <- paste0(user.od, "data.ss")
  data$styr <- catch.yy[1]
  data$endyr <- catch.yy[2]
  # data$nseas
  # data$months_per_seas
  # data$Nsubseasons
  # data$spawn_month
  data$Nages <- catch.nages
  # data$Nareas
  data$fleetnames <- c(fleet.names, survey.names)
  data$Nfleets <- 1 + length(survey.names)
  data$fleetinfo <- data.frame(
    "type" = c(
      rep(1, length(fleet.names)),
      rep(3, length(survey.names))),
    "surveytiming" = c(
      rep(-1, length(fleet.names)), # throughout year
      rep(1, length(survey.names))),# read in month
    "area" = 1,
    "units" = c(
      rep(1, length(fleet.names)),   # mt
      rep(2, length(survey.names))), # thousands of fish
    "need_catch_mult" = 0,
    "fleetname" = data$fleetnames
  )
  data$surveytiming <- c(
    rep(-1, length(fleet.names)),
    rep(1, length(survey.names)))
  # Catch in mt
  # "equilibrium" catch == first year of catches with a high se
  data$units_of_catch <- data$fleetinfo[, "units"]
  data$catch <- data.frame(
    "year" = c(-999, as.numeric(names(tot.catch))),
    "seas" = 1,
    "fleet" = 1,
    "catch" = c(tot.catch[1], tot.catch) / div2mt,
    "catch_se" = c(0.5, rep(0.01, length(tot.catch))))
  data$CPUEinfo <- data.frame(
    "Fleet" = 1:data$Nfleets,
    "Units" = 1,
    "Errtype" = 0,
    "SD_Report" = 0)
  row.names(data$CPUEinfo) <- data$fleetnames
  get_survey <- function(x, fleet, se = 0.2) {
    index <- apply(x, 1, sum, na.rm = TRUE) / div2mt
    timing <- mean(attributes(x)[["time"]]) * 12
    result <- data.frame(
      "year" = as.numeric(row.names(x)),
      "seas" = timing,
      "index" = fleet,
      "obs" = index,
      "se_log" = se)
    result <- result[result[, "obs"] != 0, ]
    return(result)
  }
  get_age <- function(x, fleet, part = 0, maxage) {
    timing <- mean(attributes(x)[["time"]]) * 12
    xx <- matrix(0, nrow = nrow(x), ncol = maxage)
    xx[, 1:ncol(x)] <- x
    result <- data.frame(
      "Yr" = as.numeric(row.names(x)),
      "seas" = timing,
      "FltSvy" = fleet,
      "Gender" = 0,
      "Part" = part,
      "Ageerr" = 1,
      "Lbin_lo" = -1,
      "Lbin_hi" = -1,
      "Nsamp" = 50,
      xx)
    result[is.na(result)] <- 0
    result <- result[
      !(apply(result[, -c(1:which(colnames(result) == "Nsamp"))], 1, sum)==0), ]
    return(result)
  }
  # todo: determine SE of CPUE data
  data$CPUE <- do.call("rbind", lapply(seq_along(surveys),
    function(xx) get_survey(surveys[[xx]], xx + 1)))
  # todo: get discard data
  # todo: get mean-body-weight data
  # todo: deal with size data
  # data$lbin_method
  # data$binwidth
  # data$minimum_size
  # data$maximum_size
  # data$use_lencomp
  data$len_info <- data.frame(
    "mintailcomp" = 0,
    "addtocomp" = 1e-07,
    "combine_M_F" =  0,
    "CompressBins" = 0,
    "CompError" = rep(0, length(data$fleetnames)),
    "ParmSelect" = 0,
    "minsamplesize" = 0.001)
  row.names(data$len_info) <- data$fleetnames
  # data$N_lbins
  # data$lbin_vector
  data$lencomp <- NULL
  data$N_agebins <- dim(cn)[2]
  data$agebin_vector <- as.numeric(colnames(cn))
  data$N_ageerror_definitions <- 1
  data$ageerror <- matrix(c(
    seq(0.5, data$Nages + 0.5, by = 1),
    rep(0, data$Nages + 1)), nrow = 2, byrow = TRUE)
  data$age_info <- data$len_info
  data$age_info$CompError <- 1
  data$age_info$ParmSelect <- 1:(length(survey.names) + 1)
  data$Lbin_method <- 3
  # todo: deal with surveys that don't sample age 0
  data$agecomp <- rbind(
    do.call("rbind", lapply(seq_along(surveys),
    function(xx) get_age(surveys[[xx]], xx + 1, maxage = max(catch.ages)))),
    get_age(cn, fleet = 1, maxage = max(catch.ages)))
  data$use_MeanSize_at_Age_obs <- 0
  data$MeanSize_at_Age_obs <- NULL
  # data$N_environ_variables
  # data$N_sizefreq_methods
  # data$do_tags
  # data$morphcomp_data
  # data$use_selectivity_priors
  # data$eof
  # data$spawn_seas
  data$Nfleet <- 1
  data$Nsurveys <- data$Nfleets - 1
  data$N_areas <- 1
  data$Ngenders <- -1
  # data$N_cpue
  # data$fleetinfo1
  # data$fleetinfo2
  # data$N_meanbodywt
  # data$comp_tail_compression
  # data$add_to_comp
  # data$max_combined_lbin
  # data$N_lbinspop
  # data$lbin_vector_pop
  r4ss::SS_writedat(data, verbose = FALSE, outfile = paste0(user.od, "data.ss"),
    overwrite = TRUE)

  #### Control
  ctl$N_Block_Designs <- 0
  ctl$blocks_per_pattern <- NULL
  ctl$Block_Design <- NULL
  ctl$EmpiricalWAA <- 1
  ctl$Growth_Age_for_L1 <- 1
  ctl$Growth_Age_for_L2 <- max(catch.ages)
  ctl$maturity_option <- 5 # disable maturity
  ctl$MainRdevYrFirst <- min(data$agecomp[data$agecomp$FltSvy != 1, "Yr"])
  ctl$MainRdevYrLast <- catch.yy[2] - min(catch.ages)
  ctl$recdev_phase <- 3
  ctl$recr_dist_method <- 4
  ctl$recdev_early_start <- catch.yy[1] - catch.nages
  ctl$recdev_early_phase <- 4
  ctl$Fcast_recr_phase <- 0
  # todo: make sure other models are not using bias adjustment
  ctl$last_early_yr_nobias_adj <- catch.yy[1] - 1
  ctl$first_yr_fullbias_adj <- min(data$agecomp[data$agecomp$FltSvy != 1, "Yr"])
  ctl$last_yr_fullbias_adj <- catch.yy[2]
  ctl$first_recent_yr_nobias_adj <- catch.yy[2] + 1
  ctl$max_bias_adj <- -1
  ctl$MG_parms <- ctl$MG_parms[-grep("RecrDist|_Mal_", rownames(ctl$MG_parms)), ]
  ctl$MG_parms[grep("L_at|VonB|CV", rownames(ctl$MG_parms)), "PHASE"] <- -6
  # Fishing mortality
  ctl$F_ballpark <- 1.0
  ctl$F_ballpark_year <- catch.yy[1]
  ctl$F_Method <- 2
  ctl$maxF <- 3
  ctl$F_iter <- NULL
  # F setup of Baranov continuous F, which is recommended in high F cases:
  # Initial F, Phase, Number of additional lines
  # Equilibrium F for each fishery: 7 parameter line
  ctl$F_setup <- c(0.5, 2, 0)
  ctl$init_F <- data.frame("LO" = 0, "HI" = 4,
    "INIT" = 0.5, "PRIOR" = 0.5, "PR_SD" = 99, "PR_type" = 6,
    "PHASE" = 2)
  # Selectivity
  ctl$size_selex_types <- setNames(data.frame(matrix(0, ncol = 4,
    nrow = nfleets + n.surveys, byrow = TRUE)),
    colnames(ctl$size_selex_types))
  ctl$age_selex_types <- setNames(data.frame(matrix(
    rep(ctl$age_selex_types[1, ], nfleets + n.surveys),
    nrow = nfleets + n.surveys, byrow = TRUE)),
    colnames(ctl$age_selex_types))
  ctl$age_selex_types[, 2:3] <- 0
  ctl$age_selex_types[, 1] <- 17
  ctl$age_selex_types[, 4] <- catch.nages
  ctl$age_selex_parms <- data.frame(
    "LO" = c(-10002, -1, rep(-10, catch.nages - 1)),
    "HI" = c(1, rep(10, catch.nages)),
    "INIT" = c(-1000, 0, rep(0.1, catch.nages - 1)),
    "PRIOR" = 0, "SD" = 0, "PR_TYPE" = 0,
    "PHASE" = c(-5, -5, rep(5, catch.nages - 1)),
    matrix(0, ncol = 7, nrow = catch.nages + 1))
  ctl$age_selex_parms <- do.call("rbind",
    replicate(length(survey.names) + 1,
      ctl$age_selex_parms, simplify = FALSE))
  if (slx == 26) { # Exponential logistic
    ctl$age_selex_types[1, 1] <- 26
    ctl$age_selex_parms <- rbind(data.frame(
      "LO" = rep(0.001, 3),
      "HI" = c(1, 1, 0.5),
      "INIT" = c(0.1, 0.5, 0.01),
      "PRIOR" = 0, "SD" = 0, "PR_TYPE" = 0, "PHASE" = 5,
      matrix(0, ncol = 7, nrow = 3)),
    ctl$age_selex_parms[-(0:catch.nages + 1), ])
  }
  dmpars <- do.call("rbind", replicate(length(survey.names) + 1,
    data.frame("LO" = -5, "HI" = 20, "INIT" = 3,
    "PRIOR" = 0, "SD" = 1.813, "PR_TYPE" = 6, "PHASE" = 6,
    0, 0, 0, 0, 0, 0, 0),
    simplify = FALSE))
  colnames(dmpars) <- colnames(ctl$size_selex_parms)
  colnames(ctl$age_selex_parms) <- colnames(ctl$size_selex_parms)
  ctl$age_selex_parms <- rbind(ctl$age_selex_parms, dmpars)
  ctl$size_selex_parms <- NULL
  if (tvslx) {
    ctl$Use_2D_AR1_selectivity <- 1
    ctl$specs_2D_AR <- data.frame("fleet" = 1,
      "ymin" = catch.yy[1],
      "ymax" = catch.yy[2],
      "amin" = min(catch.ages),
      "amax" = max(catch.ages),
      "sigma_amax" = 1, "use_rho" = 1, "len1/age2" = 2,
      "devphase" = 5, "before_range" = 0, "after_range" = 0)
    ctl$pars_2D_AR <- data.frame(
      "LO" = c(0, -1, -1),
      "HI" = c(4, 1, 1),
      "INIT" = c(1, 0, 0),
      "PRIOR" = c(1, 0, 0),
      "SD" = 0.1, "PR_TYPE" = 6, "PHASE" = -5)
  }
  # todo: implement added SD for all surveys
  ctl$Q_options <- data.frame(
    "fleet" = 2:(length(survey.names) + 1),
    "link" = 1, "link_info" = 0, "extra_se" = 1,
    "biasadj" = 0, "float" = 1)

  ctl$Q_setup <- data.frame(
    "Den_dep" = 0, "env_var" = 0,
    "extra_se" = c(0, rep(1, length(survey.names))),
    "Q_type" = c(0, rep(2, length(survey.names))))

# todo: check Q phase
  ctl$Q_parms <- data.frame(
    "LO" = rep(c(-15, 0.001), times = length(survey.names)),
    "HI" = rep(c(15, 2.0), times = length(survey.names)),
    "INIT" = rep(c(0.1, ifelse(Q_sd, 0.05, 0.0)), times = length(survey.names)),
    "PRIOR" = rep(c(0.0, 0.08), times = length(survey.names)),
    "PRIOR_SD" = rep(c(-1.0, 0.1), times = length(survey.names)),
    "PRIOR_type" = 0,
    "PHASE" = rep(c(-1, ifelse(Q_sd, 4, -4)), times = length(survey.names)),
    "env_var&link" = 0,
    "dev_link" = 0,
    "dev_minyr" = 0,
    "dev_maxyr" = 0,
    "dev_PH" = 0,
    "Block" = 0,
    "Block_Fxn" = 0
  )

  # Remove lambdas b/c using DM parameters and extra SD
  ctl$N_lambdas <- 0
  ctl$lambdas <- NULL

  # todo: this assumes a time-invariant fixed natural mortality
  ctl$natM_type <- 3
  matage <- c(mean(nm[1:catch.nyrs,1]), apply(nm[1:catch.nyrs,], 2, mean))
  ctl$natM <- as.data.frame(matage)
  ctl$MG_parms <- ctl$MG_parms[-grep("NatM", rownames(ctl$MG_parms)), ]
  ctl$MG_parms$PHASE <- -2

  # Fix steepness at 1 and sigma_R at 0.8
  ctl$SR_parms[grep("LN", rownames(ctl$SR_parms)), "HI"] <- 20
  ctl$SR_parms[grep("sigma", rownames(ctl$SR_parms)), "INIT"] <- 0.8
  ctl$SR_parms[grep("steep", rownames(ctl$SR_parms)), "INIT"] <- 1
  ctl$SR_parms[grep("steep", rownames(ctl$SR_parms)), "PHASE"] <- -1

  ctl$more_stddev_reporting <- 0
  # Report selectivity, Age-based selectivty, Ending year selectivity,
  # Number of selectivity bins #todo think about what bins (ages 1-4?),
  # growth pattern, growth bins, Area for Natage, Natage Year, Natage bins
  # todo: sort out growth vs natage stdreporting
  ctl$stddev_reporting_specs <- c(1, 2, -1, catch.nages, 1, 1, 1, -1, 1)
  ctl$stddev_reporting_selex <- seq(1, catch.nages)
  ctl$stddev_reporting_growth <- seq(1, catch.nages)
  ctl$stddev_reporting_N_at_A <- -1

  r4ss::SS_writectl(ctl, outfile = paste0(user.od, "control.ss"),
    overwrite = TRUE, verbose = FALSE)

  # Change array of wtatage to data frame with fleet
  # -2 Fecundity -1 population 0 population 1+ fleets and surveys
  waa.pop <- do.call("rbind",
    replicate(n.surveys + 2, data.frame(
    "Yr" = catch.yrs,
    "Seas" = 1, "Sex" = 1, "Bio_Pattern" = 1, "BirthSeas" = 1,
    "Fleet" = -1,
    "0" = waa.array[, 1, 3], waa.array[, , 3]), simplify = FALSE))
  waa.pop[, "Fleet"] <- rep(c(-1, 0, seq(n.surveys)+nfleets), each = catch.nyrs)
  waa.fleets <- do.call("rbind",
    replicate(nfleets, data.frame(
    "Yr" = catch.yrs,
    "Seas" = 1, "Sex" = 1, "Bio_Pattern" = 1, "BirthSeas" = 1,
    "Fleet" = 1,
    "0" = waa.array[, 1, 1], waa.array[, , 1]), simplify = FALSE))
  waa.fleets[, "Fleet"] <- rep(seq(nfleets), each = catch.nyrs)
  usemo <- matrix(0, ncol = ncol(waa.array[, , 3]),
     nrow = nrow(waa.array[, , 3]))
  colnames(usemo) <- 1:NCOL(usemo)
  usemo[, which(colnames(usemo) %in% colnames(mo))] <- mo[
    rownames(mo) %in% catch.yrs, which(colnames(mo) %in% colnames(usemo))]
  fecmat <- waa.array[, , 3] * usemo
  waa.fec <- data.frame(
    "Yr" = catch.yrs,
    "Seas" = 1, "Sex" = 1, "Bio_Pattern" = 1, "BirthSeas" = 1,
    "Fleet" = -2,
    "0" = fecmat[, 1], fecmat)
  waa.new <- rbind(waa.fec, waa.pop, waa.fleets)
  waa.new[waa.new[, "Yr"] %in% catch.yy[2], "Yr"] <-
    waa.new[waa.new[, "Yr"] %in% catch.yy[2], "Yr"] * -1
  r4ss::SS_writewtatage(mylist = waa.new, dir = user.od,
    warn = FALSE, verbose = FALSE, overwrite = TRUE)

  # copy ss executable
  ignore <- file.copy(overwrite = TRUE, file.path(
    system.file("bin", "Windows64", package = "ss3sim"),
    "ss.exe"),
    paste0(user.od, "ss.exe"))
}

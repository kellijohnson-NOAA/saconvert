#' Convert ICES input files to Stock Synthesis (SS) input files
#' 
#' A series of files are used as input to ICES models.
#' These files are standard input files for many types of models.
#' This code acts as a standardized way to convert these ICES files
#' to files that can be used in a Stock Synthesis (SS) model.
#' 
#' @param user.wd A file path to the directory containing the ICES files.
#' @param user.od A file path to a directory where the resulting files will
#' be saved.
#' @param ices.id A character value to set the ices.id because, sometimes,
#' there is no stock id at the beginning of the file names
#' @param slx A numerical value specifying the age-based selectivity type that
#' will be used in Stock Synthesis.
#' @param tvslx A logical value specifying whether or not to implement
#' time-varying selectivity.
#' @param ages Two values, a min and max age used for reporting fishing
#' mortality estimates for each year. These values will be placed in the
#' Stock Synthesis starter file for F std report ages.
#' @param forN An integer value specifying the number of forecast years for
#' the projections.
#' 
#' @examples
#' ICES2SS(user.wd = "c:/stockAssessment/MGWG/state-space/NScod", 
#'   user.od = "c:/stockAssessment/MGWG/state-space/NScod/SS",
#'   slx = 17, tvslx = FALSE, ages = c(2,4))
#' output <- r4ss::SS_output(getwd(), covar = FALSE)
#' r4ss::SS_plots(output, uncertainty = FALSE)
#'
#' @export
#' @author Kelli Faye Johnson
#'
ICES2SS <- function(user.wd, user.od, ices.id = "", 
  slx = 26, tvslx = TRUE, ages = c(1, 2), 
  forN = 1) {
  
  #### todo items
  # todo: options for how to implement selectivity in the forecasts
  
  if (!lastcharacter(user.wd, .Platform$file.sep)) {
    user.wd <- file.path(user.wd, .Platform$file.sep)
  }
  if (!lastcharacter(user.od, .Platform$file.sep)) {
    user.od <- file.path(user.od, .Platform$file.sep)
  }
  dir.create(user.od, showWarnings = FALSE, recursive = TRUE)

  cn <- read.ices(paste(user.wd,ices.id,"cn.dat",sep=""))
  attr(cn, "time") <- 0.5
  cw <- read.ices(paste(user.wd,ices.id,"cw.dat",sep=""))
  dw <- read.ices(paste(user.wd,ices.id,"dw.dat",sep=""))
  #lf <- read.ices(paste(user.wd,ices.id,"lf.dat",sep=""))
  #lw <- read.ices(paste(user.wd,ices.id,"lw.dat",sep=""))
  mo <- read.ices(paste(user.wd,ices.id,"mo.dat",sep=""))
  nm <- read.ices(paste(user.wd,ices.id,"nm.dat",sep=""))
  #propf <- read.ices(paste(user.wd,ices.id,"pf.dat",sep=""))
  pm <- read.ices(paste(user.wd,ices.id,"pm.dat",sep=""))
  sw <- read.ices(paste(user.wd,ices.id,"sw.dat",sep=""))
  surveys <- read.ices(paste(user.wd,ices.id,"survey.dat",sep=""))
  t.spawn <- pm[1,1] #assuming time/age invariant spawning time
  
  catch.yy <- as.numeric(c(min(rownames(cn)), max(rownames(cn))))  # assuming there is only one CAA matrix
  nfleets <- 1                       # thus, also assuming nfleets=1
  catch.yrs <- seq(catch.yy[1], catch.yy[2])     #assuming catch defines the start/end year
  catch.nyrs <- length(catch.yrs)
  catch.nages <- dim(cn)[2]     #assuming catch matrix defines total number of modeled ages
  # setting Freport as (catch.nages):(catch.nages-1)  ; unweighted F
  catch.ages <- range(as.numeric(colnames(cn)))
  catch.ages <- seq(catch.ages[1],catch.ages[2],1)
  asap.ages <- 1:catch.nages
  # assuming 3 WAA matrices (catch, discard, and spawning weight); since assuming 1 fleet, cw should equal lw in ASAP)
  waa.array <- array(NA, dim=c(catch.nyrs, catch.nages, 3))
  waa.array[,,1] <- cw[1:catch.nyrs,]
  waa.array[,,2] <- dw[1:catch.nyrs,]
  waa.array[,,3] <- sw[1:catch.nyrs,]
  waa.pointer.vec <- c(1, 2, 1, 2, 3, 3) #assuming spawning weight-jan-1 biomass
  
  f.sel.blks <- rep(1, catch.nyrs) # assuming 1 selectivity block for all years (catch)
  f.sel.type <-  2 # assuming logistic (1=by age; 2=logistic; 3=double logistic)
  f.peak <- get.peak.age(cn)
  f.sel.mats.c1 <- c( seq(0.1,0.9, length.out=(catch.nages)),  round((f.peak)/2,2), 0.9,
                      round((f.peak)/4,2), 0.6,  round((catch.nages)/1.5,2), 1.1)
  f.sel.mats.c1[f.peak] <-1
  f.sel.mats.c2 <- c( rep(1, catch.nages), 2,3, rep(1, 4)) # phase for estimation
  f.sel.mats.c2[f.peak] <- -1
  f.sel.mats.c3 <- rep(0, (catch.nages+6)) # lambda for sel parameters
  f.sel.mats.c4 <- rep(1, (catch.nages+6)) # CV for sel parameters (irrelevant if lambda=0)
  f.sel.mats <- cbind(f.sel.mats.c1, f.sel.mats.c2, f.sel.mats.c3, f.sel.mats.c4)
  
  rel.mort.fleet <- rep(0,nfleets) # assuming release mortality at age (discard) is 0
  rel.prop <- matrix(0,nfleets*catch.nyrs, catch.nages)
  tot.catch <- apply(cn*cw,1,sum)
  
  
  n.surveys <- length(surveys)
  units.ind <- rep(2, n.surveys) # assuming unites=number (1=biomass; 2=number)
  time.ind <- get.survey.time(surveys)
  fish.ind <- rep(-1, n.surveys) #assuming none of the indices link to a fleet (i.e. all fishery-independent indices)
  index.sel.type <-  rep(2, n.surveys) #assuming logistic for simple setup
  ind.ages <- get.survey.ages(surveys)
  #ind.age1 <- sapply(ind.ages, min)
  ind.age1 <- rep(1, length(ind.ages))
  #ind.age2 <-  sapply(ind.ages, max)
  ind.age2 <- rep(length(catch.ages), length(ind.ages))
  ind.use <-  rep(1, n.surveys)
  i.peak <- get.peak.age(surveys)
  ind.sel.mats <- setup.surv.sel(surveys, i.peak, catch.ages, ind.ages)
  ind.cv = 0.2    # assume same CV for all years, all indices to setup ASAP indices matrix
  ind.neff = 50   # assume same Effective sample size for all years, all indices to setup ASAP indices matrix
  #ind.mat <- get.index.mat(x=surveys, a=ind.ages,  cv=0.2, neff=50)  #calculate total index and append CV and Neff columns
  #ind.mat = get.index.mat(x=surveys, cv = 0.2, neff = 50, first.year = catch.yy[1], nyears = catch.nyrs, catch.ages, survey.ages)  {
 
  recr.CV <- rep(0.5, catch.nyrs)
  catch.CV <- rep(0.1, catch.nyrs)
  disc.CV <- rep(0, catch.nyrs)
  Neff.catch <- rep(100, catch.nyrs)
  Neff.disc <- rep(0, catch.nyrs)
  Fmult.y1 <- 0.1
  #Get equilibrium numbers at age under F
  naa.y1 <- (nm[1,1]/(nm[1,1]+Fmult.y1))*cn[1,]/(1-exp(-nm[1,1]-Fmult.y1))
  if(naa.y1[1]==min(naa.y1) ) naa.y1[1] <-  10*mean(naa.y1)
  naa.y1[which(naa.y1==0)] <- mean(naa.y1)
  q.y1 <-  jitter(rep(0.05, n.surveys) , 30 )
  
  proj.yr=(catch.yy[2]+2) #dummy set up for 2 year projection
  proj.specs <- matrix(NA, nrow=2, ncol=5)
  proj.specs[,1] <- c((catch.yy[2]+1), (catch.yy[2]+2))
  proj.specs[,2] <- rep(-1, 2)
  proj.specs[,3]<- c(1,3)
  proj.specs[,4] <-  c(150,-99)
  proj.specs[,5] <-  rep(0,2)
  
  fleet.names <- "fleet1"
  survey.names <- gsub("\\s", "", names(surveys))
  fleet.dir <-  rep(1,nfleets)
  disc.flag =  F
  
  #### Starter
  starter <- r4ss::SS_readstarter(verbose = FALSE,
    file = dir(utils::tail(dir(system.file("extdata", package = "r4ss"), pattern = "simple", full.names = TRUE), 1), pattern = "starter", full.names = TRUE))
  starter$sourcefile <- paste0(user.od, "starter.ss")
  starter$datfile <- "data.ss"
  starter$ctlfile <- "control.ss"
  starter$F_age_range <- ages
  starter$F_report_basis <- 0
  r4ss::SS_writestarter(starter, dir = user.od, 
    overwrite = TRUE, warn = FALSE, verbose = FALSE)

  #### Forecast
  forecast <- r4ss::SS_readforecast(verbose = FALSE,
    file = dir(utils::tail(dir(system.file("extdata", package = "r4ss"), pattern = "simple", full.names = TRUE), 1), pattern = "forecast", full.names = TRUE))
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
    file = dir(utils::tail(dir(system.file("extdata", package = "r4ss"), pattern = "simple", full.names = TRUE), 1), pattern = "data", full.names = TRUE))
  data$sourcefile <- paste0(user.od, "data.ss")
  data$styr <- catch.yy[1]
  data$endyr <- catch.yy[2]
  # data$nseas
  # data$months_per_seas
  # data$Nsubseasons
  # data$spawn_month
  # data$Nsexes
  data$Nages <- catch.nages
  # data$Nareas
  data$Nfleets <- 1 + length(survey.names)
  data$fleetinfo[2:data$Nfleets, ] <- data$fleetinfo[2, ]
  data$fleetinfo[, "fleetname"] <- data$fleetnames <- c(fleet.names, survey.names)
  data$surveytiming <- c(-1, rep(1, data$Nfleets - 1))
  # todo: determine units of catch
  data$units_of_catch[1] <- data$fleetinfo[1, "units"] <- 2
  data$catch <- data.frame(
    "year" = as.numeric(row.names(cn)),
    "seas" = 1,
    "fleet" = 1,
    "catch" = apply(cn, 1, sum),
    "catch_se" = 0.01)
  data$CPUEinfo[2:data$Nfleets, ] <- data$CPUEinfo[2, ]
  row.names(data$CPUEinfo) <- data$fleetnames
  data$CPUEinfo[, "Fleet"] <- 1:data$Nfleets
  data$CPUE <- data.frame(
    "year" = as.numeric())
  get_survey <- function(x, fleet, se = 0.5) {
    index <- apply(x, 1, sum)
    timing <- mean(attributes(x)[["time"]]) * 12
    result <- data.frame(
      "year" = as.numeric(row.names(x)),
      "seas" = timing,
      "index" = fleet,
      "obs" = index,
      "se_log" = se)
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
      cbind(xx, xx))
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
  data$len_info[1:data$Nfleets, ] <- data$len_info[1, ]
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
  data$age_info[1:data$Nfleets, ] <- data$age_info[1, ]
  row.names(data$age_info) <- data$fleetnames
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
  data$Ngenders <- data$Nsexes
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
  ctl <- r4ss::SS_readctl(verbose = FALSE,
    file = dir(utils::tail(dir(system.file("extdata", package = "r4ss"), pattern = "simple", full.names = TRUE), 1), pattern = "control", full.names = TRUE),
    Ngenders = data$Nsexes, Nfleet = 1, Nsurveys = data$Nfleets - 1)
  ctl$EmpiricalWAA <- 1
  ctl$Growth_Age_for_L1 <- 1
  ctl$Growth_Age_for_L2 <- max(catch.nages)
  # ctl$maturity_option <- 5 # disable maturity
  ctl$MainRdevYrFirst <- catch.yrs[1]
  ctl$MainRdevYrLast <- utils::tail(catch.yrs, 1)
  ctl$recdev_phase <- 1
  ctl$recr_dist_method <- 4
  ctl$MG_parms <- ctl$MG_parms[-grep("RecrDist", rownames(ctl$MG_parms)), ]
  # todo: change early start year
  ctl$recdev_early_start <- catch.nages * -3
  ctl$recdev_early_phase <- 3
  ctl$Fcast_recr_phase <- 6
  ctl$last_early_yr_nobias_adj <- catch.yrs[1] - 1
  ctl$first_yr_fullbias_adj <- catch.yrs[1]
  ctl$last_yr_fullbias_adj <- utils::tail(catch.yrs, 1)
  ctl$first_recent_yr_nobias_adj <- utils::tail(catch.yrs, 1) + 1
  ctl$F_ballpark_year <- catch.yrs[1]
  # Selectivity
  ctl$size_selex_types[, 1:4] <- 0
  ctl$age_selex_types[, 2:4] <- 0
  ctl$age_selex_types[, 1] <- 17
  ctl$age_selex_types[, 4] <- catch.nages
  ctl$age_selex_parms <- data.frame(
    "LO" = c(-10002, -1, rep(-10, catch.nages - 1)),
    "HI" = c(1, rep(10, catch.nages)),
    "INIT" = c(-1000, 0, rep(0.01, catch.nages - 1)),
    "PRIOR" = 0, "SD" = 0, "PR_TYPE" = 0, 
    "PHASE" = c(-4, -4, rep(4, catch.nages - 1)), 
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
      "PRIOR" = 0, "SD" = 0, "PR_TYPE" = 0, "PHASE" = 4, 
      matrix(0, ncol = 7, nrow = 3)),
    ctl$age_selex_parms[-(0:catch.nages + 1), ])
  }
  dmpars <- do.call("rbind", replicate(length(survey.names) + 1, 
    data.frame("LO" = -7, "HI" = 7, "INIT" = 1,
    "PRIOR" = 0, "SD" = 0, "PR_TYPE" = 0, "PHASE" = 7,
    0, 0, 0, 0, 0, 0, 0), 
    simplify = FALSE))
  colnames(dmpars) <- colnames(ctl$size_selex_parms)
  colnames(ctl$age_selex_parms) <- colnames(ctl$size_selex_parms)
  ctl$age_selex_parms <- rbind(ctl$age_selex_parms, dmpars)
  ctl$size_selex_parms <- NULL
  if (tvslx) {
    ctl$Use_2D_AR1_selectivity <- 1
    ctl$specs_2D_AR <- data.frame("fleet" = 1, 
      "ymin" = min(catch.yrs), 
      "ymax" = max(catch.yrs),
      "amin" = min(catch.ages), 
      "amax" = max(catch.ages), 
      "sigma_amax" = 1, "use_rho" = 1, "len1/age2" = 2,
      "devphase" = 5, "before_range" = 0, "after_range" = 0)
    ctl$pars_2D_AR <- data.frame(
      "LO" = c(0, -1, -1),
      "HI" = c(4, 1, 1),
      "INIT" = c(1, 0, 0),
      "PRIOR" = c(1, 0, 0),
      "SD" = 0.1, "PR_TYPE" = 6, "PHASE" = -6)
  } 
  # todo: implement added SD for all surveys
  ctl$Q_options <- data.frame(
    "fleet" = 2:(length(survey.names) + 1),
    "link" = 1, "link_info" = 0, "extra_se" = 1,
    "biasadj" = 0, "float" = 1)
  
  ctl$Q_setup <- data.frame(
    "Den_dep" = 0, "env_var" = 0, 
    "extra_se" = c(0, rep(1, length(survey.names))),
    "extra_se" = c(0, rep(2, length(survey.names))),)
  # todo: this assumes a time-invariant fixed natural mortality
  ctl$natM_type <- 3
  matage <- c(mean(nm[1:catch.nyrs,1]), apply(nm[1:catch.nyrs,], 2, mean))
  ctl$natM <- as.data.frame(rbind(matage, matage))
  ctl$MG_parms <- ctl$MG_parms[-grep("NatM", rownames(ctl$MG_parms)), ]
  
  # Fix steepness at 1 and sigma_R at 0.5
  ctl$SRparm[grep("sigma", rownames(ctl$SRparm)), "INIT"] <- 0.5
  ctl$SRparm[grep("steep", rownames(ctl$SRparm)), "INIT"] <- 1
  ctl$SRparm[grep("steep", rownames(ctl$SRparm)), "PHASE"] <- -1
  
  ctl$stddev_reporting_selex[1] <- -1
  ctl$stddev_reporting_growth[1] <- -1
  ctl$stddev_reporting_N_at_A[1] <- -1
  
  r4ss::SS_writectl(ctl, outfile = paste0(user.od, "control.ss"),
    overwrite = TRUE, verbose = FALSE)
  # Change array of wtatage to data frame with fleet
  waa.new <- do.call("rbind", 
    replicate((length(survey.names) + 3), data.frame(
    "Yr" = catch.yrs,
    # todo: check season of weight at age
    "Seas" = 1,
    "Sex" = 1, 
    "Bio_Pattern" = 1,
    "BirthSeas" = 1, 
    "Fleet" = 1,
    "0" = waa.array[, 1, 1],
    waa.array[, , 1]), simplify = FALSE))
  waa.new$Fleet <- rep(-1:(length(survey.names) + 1),
    each = length(catch.yrs))
  # waa.new[NROW(waa.new), "Yr"] <- waa.new[NROW(waa.new), "Yr"] * -1
  fecmat <- waa.array[, , 3] * mo[rownames(mo) %in% catch.yrs, ]
  waa.fec <- data.frame(
    "Yr" = catch.yrs,
    "Seas" = 1, 
    "Sex" = 1, 
    "Bio_Pattern" = 1,
    "BirthSeas" = 1,
    "Fleet" = -2,
    "0" = fecmat[, 1],
    fecmat)
  waa.new <- rbind(waa.new, waa.fec)
  waa.mal <- waa.new
  waa.mal$Sex <- 2
  waa.new <- rbind(waa.new, waa.mal)
  waa.forecast <- waa.new[waa.new$Yr == max(catch.yrs), ]
  waa.forecast$Yr <- 1 + waa.forecast$Yr
  waa.new <- rbind(waa.new, waa.forecast)
  
  r4ss::SS_writewtatage(mylist = waa.new, dir = user.od, 
    warn = FALSE, verbose = FALSE, overwrite = TRUE)

}

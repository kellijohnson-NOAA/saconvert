#' Translating ICES file format to "vanilla" ASAP input file
#'
#' @param user.wd xxx
#' @param user.od xxx
#' @param model.id = model.id
#' model.id = ''
#' @param ices.id because sometimes there is no stock id at the beginning of the file names
#' 
ICES2ASAP <- function(user.wd,user.od,model.id,ices.id){ 
  cn <- read.ices(paste(user.wd,ices.id,"cn.dat",sep=""))
  print(cn)
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
  
  catch.yy <- as.numeric(c(min(rownames(cn)), max(rownames(cn)) ))  # assuming there is only one CAA matrix
  nfleets <- 1                       # thus, also assuming nfleets=1
  catch.yrs <- seq(catch.yy[1], catch.yy[2])     #assuming catch defines the start/end year
  catch.nyrs <- length(catch.yrs)
  catch.nages <- dim(cn) [2]     #assuming catch matrix defines total number of modeled ages
  # setting Freport as (catch.nages):(catch.nages-1)  ; unweighted F
  catch.ages <- range(as.numeric(colnames(cn)))
  catch.ages = seq(catch.ages[1],catch.ages[2],1)
  asap.ages = 1:catch.nages
  # assuming 3 WAA matrices (catch, discard, and spawning weight); since assuming 1 fleet, cw should equal lw in ASAP)
  waa.array <- array(NA, dim=c(catch.nyrs, catch.nages, 3))
  waa.array[,,1] <- cw[1:catch.nyrs,]
  waa.array[,,2] <- dw[1:catch.nyrs,]
  waa.array[,,3] <- sw[1:catch.nyrs,]
  waa.pointer.vec <- c(1, 2, 1, 2, 3, 3) #assuming spawning weight-jan-1 biomass
  
  f.sel.blks <- rep(1, catch.nyrs) # assuming 1 selectivity block for all years (catch)
  f.sel.type <-  2 # assuming logistic (1=by age; 2=logistic; 3=double logistic)
  f.peak <- get.peak.age(cn)
  print("f.peak")
  print(f.peak)
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
  print(surveys)
  i.peak <- get.peak.age(surveys)
  print("i.peak")
  print(i.peak)
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
  survey.names <- names(surveys)
  fleet.dir <-  rep(1,nfleets)
  disc.flag =  F
  
  # call the function to setup ASAP
  #phases indicated by p.(param.name) have been set at simple default
  #by default, steepness is fixed at 1 (estimates mean recruitment with deviations)
  #  to change this default, set "h.guess" to a value in [0.21, 0.99] and set p.h to positive integer
  setup.asap.w(wd=user.wd, od=user.od, model.id=model.id, nyears=catch.nyrs,
               first.year=catch.yy[1], asap.nages=catch.nages, nfleets=nfleets,
               nselblks=nfleets, n.ind.avail=n.surveys, M.mat=nm[1:catch.nyrs,],
               fec.opt=0, t.spawn=t.spawn, mat.mat=mo[1:catch.nyrs,],
               n.waa.mats=3, waa.array=waa.array, waa.pointer.vec=waa.pointer.vec,
               sel.blks=f.sel.blks, sel.types=f.sel.type, sel.mats=f.sel.mats,
               fleet.age1=asap.ages[1], fleet.age2=asap.ages[length(asap.ages)],
               F.report.ages=c((catch.nages-1),catch.nages), F.report.opt=1,
               like.const=0, rel.mort.fleet=rel.mort.fleet, caa.mats=cbind(cn, tot.catch), daa.mats=cbind(cn*0, 0*tot.catch),
               rel.prop=rel.prop, units.ind=units.ind, time.ind=time.ind,
               fish.ind=fish.ind, sel.ind=index.sel.type,
               ind.age1, ind.age2, ind.use, ind.sel.mats, ind.mat=surveys, 
               ind.cv=ind.cv, ind.neff=ind.neff,
               p.Fmult1=1, p.Fmult.dev=3, p.recr.dev=3, p.N1=2, p.q1=1, p.q.dev=-1, p.SR=1, p.h=-2,
               recr.CV=rep(0.5,catch.nyrs), lam.ind=rep(1,n.surveys),
               lam.c.wt=rep(1,nfleets), lam.disc=rep(0,nfleets), catch.CV=catch.CV, disc.CV=disc.CV,
               Neff.catch, Neff.disc, lam.Fmult.y1=rep(0, nfleets),
               CV.Fmult.y1=rep(1, nfleets), lam.Fmult.dev=rep(0,nfleets), CV.Fmult.dev=rep(1,nfleets),
               lam.N1.dev=0, CV.N1.dev=1, lam.recr.dev=1,
               lam.q.y1=rep(0, n.surveys), CV.q.y1=rep(1, n.surveys), lam.q.dev=rep(0, n.surveys),
               CV.q.dev=rep(1, n.surveys), lam.h=0, CV.h=1, lam.SSB0=0, CV.SSB0=1,
               naa.y1, Fmult.y1, q.y1, SSB0=1e7, h.guess=1.0, F.max=5, ignore.guess=0,
               do.proj=0, fleet.dir=fleet.dir, proj.yr=proj.yr, proj.specs=proj.specs,
               do.mcmc=0, mcmc.nyr.opt=0, mcmc.nboot=1000, mcmc.thin=200, mcmc.seed=5230547,
               recr.agepro=0, recr.start.yr=(catch.yy[2]-12), recr.end.yr=(catch.yy[2]-2),
               test.val=-23456, fleet.names, survey.names, disc.flag, catch.ages = catch.ages, survey.ages = ind.ages )
}

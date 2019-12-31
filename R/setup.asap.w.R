# Code to take ICES format and convert to ASAP
#   for ICES-WGMG projects
# assumptions for call to setup.asap begin ~line 742
# Liz Brooks
# Version 1.0
# also uses : SAM read.ices fn modified by Dan Hennen (starting line 422)
##


#rm(list=ls(all.names=F))
#graphics.off()


#==============================================================
## User specify below

#-------------------
#user.wd <- ""   #user: specify path to working directory where ICES files are
#user.od <- ""   #user: specify path to output directory
#model.id <- "CCGOMyt_"   # user: specify prefix found on ICES files (will create same name for ASAP case)
#-------------------
#user.wd <- "C:/liz/SAM/GBhaddock/"  # user: specify path to working directory where ICES files are
#user.od <- "C:/liz/SAM/GBhaddock/"  # user: specify path to output directory
#model.id <- "GBhaddock_"  # user: specify prefix found on ICES files (will create same name for ASAP case)
#-------------------
#user.wd <- "C:/liz/SAM/GBwinter/"  # user: specify path to working directory where ICES files are
#user.od <- "C:/liz/SAM/GBwinter/"  # user: specify path to output directory
#model.id <- "GBwinter_"  # user: specify prefix found on ICES files (will create same name for ASAP case)
#-------------------
#user.wd <- "C:/liz/SAM/Plaice/"  # user: specify path to working directory where ICES files are
#user.od <- "C:/liz/SAM/Plaice/"  # user: specify path to output directory
#model.id <- "Plaice_"  # user: specify prefix found on ICES files (will create same name for ASAP case)
#-------------------
#user.wd <- "C:/liz/SAM/NScod/"  # user: specify path to working directory where ICES files are
#user.od <- "C:/liz/SAM/NScod/"  # user: specify path to output directory
#model.id <- "ICEHerr_"  # user: specify prefix found on ICES files (will create same name for ASAP case)
   ## *** Notes: had to append "NScod_" to all ICES filenames
#-------------------
#user.wd <- "C:/liz/SAM/ICEherring/"  # user: specify path to working directory where ICES files are
#user.od <- "C:/liz/SAM/ICEherring/"  # user: specify path to output directory
#model.id <- "ICEherring_"  # user: specify prefix found on ICES files (will create same name for ASAP case)
 # *** Notes: only VPA files available now; need to convert to ICES format before running this
#-------------------

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Function to set-up asap3 "west coast style"
# Liz Brooks
# Version 1.0
# Created 30 September 2010
# Last Modified: 18 September 2013
#                16 November 2017 for ices-wgmg
#                21 November 2017: tested & works on CCGOMyt, GBhaddock, GBwinter, Plaice, NScod
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

#' @param wd         working directory path (where files are read from)
#' @param od         output directory path (where files are written)
#' @param model.id    model identifier
#' @param nyears      total number of years of data
#' @param first.year   first year of data
#' @param asap.nages       number of age classes (age 1 is first age class by default)
#' @param nfleets      number of fishing fleets
#' @param nselblks     total number of selectivity blocks (sum for all fleets)
#' @param n.ind.avail   number of available indices (whether or you "turn them on" to be used)
#' @param M.mat         matrix of natural mortality by age (col) and year (row)
#' @param fec.opt      0(use WAA*mat.age)  or 1 (use empirical fecundity at age values)
#' @param t.spawn      fraction of year elapsed prior to ssb calcs
#' @param mat.mat       maturity matrix by age (col) and year (row)
#' @param n.waa.mats xxx
#' @param waa.array xxx
#' @param waa.pointer.vec xxx
#' @param sel.blks    a vertical vector of nselblks*nyears
#' @param sel.types   vector of length nselblks (1=by age; 2= logistic; 3= double logistic)
#' @param sel.mats    nselblks X matrix(sel.specs, nrow= nages+6, ncol=4)
#' @param fleet.age1   starting age for selectivity by fleet
#' @param fleet.age2    ending age for selectivity by fleet
#' @param F.report.ages  vector of 2 ages for summarizing F trend
#' @param F.report.opt  option to report F as unweighted(1), Nweighted(2), Bweighted(3)
#' @param like.const    flag to use(1) or not(0) likelihood constants
#' @param rel.mort.fleet  flag for whether there is release mortality by fleet (nfleets entries)
#' @param caa.mats      nfleets X cbind(matrix(caa, nyears,nages), tot.cat.biomass)
#' @param daa.mats      nfleets X cbind(matrix(disc.aa, nyears, nages), tot.disc.biomass)
#' @param rel.prop      nfleets X matrix(release.prop.aa, nyears, nages)
#' @param units.ind     n.ind.avail vector for units (1=biomass, 2=number)
#' @param time.ind     n.ind.avail vector for month index sampled
#' @param fish.ind      link to fleet (-1 if no link, fleet.number otherwise)
#' @param sel.ind       functional form for indices (n.ind.avail)
#' @param ind.age1      first age each index selects (n.ind.avail)
#' @param ind.age2      last age each index selects (n.ind.avail)
#' @param ind.use       flag to use(1) or not(0) each index
#' @param ind.sel.mats  n.ind.avail X matrix(sel.specs, nrow= nages+6, ncol=4)
#' the 6 additional are: Units, month, sel.link.to.fleet, sel.start.age, sel.end.age, use.ind
#' @param ind.mat       n.ind.avail X matrix(index.stuff, nyears, ncol=nages+4)
#' ICES one-offs (calls function get.index.mat)
#' @param ind.cv        one-off for ICES (CV assumed for all indices, all years)
#' @param ind.neff      one-off for ICES (Effectice Number assumed for all indices, all years)
#' end ICES one-offs
#' @param p.Fmult1      phase for estimating F mult in 1st year
#' @param p.Fmult.dev   phase for estimating devs for Fmult
#' @param p.recr.dev    phase for estimating recruitment deviations
#' @param p.N1          phase for estimating N in 1st year
#' @param p.q1          phase for estimating q in 1st year
#' @param p.q.dev       phase for estimating q deviations
#' @param p.SR          phase for estimating SR relationship
#' @param p.h           phase for estimating steepness
#' @param recr.CV       vertical vector of CV on recruitment per year
#' @param lam.ind    lambda for each index
#' @param lam.c.wt      lambda for total catch in weight by fleet
#' @param lam.disc      lambda for total discards at age by fleet
#' @param catch.CV      matrix(CV.fleet, nyears, nfleets)
#' @param disc.CV       matrix(CV.fleet, nyears, nfleets)
#' @param Neff.catch    input effective sample size for CAA (matrix(Neff, nyears, nfleets)
#' @param Neff.disc     input effective sample size for disc.AA (matrix(Neff, nyears, nfleets)
#' @param lam.Fmult.y1  lambda for Fmult in first year by fleet (nfleets)
#' @param CV.Fmult.y1   CV for Fmult in first year by fleet (nfleets)
#' @param lam.Fmult.dev  lambda for Fmult devs by fleet (nfleets)
#' @param CV.Fmult.dev  CV for Fmult deviations by fleet (nfleets)
#' @param lam.N1.dev     lambda for N in 1st year devs
#' @param CV.N1.dev     CV for N in 1st year devs
#' @param lam.recr.dev  lambda for recruitment devs
#' @param lam.q.y1      lambda for q in 1st yr by index (n.ind.avail)
#' @param CV.q.y1       CV for q in 1st yr by index (n.ind.avail)
#' @param lam.q.dev      lambda for q devs (n.ind.avail)
#' @param CV.q.dev       CV for q devs (n.ind.avail)
#' @param lam.h          lambda for deviation from initial steepness
#' @param CV.h           CV for deviation from initial steepness
#' @param lam.SSB0       lambda for deviation from SSB0
#' @param CV.SSB0        CV for deviation from SSB0
#' @param naa.y1          vector(nages) of initial stock size
#' @param Fmult.y1       initial guess for Fmult in yr1 (nfleets)
#' @param q.y1           q in 1st year vector(n.ind.avail)
#' @param SSB0           initial unexploited stock size
#' @param h.guess        guess for initial steepness
#' @param F.max          upper bound on Fmult
#' @param ignore.guess   flag to ignore(1) or not(0) initial guesses
#' @param do.proj        flag to do(1) or not(0) projections
#' @param fleet.dir      rep(1,nfleets)
#' @param proj.yr        (nyears+2)
#' @param proj.specs     matrix(proj.dummy, nrow=2, ncol=5)
#' @param do.mcmc        0(no) or 1(yes)
#' @param mcmc.nyr.opt   0(use.NAA.last.yr), 1(use.NAA.T+1)
#' @param mcmc.nboot     number of mcmc iterations
#' @param mcmc.thin      thinning rate for mcmc
#' @param mcmc.seed      random number seed for mcmc routine
#' @param recr.agepro     0(use NAA), 1 (use S-R), 2(use geometric mean of previous years)
#' @param recr.start.yr  starting year for calculation of R
#' @param recr.end.yr    ending year for calculation of R
#' @param test.val       -23456
#' @param fleet.names xxx
#' @param survey.names xxx
#' @param disc.flag       T if discards present, F otherwise
#' @param catch.ages xxx
#' @param survey.ages xxx
setup.asap.w <-function(wd, od, model.id, nyears, first.year, asap.nages, nfleets,
        nselblks, n.ind.avail, M.mat, fec.opt, t.spawn, mat.mat, n.waa.mats, waa.array, waa.pointer.vec,
        sel.blks, sel.types, sel.mats, fleet.age1, fleet.age2, F.report.ages, F.report.opt,
        like.const, rel.mort.fleet, caa.mats, daa.mats, rel.prop, units.ind, time.ind,
        fish.ind, sel.ind, ind.age1, ind.age2, ind.use, ind.sel.mats, ind.mat, ind.cv, ind.neff,
        p.Fmult1, p.Fmult.dev, p.recr.dev, p.N1, p.q1, p.q.dev, p.SR, p.h, recr.CV, lam.ind,
        lam.c.wt, lam.disc, catch.CV, disc.CV, Neff.catch, Neff.disc, lam.Fmult.y1,
        CV.Fmult.y1, lam.Fmult.dev, CV.Fmult.dev, lam.N1.dev, CV.N1.dev, lam.recr.dev,
        lam.q.y1, CV.q.y1, lam.q.dev, CV.q.dev, lam.h, CV.h, lam.SSB0, CV.SSB0,
        naa.y1, Fmult.y1, q.y1, SSB0, h.guess, F.max, ignore.guess,
        do.proj, fleet.dir, proj.yr, proj.specs,
        do.mcmc, mcmc.nyr.opt, mcmc.nboot, mcmc.thin, mcmc.seed,
        recr.agepro, recr.start.yr, recr.end.yr, test.val,
        fleet.names, survey.names, disc.flag, catch.ages, survey.ages )    {


# c.waa         catch weight at age (col) and year (row)
# ssb.waa       ssb weight at age (col) and year (row)
# jan1.waa      jan-1 weight at age (col) and year (row)


#---------------------------------------------------------------------
####   SET-UP ASAP FILE

#_________________________________________________________________

out.file = paste(od,"ASAP_", model.id, ".dat", sep="")
write('# ASAP VERSION 3.0 setup by convert_ICES_asap.r', file=out.file, append=F)
write(paste('# MODEL ID ', model.id, sep=''),file=out.file,append=T)
write( '# Number of Years' , file=out.file,append=T)
write(nyears, file=out.file,append=T )
write('# First year', file=out.file,append=T)  #proportion F before spawning
write(first.year, file=out.file,append=T )  #proportion M before spawning
write('# Number of ages', file=out.file,append=T)  #single value for M
write(asap.nages, file=out.file,append=T )  #last year of selectivity
write('# Number of fleets', file=out.file,append=T)     #last year of maturity
write(nfleets, file=out.file,append=T )  #last year of catch WAA
write('# Number of selectivity blocks', file=out.file,append=T)    #last year of stock biomass
write(nselblks, file=out.file,append=T )  #number of F grid values
write('# Number of available indices', file=out.file,append=T)  #
write(n.ind.avail, file=out.file,append=T )  #specifies BH or Ricker
write( '# M matrix' , file=out.file,append=T)   #, ncolumns=(nyears))
write(t(M.mat), file=out.file,append=T, ncolumns=asap.nages)
write('# Fecundity option', file=out.file,append=T)  #specifies normal or lognormal error
write(fec.opt, file=out.file,append=T)  #
write('# Fraction of year elapsed before SSB calculation', file=out.file,append=T)  #
write(t.spawn , file=out.file,append=T)  #
write( '# MATURITY matrix' , file=out.file,append=T)   #, ncolumns=(nyears))
write(t(mat.mat), file=out.file,append=T, ncolumns=asap.nages)
write( '# Number of WAA matrices' , file=out.file,append=T)   #, ncolumns=(nyears))
write(n.waa.mats, file=out.file,append=T, ncolumns=asap.nages)
write( '# WAA matrix-1' , file=out.file,append=T)   #, ncolumns=(nyears))
write(t(waa.array[,,1]), file=out.file,append=T, ncolumns=asap.nages)
if (n.waa.mats>1)  {
    for (j in 2:n.waa.mats)  {
    write(paste('# WAA matrix-',j, sep=""), file=out.file,append=T, ncolumns=asap.nages)
    write(t(waa.array[,,j]), file=out.file,append=T, ncolumns=asap.nages)

    } # end loop over j (for WAA matrices)
}  # end if-test for n.waa.mat
#write('# test', file=out.file,append=T)
write( '# WEIGHT AT AGE POINTERS' , file=out.file,append=T)   #, ncolumns=(nyears))
write(waa.pointer.vec, file=out.file,append=T, ncolumns=1)
write( '# Selectivity blocks (blocks within years)' , file=out.file,append=T)   #, ncolumns=(nyears))
for(i in 1:nfleets) 
{
  write(paste0('# Fleet ', i, ' Selectivity Block Assignment') , file=out.file,append=T)   #, ncolumns=(nyears))
  write(sel.blks[(i-1)*nyears + 1:nyears], file=out.file,append=T, ncolumns=1)
}
write( '# Selectivity options for each block' , file=out.file,append=T)   #, ncolumns=(nyears))
write(t(sel.types), file=out.file,append=T, ncolumns=nselblks)
temp = t(sel.mats)
temp = sel.mats
x = asap.nages+6
for(i in 1:nselblks)
{
  write(paste0('# Selectivity Block #', i, " Data") , file=out.file,append=T)   #, ncolumns=(nyears))
  write(t(temp[(i-1)*x + 1:x,]), file=out.file,append=T, ncolumns=4)
}
write( '# Selectivity start age by fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
write(fleet.age1, file=out.file,append=T, ncolumns=nfleets )
write( '# Selectivity end age by fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
write(fleet.age2, file=out.file,append=T, ncolumns=nfleets )
write( '# Age range for average F' , file=out.file, append=T)   #, ncolumns=(nyears))
write(F.report.ages, file=out.file,append=T, ncolumns=2)
write( '# Average F report option ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(F.report.opt, file=out.file,append=T, ncolumns=2)
write( '# Use likelihood constants?' , file=out.file,append=T)   #, ncolumns=(nyears))
write(like.const, file=out.file, append=T )
write( '# Release Mortality by fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
write( rel.mort.fleet, file=out.file,append=T, ncolumns=nfleets)
#write( '# Catch at age matrices (nyears*nfleets rows)' , file=out.file,append=T)   #, ncolumns=(nyears))
write( '# Catch Data', file=out.file,append=T)   #, ncolumns=(nyears))
for(i in 1:nfleets)
{
  write(paste0("# Fleet-", i, " Catch Data"), file=out.file,append=T)
  write(t(caa.mats[(i-1)*nyears + 1:nyears,]), file=out.file,append=T, ncolumns= (asap.nages+1) )
}
write( '# Discards at age by fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
for(i in 1:nfleets)
{
  write(paste0("# Fleet-", i, " Discards Data"), file=out.file,append=T)
  write(t(daa.mats[(i-1)*nyears + 1:nyears,]), file=out.file,append=T, ncolumns= (asap.nages+1) )
}
write( '# Release proportion at age by fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
for(i in 1:nfleets)
{
  write(paste0("# Fleet-", i, " Release Data"), file=out.file,append=T)
  write(t(rel.prop[(i-1)*nyears + 1:nyears,]), file=out.file,append=T, ncolumns= asap.nages )
}
write( '# Survey Index Data' , file=out.file,append=T)   #, ncolumns=(nyears))
write( '# Index units' , file=out.file,append=T)   #, ncolumns=(nyears))
write(units.ind, file=out.file,append=T, ncolumns=n.ind.avail )
write( '# Index Age comp. units' , file=out.file,append=T)   #, ncolumns=(nyears))
write(units.ind, file=out.file,append=T, ncolumns=n.ind.avail )
write( '# Index WAA matrix' , file=out.file,append=T)   #, ncolumns=(nyears))
write((rep(1,n.ind.avail)), file=out.file,append=T, ncolumns=n.ind.avail )

write( '# Index month' , file=out.file, append=T)   #, ncolumns=(nyears))
write(time.ind, file=out.file,append=T, ncolumns=n.ind.avail )
write( '# Index link to fleet? ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(fish.ind, file=out.file,append=T, ncolumns=n.ind.avail)
write( '# Index selectivity option ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(sel.ind, file=out.file,append=T, ncolumns=n.ind.avail)
write( '# Index start age' , file=out.file,append=T)   #, ncolumns=(nyears))
write(ind.age1, file=out.file, append=T, ncolumns=n.ind.avail )
write( '# Index end age' , file=out.file,append=T)   #, ncolumns=(nyears))
write(ind.age2, file=out.file, append=T, ncolumns=n.ind.avail )

write( '# Index Estimate Proportion (YES=1)' , file=out.file,append=T)   #, ncolumns=(nyears))
write(t(rep(1,n.ind.avail)), file=out.file, append=T, ncolumns=n.ind.avail )
write( '# Use Index' , file=out.file,append=T)   #, ncolumns=(nyears))
write(ind.use, file=out.file, append=T, ncolumns=n.ind.avail )
x = asap.nages+6
for(i in 1:n.ind.avail)
{
  write(paste0('# Index-', i, ' Selectivity Data') , file=out.file,append=T)   #, ncolumns=(nyears))
  write(t(ind.sel.mats[(i-1)*x + 1:x,]), file=out.file,append=T, ncolumns=4)
}
write( '# Index data matrices (n.ind.avail.*nyears)' , file=out.file,append=T)   #, ncolumns=(nyears))

# ----------one-off for ICES to ASAP
  for ( kk in 1:length(ind.use))  {
     if (ind.use[kk]==1) {
write( paste0('# Index   ', survey.names[kk]) , file=out.file,append=T)   #, ncolumns=(nyears))
        tmp.s <- ind.mat[[kk]]
        ind.mat2 <- get.index.mat(tmp.s, ind.cv, ind.neff, first.year, nyears, catch.ages, survey.ages[[kk]])
write(t(ind.mat2), file=out.file,append=T, ncolumns=(asap.nages + 4) )
        } # end ind.use test
  } #end kk loop

# ----------one-off for ICES to ASAP
write( '#########################################' , file=out.file,append=T)   #, ncolumns=(nyears))
write( '# Phase data' , file=out.file,append=T)   #, ncolumns=(nyears))
write( '# Phase for Fmult in 1st year' , file=out.file,append=T)   #, ncolumns=(nyears))
write(p.Fmult1, file=out.file,append=T  )
write( '# Phase for Fmult deviations' , file=out.file, append=T)   #, ncolumns=(nyears))
write(p.Fmult.dev, file=out.file,append=T  )
write( '# Phase for recruitment deviations ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(p.recr.dev, file=out.file,append=T )
write( '# Phase for N in 1st year ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(p.N1, file=out.file,append=T )
write( '# Phase for catchability in 1st year' , file=out.file,append=T)   #, ncolumns=(nyears))
write(p.q1, file=out.file, append=T  )
write( '# Phase for catchability deviations' , file=out.file,append=T)   #, ncolumns=(nyears))
write(p.q.dev, file=out.file, append=T )
write( '# Phase for stock recruit relationship' , file=out.file,append=T)   #, ncolumns=(nyears))
write(p.SR, file=out.file, append=T  )
write( '# Phase for steepness' , file=out.file,append=T)   #, ncolumns=(nyears))
write(p.h, file=out.file,append=T  )
write( '#########################################' , file=out.file,append=T)   #, ncolumns=(nyears))
write( '# Lambdas and CVs' , file=out.file,append=T)   #, ncolumns=(nyears))
write( '# Recruitment CV by year' , file=out.file,append=T)   #, ncolumns=(nyears))
write(recr.CV, file=out.file,append=T , ncolumns=1 )
write( '# Lambda for each index' , file=out.file,append=T)   #, ncolumns=(nyears))
write(lam.ind, file=out.file,append=T, ncolumns=n.ind.avail  )
write( '# Lambda for Total catch in weight by fleet' , file=out.file, append=T)   #, ncolumns=(nyears))
write(lam.c.wt, file=out.file,append=T, ncolumns=nfleets  )
write( '# Lambda for total discards at age by fleet ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(lam.disc, file=out.file,append=T, ncolumns=nfleets )
write( '# Catch Total CV by year and fleet ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(catch.CV, file=out.file,append=T, ncolumns=nfleets )
write( '# Discard total CV by year and fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
write(disc.CV, file=out.file, append=T, ncolumns=nfleets  )
write( '# Input effective sample size for catch at age by year and fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
write(Neff.catch, file=out.file, append=T, ncolumns=nfleets )
write( '# Input effective sample size for discards at age by year and fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
write(Neff.disc, file=out.file, append=T , ncolumns=nfleets )
write( '# Lambda for Fmult in first year by fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
write(lam.Fmult.y1, file=out.file,append=T, ncolumns=nfleets  )
write( '# CV for Fmult in first year by fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
write(CV.Fmult.y1, file=out.file,append=T, ncolumns=nfleets  )
write( '# Lambda for Fmult deviations' , file=out.file,append=T)   #, ncolumns=(nyears))
write(lam.Fmult.dev, file=out.file,append=T, ncolumns=nfleets  )
write( '# CV for Fmult deviations' , file=out.file,append=T)   #, ncolumns=(nyears))
write(CV.Fmult.dev, file=out.file,append=T, ncolumns=nfleets  )
write( '# Lambda for N in 1st year deviations ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(lam.N1.dev, file=out.file,append=T )
write( '# CV for N in 1st year deviations ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(CV.N1.dev, file=out.file,append=T  )
write( '# Lambda for recruitment deviations' , file=out.file,append=T)   #, ncolumns=(nyears))
write(lam.recr.dev, file=out.file, append=T  )
write( '# Lambda for catchability in first year by index' , file=out.file,append=T)   #, ncolumns=(nyears))
write(lam.q.y1, file=out.file, append=T, ncolumns=n.ind.avail )
write( '# CV for catchability in first year by index' , file=out.file,append=T)   #, ncolumns=(nyears))
write(CV.q.y1, file=out.file, append=T , ncolumns=n.ind.avail )
write( '# Lambda for catchability deviations by index' , file=out.file,append=T)   #, ncolumns=(nyears))
write(lam.q.dev, file=out.file,append=T, ncolumns=n.ind.avail  )
write( '# CV for catchability deviations by index' , file=out.file,append=T)   #, ncolumns=(nyears))
write(CV.q.dev, file=out.file,append=T  )
write( '# Lambda for deviation from initial steepness' , file=out.file,append=T)   #, ncolumns=(nyears))
write(lam.h, file=out.file,append=T   )
write( '# CV for deviation from initial steepness' , file=out.file,append=T)   #, ncolumns=(nyears))
write(CV.h, file=out.file,append=T  )
write( '# Lambda for deviation from initial SSB0 ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(lam.SSB0, file=out.file,append=T )
write( '# CV for deviation from initial SSB0 ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(CV.SSB0, file=out.file,append=T  )

write( '# NAA Deviations flag (1=   , 0=  ) ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(1, file=out.file,append=T  )

write('###########################################', file=out.file, append=T)
write('###  Initial Guesses', file=out.file, append=T)
write( '# NAA for year1' , file=out.file,append=T)   #, ncolumns=(nyears))
write(naa.y1, file=out.file, append=T, ncolumns=asap.nages  )
write( '# Fmult in 1st year by fleet' , file=out.file,append=T)   #, ncolumns=(nyears))
write(Fmult.y1, file=out.file, append=T, ncolumns=nfleets )
write( '# Catchability in 1st year by index' , file=out.file,append=T)   #, ncolumns=(nyears))
write(q.y1, file=out.file, append=T  )

write( '# S-R Unexploited specification (1=   0=)' , file=out.file,append=T)   #, ncolumns=(nyears))
write(1, file=out.file,append=T, ncolumns=n.ind.avail  )

write( '# Unexploited initial guess' , file=out.file,append=T)   #, ncolumns=(nyears))
write(SSB0, file=out.file,append=T, ncolumns=n.ind.avail  )
write( '# Steepness initial guess' , file=out.file,append=T)   #, ncolumns=(nyears))
write(h.guess, file=out.file,append=T  )
write( '# Maximum F (upper bound on Fmult)' , file=out.file,append=T)   #, ncolumns=(nyears))
write(F.max, file=out.file,append=T  )
write( '# Ignore guesses' , file=out.file,append=T)   #, ncolumns=(nyears))
write(ignore.guess, file=out.file,append=T  )
write('###########################################', file=out.file, append=T)
write('###  Projection Control data', file=out.file, append=T)
write( '# Do projections' , file=out.file,append=T)   #, ncolumns=(nyears))
write(do.proj, file=out.file, append=T   )
write( '# Fleet directed flag' , file=out.file,append=T)   #, ncolumns=(nyears))
write(fleet.dir, file=out.file, append=T, ncolumns=nfleets )
write( '# Final year of projections' , file=out.file,append=T)   #, ncolumns=(nyears))
write(proj.yr, file=out.file, append=T   )
write( '# Year, projected recruits, what projected, target, non-directed Fmult ' , file=out.file,append=T)   #, ncolumns=(nyears))
write(t(proj.specs), file=out.file,append=T, ncolumns=5 )
write('###########################################', file=out.file, append=T)
write('###  MCMC Control data', file=out.file, append=T)
write( '# do mcmc' , file=out.file,append=T)   #, ncolumns=(nyears))
write(do.mcmc, file=out.file,append=T  )
write( '# MCMC nyear option' , file=out.file,append=T)   #, ncolumns=(nyears))
write(mcmc.nyr.opt, file=out.file,append=T  )
write( '# MCMC number of saved iterations desired' , file=out.file,append=T)   #, ncolumns=(nyears))
write(mcmc.nboot, file=out.file,append=T  )
write( '# MCMC thinning rate' , file=out.file,append=T)   #, ncolumns=(nyears))
write(mcmc.thin, file=out.file,append=T  )
write( '# MCMC random number seed' , file=out.file,append=T)   #, ncolumns=(nyears))
write(mcmc.seed, file=out.file,append=T  )
write('###########################################', file=out.file, append=T)
write('###  A few AGEPRO specs', file=out.file, append=T)
write( '# R in agepro.bsn file' , file=out.file,append=T)   #, ncolumns=(nyears))
write(recr.agepro, file=out.file,append=T  )
write( '# Starting year for calculation of R' , file=out.file,append=T)   #, ncolumns=(nyears))
write(recr.start.yr, file=out.file,append=T  )
write( '# Ending year for calculation of R' , file=out.file,append=T)   #, ncolumns=(nyears))
write(recr.end.yr, file=out.file,append=T  )

write( '# Export to R flag (1=  0=)' , file=out.file,append=T)   #, ncolumns=(nyears))
write(1, file=out.file,append=T  )

write( '# test value' , file=out.file,append=T)   #, ncolumns=(nyears))
write(test.val, file=out.file,append=T  )
write('###########################################', file=out.file, append=T)
write('###### FINIS ######', file=out.file, append=T)
write( '# Fleet Names', file=out.file, append=T)
write(fleet.names, file=out.file, append=T, ncolumns=1)
write( '# Survey Names', file=out.file, append=T)
write(survey.names, file=out.file, append=T, ncolumns=1)



 }     # end asap setup function
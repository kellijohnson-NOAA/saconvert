#' 
#' 
get.index.mat<- function(x, cv, neff, first.year, nyears, catch.ages, survey.ages)  {
   n.ages = length(catch.ages)
   last.yr <- first.year+nyears - 1
   
   tmp.yrs <- as.numeric(rownames(x))
   all.years = first.year-1 + 1:nyears
   years.use.ind = which(tmp.yrs %in% all.years)
   #if (tmp.yrs[length(tmp.yrs)]>last.yr)  tmp.yrs <- tmp.yrs[-which(tmp.yrs>last.yr)]
   tmp.ages <- as.numeric(colnames(x))
   tmp.ages = catch.ages
   survey.ages.index = which(catch.ages %in% survey.ages)
   i.mat <- matrix(0, nyears, (n.ages + 4))
   i.mat[,1] <- all.years
   rownames(x) <- c()
   colnames(x) <- c()
   x[is.na(x)] <- 0
   print(dim(x))
   print(dim(i.mat))
   print(survey.ages.index)
   print(tmp.yrs)
   print(sum(all.years %in% tmp.yrs))
   print(all.years)
   print(x)
   tmp.ind.total <- apply(x[years.use.ind,], 1, sum)
   i.mat.ind = which(all.years %in% tmp.yrs[years.use.ind])
   i.mat[i.mat.ind,2:3] <- cbind(tmp.ind.total, rep(cv, length(years.use.ind)))
   i.mat[i.mat.ind, (3+survey.ages.index)]  <- x[years.use.ind,]
   i.mat[i.mat.ind, (n.ages+4)]  <- rep(neff, length(years.use.ind))

 return(i.mat)

}

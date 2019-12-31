setup.surv.sel = function(x, i.peak, catch.ages,survey.ages)  {   #set up matrix specs for index selectivities
 n.ind <- length(x)
 sel.c1 <- rep()
 sel.c2 <- rep()
 sel.c3 <- rep()
 sel.c4 <- rep()
 n.ages = length(catch.ages)
  for (i in 1:n.ind) {
#  tmp.nages <- ind.ages[2,i]-ind.ages[1,i]+1
   tmp.c1 = rep(0,n.ages)
   if(sum(!(survey.ages[[i]] %in% catch.ages))) stop("some survey ages are not in catch.ages")
   ind = which(catch.ages %in% survey.ages[[i]])
   peak.age.class = ind[i.peak[i]] #not necessarily the peak age
   tmp.c1[ind] = seq(0.1,0.9, length.out=length(ind))
   tmp.c1[peak.age.class] <-1
   tmp.c1 <-  c( tmp.c1,  round((peak.age.class)/2,2), 0.9,
              round((peak.age.class)/4,2), 0.6,  round(n.ages/1.5,2), 1.1)
   sel.c1 <- c(sel.c1, tmp.c1)
   tmp.c2 = rep(-1,n.ages)
   tmp.c2[ind] = 1
   tmp.c2[peak.age.class] = -1
   tmp.c2 = c(tmp.c2, 2,3, rep(1,4))
   sel.c2 <- c( sel.c2, tmp.c2)#, 2,3, rep(1, 4)) # phase for estimation
   #sel.c2[peak.age.class+(i-1)*(n.ages+6)] <- -1
   sel.c3 <- c(sel.c3, rep(0, (n.ages+6))  )# lambda for sel parameters
   sel.c4 <- c(sel.c4, rep(1, (n.ages+6)) )# CV for sel parameters (irrelevant if lambda=0)

  }#end i loop

 sel.mats <- cbind(sel.c1, sel.c2, sel.c3, sel.c4)

 return(sel.mats)

}

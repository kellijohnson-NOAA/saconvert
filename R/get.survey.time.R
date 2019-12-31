get.survey.time=function (x)    #func to grab survey timing from ICES surveys object
{
n.surv <- length(x)
tt <- rep(NA, n.surv)
for (i in 1:n.surv) {
    tt[i] <- attr(x[[i]], 'time') [1]
}
    tt <- round(tt*12, 0)
return(tt)

 }

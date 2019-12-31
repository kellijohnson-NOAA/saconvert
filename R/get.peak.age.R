get.peak.age=function (x)    #grab peak age from first couple of rows
{
if (class(x)== "matrix") {
   t.age <- rep(NA, 5)
   for (i in 1:5) {
     if(all(is.na(x[i,]))) t.age[i] = NA
     else t.age[i] <- which(x[i,]==max(x[i,], na.rm=T))
   }
   peak <- round(mean(t.age, na.rm = TRUE), 0)
    } #end matrix class

if (class(x)== "list") {
 n.mats <- length(x)
 peak <-rep(NA, n.mats)

 for (i in 1:n.mats) {
    t.mat <- x[[i]] [1:5,]
    t.age <- rep(NA,5)
    for (j in 1:5) {
     if(all(is.na(t.mat[j,]))) t.age[j] = NA
     else t.age[j] <- which(t.mat[j,]==max(t.mat[j,], na.rm=T))
     #t.age[j] <- which(t.mat[j,]==max(t.mat[j,], na.rm=T))
    } # end j loop
    peak[i] <- round(mean(t.age, na.rm = TRUE), 0)
        }  # end i loop

    } # end list class

return(peak)

 }

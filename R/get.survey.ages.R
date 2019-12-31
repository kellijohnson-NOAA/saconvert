get.survey.ages=function (x)    #func to grab survey ages from ICES surveys object
{
  lapply(x, function(y) 
  {
    r = range(as.numeric(colnames(y)))
    seq(r[1],r[2],1)
  })
}

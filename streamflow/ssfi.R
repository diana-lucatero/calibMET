## Function that computes ssfi

ssfi <- function(obs,nyears){
  ## Compute start and end of months
  days <- days_in_month(seq(from=as.Date("1990/1/1"),to=as.Date("2014/12/1"), by = "month"))
  days <- matrix(days,nrow = nyears,ncol = 12,byrow = T)
  s1 <- t(apply(days,1, cumsum)) # end day
  s0 <- s1 - (days-1)     # start day
  ## Computes percentage of dry days
  accu = array(NaN,dim=c(nyears,12))
  for (im in 1:12){
    accu0 = sapply(1:nyears, function(iy) mean(obs[iy,s0[iy,im]:s1[iy,im]]))
    accu[,im] = accu0
  }
  h <- matrix(t(accu),nrow = 25*12,ncol = 1)
  fit<-fitSCI(h,1,1,'weibull',p0=F)
  ssfi <- transformSCI(h,first.mon=1,obj=fit)
  return(ssfi) ## A nyearxnmonth matrix with percentages
}
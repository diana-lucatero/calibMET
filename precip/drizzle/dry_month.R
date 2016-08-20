## Function that computes number of dry day per month

dry_month <- function(obs,nyears){
  ## Compute start and end of months
  days <- days_in_month(seq(from=as.Date("1990/1/1"),to=as.Date("2014/12/1"), by = "month"))
  days <- matrix(days,nrow = nyears,ncol = 12,byrow = T)
  s1 <- t(apply(days,1, cumsum)) # end day
  s0 <- s1 - (days-1)     # start day
  ## Computes percentage of dry days
  dry_p = array(NaN,dim=c(nyears,12))
  for (im in 1:12){
    dry_d = sapply(1:nyears, function(iy) sum(obs[iy,s0[iy,im]:s1[iy,im]]==0)/days[iy,im])
    dry_p[,im] = dry_d*100
  }
  return(dry_p) ## A nyearxnmonth matrix with percentages
}
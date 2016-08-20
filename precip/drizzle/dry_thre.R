## Function that computes thresholds at which the % of forecasted
## dry days matches the % of the observed dry days on a given month

dry_thre <- function(ens,nmonth,dry,ngrid,no_ens,nyears){
  s0 <- c(1,31,61,91,121,151,181)
  s1 <- c(30,60,90,120,150,180,215)
  d  <- c(rep(30,6),35)
  date0 <- paste0('1990/',nmonth,'/01')
  m <- month(seq(from=as.Date(date0),by='months',length = 7))
  dryP <- dry[ngrid,,m]
  thre <- array(NaN,dim = c(no_ens,nyears,7))
  for (iy in 0:(nyears-1)){
    for (im in 1:7){
      ens0 <- ens[((iy*no_ens)+1):((iy*no_ens)+no_ens),s0[im]:s1[im]]
      ensR <- t(sapply(1:no_ens, function(ie) sort(ens0[ie,],na.last = NA)))
      per0 <- dryP[(iy+1),im];per <- round(d[im]*(per0/100))
      if (per == 0) {thre[,(iy+1),im] <- array(NaN,dim = c(no_ens))}
      else {thre[,(iy+1),im] <- ensR[,per]}
    } # end month
  } # end year
  return(thre)
} ## End function
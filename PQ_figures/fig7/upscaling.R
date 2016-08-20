## Temporal and spatial upscaling

Dir <- 'E:/calibMET/'
km <- c(10,20,40)
nsqr1 <- c(662,70,24)
nsqr2 <- c(724,75,31)

month <- 5
#pcp1 <- array(NaN,dim = c(662,4,length(km)))
#temp1 <- array(NaN,dim = c(724,4,length(km)))
#evap1 <- array(NaN,dim = c(724,4,length(km)))

pcp1 <- array(NaN,dim = c(length(km),4))
temp1 <- array(NaN,dim = c(length(km),4))
evap1 <- array(NaN,dim = c(length(km),4))

for (ik in 1:length(km)){
  ## PRECIPITATION
  name_file <- paste0(Dir,'precip/verification_daily/crpss_',km[ik],'_m',month,'.txt')
  p0 <- t(matrix(scan(name_file),nrow = 4,ncol = 662))
  #pcp1[1:nsqr1[ik],,ik] <- read.table(name_file)
  pcp1[ik,] <- colMeans(read.table(name_file))
  ## TEMPERATURE
  name_file <- paste0(Dir,'temp/verification_daily/crpss_',km[ik],'_m',month,'.txt')
  temp1[ik,] <- colMeans(read.table(name_file))
  ## EVAPOTRANSPIRATION
  name_file <- paste0(Dir,'ref_evap/verification_daily/crpss_',km[ik],'_m',month,'.txt')
  evap1[ik,] <- colMeans(read.table(name_file))
}

month <- 11
#pcp2 <- array(NaN,dim = c(662,4,length(km)))
#temp2 <- array(NaN,dim = c(724,4,length(km)))
#evap2 <- array(NaN,dim = c(724,4,length(km)))

pcp2 <- array(NaN,dim = c(length(km),4))
temp2 <- array(NaN,dim = c(length(km),4))
evap2 <- array(NaN,dim = c(length(km),4))

for (ik in 1:length(km)){
  ## PRECIPITATION
  name_file <- paste0(Dir,'precip/verification_daily/crpss_',km[ik],'_m',month,'.txt')
  pcp2[ik,] <- colMeans(read.table(name_file))
  ## TEMPERATURE
  name_file <- paste0(Dir,'temp/verification_daily/crpss_',km[ik],'_m',month,'.txt')
  temp2[ik,] <- colMeans(read.table(name_file))
  ## EVAPOTRANSPIRATION
  name_file <- paste0(Dir,'ref_evap/verification_daily/crpss_',km[ik],'_m',month,'.txt')
  evap2[ik,] <- colMeans(read.table(name_file))
}


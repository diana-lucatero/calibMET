## Standarized Streamflow Index (SSFI)

library(lubridate);library(SCI)

Dir <- 'E:/calibMET/streamflow/'

source0 = paste0(Dir,'/ssfi.r')
source(source0)

nyears <- 25
station <- c(20,21,82)
name_file <- paste0(Dir,'obs_station_',station[1],'.txt')
obs <- t(matrix(scan(name_file),nrow = 366,ncol = 25))

ssfi <- ssfi(obs,nyears)

Date <- seq(from=as.Date("1990/1/1"),to=as.Date("2014/12/1"), by = "month")
plot(Date,ssfi,type = 'l')
abline(h=-1,col='red')
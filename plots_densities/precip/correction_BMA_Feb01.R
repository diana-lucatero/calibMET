## Correction BMA

ptm = proc.time()
Nlead = 215 #Lead time in days
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
library(ensembleBMA)

ObsDir = 'E:/calibMET/precip/'
grid_filename = paste(ObsDir,'obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
info_filename = paste(ObsDir,'ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))

nmonth = 2

for (ngrid in 269:269){ # Loop for space :dim(grid)[1]
  print(ngrid)
  obs_filename = paste(ObsDir,'observation/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_filename),nrow=216,ncol=293)); obs[obs<0]=0;
  ens_filename = paste(ObsDir,'ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322)); ens[ens<0]=0;
  index     = which(ens_info[1:293,3] == nmonth)
  index_out = which(ens_info[index,4] != Nens[nmonth])
  index     = index[-(index_out)]
  ens_index = ens_info[index,c(31,32)]
  no_ens    = Nens[nmonth]
  dimIndex  = length(index)
  no_ens_clim = dimIndex[1]-1
  index_f   = matrix(nrow=1, ncol=dimIndex*no_ens)
  index_f   = sapply(1:dimIndex, function(ix) c(ens_index[ix,1]:ens_index[ix,2]))
  # Extract ensembles with start date: nmonth [Jan(1) to Dec(12)]
  Pens      = ens[index_f,3:217]
  # Extract observations with start date: nmonth [Jan(1) to Dec(12)]
  Pobs = obs[index,2:216]
  ## Start of Postprocessing     
  #Pbma  = matrix(nrow=dimIndex[1]*no_ens, ncol=Nlead)
  exc = rep(1,no_ens)
  for (nyear in 12:12){  #  0:(dimIndex[1]-1) Loop for year of validation
    index_nyear = ((nyear*no_ens)+1):((nyear*no_ens)+no_ens)
    thres = seq(from=0,to=100,by=1)
    probsample=runif(no_ens)
    for (nlead in 30:30){#1:Nlead
      Xpredictor  = (Pens[-(index_nyear),nlead])
      Xpredict    = (Pens[index_nyear,nlead])
      Ypredictand = Pobs[-(nyear+1),nlead]
      train = t(matrix(Xpredictor,nrow = no_ens,ncol=dimIndex[1]-1))
      ens_data = ensembleData(forecasts=train,observations=Ypredictand,verbose=FALSE,dates=NULL,forecastHour = NULL,initializationTime = NULL)
      calib = fitBMAgamma0(ens_data,exchangeable = exc)
      fore = ensembleData(t(Xpredict),verbose=FALSE,dates=NULL,forecastHour = NULL,initializationTime = NULL)
#     #imax = which(fore==max(fore)); mmax = calib$biasCoefs[1,imax]+calib$biasCoefs[2,imax]*(fore[imax])^(1/3)
#     #if (mmax < 0) fore[imax] = fore[imax]^(1/3)
      Probcurve = 1 - cdf.fitBMAgamma0(calib, fore, values = thres)
      #Pblm = t(approx(Probcurve, thres, probsample, rule=2, method="linear")$y)
      Pbma = t(approx(Probcurve, thres, probsample, rule=2, method="linear")$y)
      #Pbma[((nyear*no_ens)+1):((nyear*no_ens)+no_ens),nlead] = t(approx(Probcurve, thres, probsample, rule=2, method="linear")$y)
    } # End lead time
  } # End cross-validation
  ## Save ensembles 
  #name_file = paste(ObsDir,'ensembles_c/ens_bma_g',grid[ngrid,2],'_m',nmonth,'.txt',sep='')
  #write.table(Pbma, file=name_file, row.names=FALSE, col.names=FALSE)
  #Deletes files for grid ngrid
  #rm(obs) ;  rm(ens) ; rm(Pens) ; rm(Pobs)
} # End loop grid

proc.time() - ptm
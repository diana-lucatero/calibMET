## ## Correction Empirical Quantile Mapping Precipitation
#### PART 1. Extract info

Nlead = 215 #Lead time in days
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
library(Hmisc)#;library(leaps);library(locfit);library(relaimpo)
library(matrixStats);library(boot)

ObsDir = 'E:/calibMET/precip/'
grid_filename = paste(ObsDir,'obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
info_filename = paste(ObsDir,'ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))

## Compute corrected monthly average
for (ngrid in 269:269){ # Loop for space :dim(grid)[1]
  print(ngrid)
  obs_filename = paste(ObsDir,'observation/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_filename),nrow=216,ncol=293));obs[obs<0]=0;
  ens_filename = paste(ObsDir,'ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322)); ens[ens<0]=0;
  for (nmonth in 2:2){ # Loop for forecast start month
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
    Pqm       = matrix(nrow=dimIndex[1]*no_ens, ncol=Nlead)
    for (nyear in 12:12){  # 0:(dimIndex[1]-1) for 2003 02 
      index_nyear = ((nyear*no_ens)+1):((nyear*no_ens)+no_ens)
      ## Second method - Quantile mapping
      for (nlead in 30:30){#1:Nlead
        train = Pens[-(index_nyear),nlead]
        calib = Pens[index_nyear,nlead]
        Pcalib = array(NaN,dim=c(no_ens))
        for (iens in 1:no_ens){ # Start ensemble loop
          iiens = seq(iens,length(train),no_ens)
          modt = sort(train[iiens],decreasing=TRUE); Nexcm = length(modt)
          Pexcm = (1:Nexcm)/(Nexcm+1) # empirical ensembles
          percentile = approxExtrap(modt,Pexcm,calib[iens],method="linear")
          percentile$y[percentile$y < 0] = 0.0001; percentile$y[percentile$y > 1] = 0.9999 # Correct in case extrapolations are too ugly
          Pexc = (1:dimIndex[1])/(dimIndex[1]+1) # empirical observations
          Pcalib[iens] = approxExtrap(Pexc,sort(Pobs[-(nyear+1),nlead],decreasing = TRUE),
                                                                             percentile$y[!is.na(percentile$y)],
                                                                             method = "linear")$y
        } # End ensemble loop
        #Pqm[((nyear*no_ens)+1):((nyear*no_ens)+no_ens),nlead] = Pcalib
      } # End lead time loop in qm
    } # Finish loop validation years
  } # End loop month
  #Deletes files for grid ngrid
  #rm(obs) ;  rm(ens)
} # End loop grid
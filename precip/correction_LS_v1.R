## Correction Linear Scaling Precipitation CF = obs/ens
#### PART 1. Extract info

Nlead = 215 #Lead time in days
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
library(Hmisc);library(leaps);library(locfit)
library(matrixStats);library(boot)

ObsDir = 'E:/calibMET/precip/'
grid_filename = paste(ObsDir,'obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
info_filename = paste(ObsDir,'ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))

## Compute corrected monthly agregated
for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  obs_filename = paste(ObsDir,'observation/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_filename),nrow=216,ncol=293));obs[obs<0]=0;
  ens_filename = paste(ObsDir,'ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322)); ens[ens<0]=0;
  for (nmonth in 9:12){ # Loop for forecast start month
    index     = which(ens_info[1:293,3] == nmonth)
    #index_out = which(ens_info[index,4] != Nens[nmonth])
    #index     = index[-(index_out)]
    ens_index = ens_info[index,c(31,32)];no_ens    = Nens[nmonth]
    dimIndex  = length(index)
    index_f   = matrix(nrow=1, ncol=dimIndex*no_ens)
    index_f   = sapply(1:dimIndex, function(ix) c(ens_index[ix,1]:ens_index[ix,2]))
    # Extract ensembles with start date: nmonth [Jan(1) to Dec(12)]
    Pens      = ens[index_f,3:217]
    # Extract observations with start date: nmonth [Jan(1) to Dec(12)]
    Pobs = obs[index,2:216]
    ## Start of Postprocessing     
    factor = matrix(nrow=dimIndex[1], ncol=Nlead)
    Pdc       = matrix(nrow=dimIndex[1]*no_ens, ncol=Nlead)
    for (nyear in 0:(dimIndex[1]-1)){  # Loop for year of validation
      index_nyear = ((nyear*no_ens)+1):((nyear*no_ens)+no_ens)
      mens = colMeans(Pens[-(index_nyear),],na.rm = TRUE) # Mean of all = Mean of the mean
      mobs = colMeans(Pobs[-(nyear+1),],na.rm = TRUE)
      ## First method - Bias correction linear scaling
      BC_factor = mobs/mens 
      factor[(nyear+1),] = BC_factor # Save factor
      Pdc[((nyear*no_ens)+1):((nyear*no_ens)+no_ens),] = t(BC_factor * t(Pens[index_nyear,]))
    } # Finish loop validation years
    ## Save ensembles plus correction factor
    name_file = paste(ObsDir,'ensembles_c/ens_dc_g',grid[ngrid,2],'_m',nmonth,'.txt',sep='')
    write.table(Pdc, file=name_file, row.names=FALSE, col.names=FALSE)
    name_file = paste(ObsDir,'ensembles_c/factor_g',grid[ngrid,2],'_m',nmonth,'.txt',sep='')
    write.table(factor, file=name_file, row.names=FALSE, col.names=FALSE)
    ## Deletes ensembles
    rm(Pens);rm(Pobs)
    rm(factor);rm(Pdc)
  } # End loop month
  #Deletes files for grid ngrid
  rm(obs) ;  rm(ens)
} # End loop grid
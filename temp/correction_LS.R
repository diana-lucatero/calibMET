## Correction Linear Scaling Temperature CF = obs- ens
#### PART 1. Extract info

ptm = proc.time()
Nlead = 215 #Lead time in days
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)

ObsDir = 'E:/calibMET/temp/'
grid_filename = paste(ObsDir,'obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=5,ncol=724))
info_filename = paste(ObsDir,'ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))


s0 <- c(1,31,61,91,121,151,181)
s1 <- c(30,60,90,120,150,180,215)
d  <- c(rep(30,6),35)

## Compute corrected monthly agregated
for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  obs_filename = paste(ObsDir,'observations/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_filename),nrow=216,ncol=293))
  ens_filename = paste(ObsDir,'ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322))
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
    factor = matrix(nrow=dimIndex[1], ncol=7)
    Pdc       = matrix(nrow=dimIndex[1]*no_ens, ncol=Nlead)
    for (nyear in 0:(dimIndex[1]-1)){  # Loop for year of validation
      index_nyear = ((nyear*no_ens)+1):((nyear*no_ens)+no_ens)
      # Train
      for (nt in 1:7){#7
        Xpredictor  = Pens[-(index_nyear),s0[nt]:s1[nt]]
        Ypredictand = Pobs[-(nyear+1),s0[nt]:s1[nt]]
        mens = mean(Xpredictor,na.rm = TRUE) # Mean of all = Mean of the mean
        mobs = mean(Ypredictand,na.rm = TRUE)
        ## First method - Bias correction linear scaling
        BC_factor = mobs - mens 
        factor[(nyear+1),nt] = BC_factor # Save factor
        for (nlead in s0[nt]:s1[nt]){
          Pdc[((nyear*no_ens)+1):((nyear*no_ens)+no_ens),nlead] = t(BC_factor + t(Pens[index_nyear,nlead]))
        }# End lead time
      }# End month
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

proc.time() - ptm
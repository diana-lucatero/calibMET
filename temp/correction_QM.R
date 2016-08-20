## ## Correction Empirical Quantile Mapping
#### PART 1. Extract info

ptm = proc.time()
Nlead = 215 #Lead time in days
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)

library(Hmisc)

ObsDir = 'E:/calibMET/temp/'
grid_filename = paste(ObsDir,'obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=5,ncol=724))
info_filename = paste(ObsDir,'ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))

s0 <- c(1,31,61,91,121,151,181)
s1 <- c(30,60,90,120,150,180,215)
d  <- c(rep(30,6),35)

## Compute corrected monthly average
for (ngrid in 1:1){ # Loop for space :dim(grid)[1]
  print(ngrid)
  obs_filename = paste(ObsDir,'observations/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_filename),nrow=216,ncol=293))
  ens_filename = paste(ObsDir,'ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322))
  for (nmonth in 1:1){ # Loop for forecast start month
    index     = which(ens_info[1:293,3] == nmonth)
    #index_out = which(ens_info[index,4] != Nens[nmonth])
    #index     = index[-(index_out)]
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
    for (nyear in 0:(dimIndex[1]-1)){  # Loop for year of validation 0:(dimIndex[1]-1)
      index_nyear = ((nyear*no_ens)+1):((nyear*no_ens)+no_ens)
      ## Second method - Quantile mapping
      # Train
      for (nt in 1:7){#7
        Xpredictor  = Pens[-(index_nyear),s0[nt]:s1[nt]]
        Xpredictor0 <- matrix(Xpredictor,ncol = 1,nrow = dim(Xpredictor)[1]*dim(Xpredictor)[2])
        Ypredictand = Pobs[-(nyear+1),s0[nt]:s1[nt]]
        Ypredictand0 <- matrix(Ypredictand,ncol = 1,nrow = dim(Ypredictand)[1]*dim(Ypredictand)[2])
        for (nlead in s0[nt]:s1[nt]){#s0[nt]:s1[nt]
          calib = Pens[index_nyear,nlead]
          Pcalib = array(NaN,dim=c(no_ens))
          for (iens in 1:no_ens){ # Start ensemble loop no_ens
            iiens = seq(iens,length(Xpredictor0),no_ens)
            modt = sort(Xpredictor0[iiens],decreasing=TRUE); Nexcm = length(modt)
            Pexcm = (1:Nexcm)/(Nexcm+1) # empirical ensembles
            percentile = approxExtrap(modt,Pexcm,calib[iens],method="linear")
            percentile$y[percentile$y < 0] = 0.0001; percentile$y[percentile$y > 1] = 0.9999 # Correct in case extrapolations are too ugly
            modt0 <- sort(Ypredictand0,decreasing = TRUE); Nexcm0 = length(modt0)
            Pexc = (1:Nexcm0)/(Nexcm0+1) # empirical observations
            Pcalib[iens] = approxExtrap(Pexc,modt0,
                                      percentile$y[!is.na(percentile$y)],
                                      method = "linear")$y
          } # End ensemble loop
          Pqm[((nyear*no_ens)+1):((nyear*no_ens)+no_ens),nlead] = Pcalib
        } # End lead time loop in qm
      }# End month
    } # Finish loop validation years
    ## Save ensembles plus 
    name_file = paste(ObsDir,'ensembles_c/ens_qm_g',grid[ngrid,2],'_m',nmonth,'.txt',sep='')
    write.table(Pqm, file=name_file, row.names=FALSE, col.names=FALSE)
    ## Deletes ensembles
    rm(Pens);rm(Pobs);rm(Pqm)
  } # End loop month
  #Deletes files for grid ngrid
  rm(obs) ;  rm(ens)
} # End loop grid

proc.time() - ptm
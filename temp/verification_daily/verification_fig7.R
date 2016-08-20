## Computes verification of upscaling

##Computes verification precipitation daily

ptm = proc.time()

ObsDir = 'E:/calibMET/temp'
month = seq(1,12,1);mlead = seq(1,7,1)
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
Nup = 4 # 1 <- daily, 2 <- weekly, 3 <- monthly, 4 <- season
nop <- 1 # 1 <- mean, 2 <- sum
library(SpecsVerification);library(Metrics)
library(matrixStats)
source(paste0(ObsDir,'/verification_daily/upscale.R'))

s0 <- c(31,61,91)
s1 <- c(60,90,120)
d  <- c(rep(30,3))

grid_filename = paste(ObsDir,'/obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=5,ncol=724))
info_filename = paste(ObsDir,'/ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))

nmonth <- 5
index     = which(ens_info[1:293,3] == month[nmonth])
#index_out = which(ens_info[index,4] != Nens[nmonth]) ## For Feb and Aug
#index     = index[-(index_out)]
ens_index = ens_info[index,c(31,32)]
no_ens    = Nens[nmonth];dimIndex  = length(index)
no_ens_clim = dimIndex[1]-1
index_f   = matrix(nrow=1, ncol=dimIndex*no_ens)
index_f   = sapply(1:dimIndex, function(ix) c(ens_index[ix,1]:ens_index[ix,2]))
Pclim = array(NaN, dim=c(dimIndex[1],215))
# Verification arrays raw
crpssr_clim = array(NaN,dim=c(dim(grid)[1],Nup))

## Upscale to 10 km
for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  ## Observations
  obs_namefile = paste(ObsDir,'/observations/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_namefile),nrow=216,ncol=293))
  obs = obs[index,2:216]
  ## Climatology
  for (nyear in 1:dimIndex[1]){
    Pclim[nyear,] = colMeans(obs[-(nyear),1:215])
  }
  ## Ensembles Raw
  ens_filename0 = paste(ObsDir,'/ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens0 = t(matrix(scan(ens_filename0),nrow=217,ncol=8322))
  Pens00 = ens0[index_f,3:217];rm(ens0)
  ## Verification raw
  Pens0 = Pens00[,s0[1]:s1[3]]
  gobs = obs[,s0[1]:s1[3]]
  gclim = Pclim[,s0[1]:s1[3]]
  ## Daily
  Pens <- array(Pens0,dim = dim(Pens0)[1]*dim(Pens0)[2])
  gobs0 <- matrix(gobs,nrow  = dimIndex*dim(Pens0)[2],ncol = 1)
  gclim0 <- matrix(gclim,nrow = dimIndex*dim(Pens0)[2],ncol = 1)
  ens_vr = sapply(0:(dimIndex*dim(Pens0)[2]-1), function(nind) Pens[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
  MAE_clim = mae(gobs0,gclim0)
  crpsr = round(mean(EnsCrps(ens_vr,gobs0)),digits=2)
  crpssr_clim[ngrid,1] = round(1-crpsr/MAE_clim,digits = 2)
  for (iup in 2:3){#Nup
    Pens1 <- upscale(Pens0,nop,iup)
    gobs1 <- upscale(gobs,nop,iup)
    gclim1 <- upscale(gclim,nop,iup)
    Pens <- array(Pens1,dim = dim(Pens1)[1]*dim(Pens1)[2])
    gobs1 <- matrix(gobs1,nrow  = dimIndex*dim(Pens1)[2],ncol = 1)
    gclim1 <- matrix(gclim1,nrow = dimIndex*dim(Pens1)[2],ncol = 1)
    ens_vr = sapply(0:(dimIndex*dim(Pens1)[2]-1), function(nind) Pens[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
    MAE_clim = mae(gobs1,gclim1)
    crpsr = round(mean(EnsCrps(ens_vr,gobs1)),digits=2)
    crpssr_clim[ngrid,iup] = round(1-crpsr/MAE_clim,digits = 2)
  }# End lead time raw
  ## Seasonal
  Pens1 <- upscale(Pens0,nop,4)
  gobs1 <- upscale(gobs,nop,4)
  gclim1 <- upscale(gclim,nop,4)
  ens_vr = sapply(0:(dimIndex-1), function(nind) Pens1[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
  MAE_clim = mae(gobs1,gclim1)
  crpsr = round(mean(EnsCrps(ens_vr,gobs1)),digits=2)
  crpssr_clim[ngrid,4] = round(1-crpsr/MAE_clim,digits = 2)
} #End grid
#Save verification 662(ngrid) x 4 upscaling scales
write(t(crpssr_clim),file=paste0(ObsDir,'/verification_daily/crpss_10_m',nmonth,'.txt'),ncolumns=4)




## Upscale to 20 km
name_file <- paste0(ObsDir,'/verification_daily/temp_20.txt')
squares <- read.table(name_file)
# Read grid points in squares
Pclim = array(NaN, dim=c(dimIndex[1],215))
# Verification arrays raw
crpssr_clim = array(NaN,dim=c(75,Nup))

for (nsqr in 1:75){ # Loop for space :75
  print(nsqr)
  st0 <- (nsqr-1)*3+1; st1 <- (nsqr-1)*3+3
  igrd <- squares[st0:st1,];Ngrid <- sum(!is.na(igrd))
  grid <- which(!is.na(igrd),arr.ind = T)
  Obs<- array(NaN,dim = c(293,215,Ngrid));Ens <- array(NaN,dim = c(8322,215,Ngrid))
  for (ngrid in 1:Ngrid){
    ## Observations
    obs_namefile = paste(ObsDir,'/observations/obs_station_',igrd[grid[ngrid,1],grid[ngrid,2]],'.txt',sep='')
    obs = t(matrix(scan(obs_namefile),nrow=216,ncol=293));obs[obs<0]=0;
    Obs[,,ngrid] = obs[,2:216];rm(obs)
    ## Ensembles Raw
    ens_filename0 = paste(ObsDir,'/ensembles_idw/ens_station_',igrd[grid[ngrid,1],grid[ngrid,2]],'.txt',sep='')
    ens0 = t(matrix(scan(ens_filename0),nrow=217,ncol=8322));ens0[ens0<0]=0;
    Ens[,,ngrid] = ens0[,3:217];rm(ens0)
  }# End square
  # Average square
  obs <- apply(Obs,c(1,2),mean);rm(Obs);obs <- obs[index,]
  Pens00 <- apply(Ens,c(1,2),mean);rm(Ens); Pens00 <- Pens00[index_f,]
  ## Climatology
  for (nyear in 1:dimIndex[1]){
    Pclim[nyear,] = colMeans(obs[-(nyear),1:215])
  }
  ## Verification raw
  Pens0 = Pens00[,s0[1]:s1[3]]
  gobs = obs[,s0[1]:s1[3]]
  gclim = Pclim[,s0[1]:s1[3]]
  ## Daily
  Pens <- array(Pens0,dim = dim(Pens0)[1]*dim(Pens0)[2])
  gobs0 <- matrix(gobs,nrow  = dimIndex*dim(Pens0)[2],ncol = 1)
  gclim0 <- matrix(gclim,nrow = dimIndex*dim(Pens0)[2],ncol = 1)
  ens_vr = sapply(0:(dimIndex*dim(Pens0)[2]-1), function(nind) Pens[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
  MAE_clim = mae(gobs0,gclim0)
  crpsr = round(mean(EnsCrps(ens_vr,gobs0)),digits=2)
  crpssr_clim[nsqr,1] = round(1-crpsr/MAE_clim,digits = 2)
  for (iup in 2:3){#Nup
    Pens1 <- upscale(Pens0,nop,iup)
    gobs1 <- upscale(gobs,nop,iup)
    gclim1 <- upscale(gclim,nop,iup)
    Pens <- array(Pens1,dim = dim(Pens1)[1]*dim(Pens1)[2])
    gobs1 <- matrix(gobs1,nrow  = dimIndex*dim(Pens1)[2],ncol = 1)
    gclim1 <- matrix(gclim1,nrow = dimIndex*dim(Pens1)[2],ncol = 1)
    ens_vr = sapply(0:(dimIndex*dim(Pens1)[2]-1), function(nind) Pens[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
    MAE_clim = mae(gobs1,gclim1)
    crpsr = round(mean(EnsCrps(ens_vr,gobs1)),digits=2)
    crpssr_clim[nsqr,iup] = round(1-crpsr/MAE_clim,digits = 2)
  }# End lead time raw
  ## Seasonal
  Pens1 <- upscale(Pens0,nop,4)
  gobs1 <- upscale(gobs,nop,4)
  gclim1 <- upscale(gclim,nop,4)
  ens_vr = sapply(0:(dimIndex-1), function(nind) Pens1[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
  MAE_clim = mae(gobs1,gclim1)
  crpsr = round(mean(EnsCrps(ens_vr,gobs1)),digits=2)
  crpssr_clim[nsqr,4] = round(1-crpsr/MAE_clim,digits = 2)
} #End grid

#Save verification 70(nsqr) x 4 upscaling scales
write(t(crpssr_clim),file=paste0(ObsDir,'/verification_daily/crpss_20_m',nmonth,'.txt'),ncolumns=4)



## Upscale to 40 km
name_file <- paste0(ObsDir,'/verification_daily/temp_40.txt')
squares <- read.table(name_file)
# Read grid points in squares
Pclim = array(NaN, dim=c(dimIndex[1],215))
# Verification arrays raw
crpssr_clim = array(NaN,dim=c(31,Nup))

for (nsqr in 1:31){ # Loop for space :31
  print(nsqr)
  st0 <- (nsqr-1)*5+1; st1 <- (nsqr-1)*5+5
  igrd <- squares[st0:st1,];Ngrid <- sum(!is.na(igrd))
  grid <- which(!is.na(igrd),arr.ind = T)
  Obs<- array(NaN,dim = c(293,215,Ngrid));Ens <- array(NaN,dim = c(8322,215,Ngrid))
  for (ngrid in 1:Ngrid){
    ## Observations
    obs_namefile = paste(ObsDir,'/observations/obs_station_',igrd[grid[ngrid,1],grid[ngrid,2]],'.txt',sep='')
    obs = t(matrix(scan(obs_namefile),nrow=216,ncol=293));obs[obs<0]=0;
    Obs[,,ngrid] = obs[,2:216];rm(obs)
    ## Ensembles Raw
    ens_filename0 = paste(ObsDir,'/ensembles_idw/ens_station_',igrd[grid[ngrid,1],grid[ngrid,2]],'.txt',sep='')
    ens0 = t(matrix(scan(ens_filename0),nrow=217,ncol=8322));ens0[ens0<0]=0;
    Ens[,,ngrid] = ens0[,3:217];rm(ens0)
  }# End square
  # Average square
  obs <- apply(Obs,c(1,2),mean);rm(Obs);obs <- obs[index,]
  Pens00 <- apply(Ens,c(1,2),mean);rm(Ens); Pens00 <- Pens00[index_f,]
  ## Climatology
  for (nyear in 1:dimIndex[1]){
    Pclim[nyear,] = colMeans(obs[-(nyear),1:215])
  }
  ## Verification raw
  Pens0 = Pens00[,s0[1]:s1[3]]
  gobs = obs[,s0[1]:s1[3]]
  gclim = Pclim[,s0[1]:s1[3]]
  ## Daily
  Pens <- array(Pens0,dim = dim(Pens0)[1]*dim(Pens0)[2])
  gobs0 <- matrix(gobs,nrow  = dimIndex*dim(Pens0)[2],ncol = 1)
  gclim0 <- matrix(gclim,nrow = dimIndex*dim(Pens0)[2],ncol = 1)
  ens_vr = sapply(0:(dimIndex*dim(Pens0)[2]-1), function(nind) Pens[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
  MAE_clim = mae(gobs0,gclim0)
  crpsr = round(mean(EnsCrps(ens_vr,gobs0)),digits=2)
  crpssr_clim[nsqr,1] = round(1-crpsr/MAE_clim,digits = 2)
  for (iup in 2:3){#Nup
    Pens1 <- upscale(Pens0,nop,iup)
    gobs1 <- upscale(gobs,nop,iup)
    gclim1 <- upscale(gclim,nop,iup)
    Pens <- array(Pens1,dim = dim(Pens1)[1]*dim(Pens1)[2])
    gobs1 <- matrix(gobs1,nrow  = dimIndex*dim(Pens1)[2],ncol = 1)
    gclim1 <- matrix(gclim1,nrow = dimIndex*dim(Pens1)[2],ncol = 1)
    ens_vr = sapply(0:(dimIndex*dim(Pens1)[2]-1), function(nind) Pens[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
    MAE_clim = mae(gobs1,gclim1)
    crpsr = round(mean(EnsCrps(ens_vr,gobs1)),digits=2)
    crpssr_clim[nsqr,iup] = round(1-crpsr/MAE_clim,digits = 2)
  }# End lead time raw
  ## Seasonal
  Pens1 <- upscale(Pens0,nop,4)
  gobs1 <- upscale(gobs,nop,4)
  gclim1 <- upscale(gclim,nop,4)
  ens_vr = sapply(0:(dimIndex-1), function(nind) Pens1[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
  MAE_clim = mae(gobs1,gclim1)
  crpsr = round(mean(EnsCrps(ens_vr,gobs1)),digits=2)
  crpssr_clim[nsqr,4] = round(1-crpsr/MAE_clim,digits = 2)
} #End grid

#Save verification 24(nsqr) x 4 upscaling scales
write(t(crpssr_clim),file=paste0(ObsDir,'/verification_daily/crpss_40_m',nmonth,'.txt'),ncolumns=4)
proc.time() - ptm


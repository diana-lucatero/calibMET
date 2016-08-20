##Computes verification 

ptm = proc.time()


ObsDir = 'E:/calibMET/precip'
month = seq(1,12,1);mlead = seq(1,7,1)
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
Nlead = 7
library(SpecsVerification);library(Metrics)
#library(matrixStats)
library(RNetCDF)

s0 <- c(1,31,61,91,121,151,181)
s1 <- c(30,60,90,120,150,180,215)
d  <- c(rep(30,6),35)

grid_filename = paste(ObsDir,'/obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
info_filename = paste(ObsDir,'/ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))

nmonth = 12
index     = which(ens_info[1:293,3] == month[nmonth])
#index_out = which(ens_info[index,4] != Nens[nmonth]) ## For Feb and Aug
#index     = index[-(index_out)]
ens_index = ens_info[index,c(31,32)]
no_ens    = Nens[nmonth];dimIndex  = length(index)
no_ens_clim = dimIndex[1]-1
index_f   = matrix(nrow=1, ncol=dimIndex*no_ens)
index_f   = sapply(1:dimIndex, function(ix) c(ens_index[ix,1]:ens_index[ix,2]))
Pclim = array(NaN, dim=c(dimIndex[1],215))

# Load thresholds
file <- paste0(ObsDir,'/drizzle/percentage_m',nmonth,'.nc')
nc <- open.nc(file)
thre <- var.get.nc(nc,'per')

# Verification arrays raw
crpssr_clim = array(NaN,dim=c(dim(grid)[1]))
len = array(NaN,dim=c(dim(grid)[1]))

for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  ## Observations
  obs_namefile = paste(ObsDir,'/observation/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_namefile),nrow=216,ncol=293));obs[obs<0]=0;
  obs = obs[index,2:216]
  ## Climatology
  for (nyear in 1:dimIndex[1]){
    Pclim[nyear,] = colMeans(obs[-(nyear),1:215])
  }
  ## Ensembles Raw
  ens_filename0 = paste(ObsDir,'/ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens0 = t(matrix(scan(ens_filename0),nrow=217,ncol=8322));ens0[ens0<0]=0;
  Pens00 = ens0[index_f,3:217];rm(ens0)
  ## Verification raw
  for (ilead in 1:1){
    thre0 <- thre[,,ilead,ngrid]; thre1 <- matrix(thre0,nrow = dim(thre0)[1]*dim(thre0)[2])
    thre1[which(is.na(thre1))]=0
    thre2 <-matrix(rep(thre1,30),ncol = 30,nrow=dim(thre0)[1]*dim(thre0)[2])
    rm(thre0);rm(thre1)
    Pens0 = Pens00[,s0[ilead]:s1[ilead]]
    ind <- which(Pens0 <= thre2)
    Pens0[ind]=0
    gobs = obs[,s0[ilead]:s1[ilead]]
    gclim = Pclim[,s0[ilead]:s1[ilead]]
    Pens <- array(Pens0,dim = dim(Pens0)[1]*dim(Pens0)[2])
    gobs0 <- matrix(gobs,nrow  = dimIndex*dim(Pens0)[2],ncol = 1)
    gclim0 <- array(gclim,dim = dimIndex*dim(Pens0)[2])
    ens_vr = sapply(0:(dimIndex*dim(Pens0)[2]-1), function(nind) Pens[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
    i0 <- which(gobs0==0) #i0 <- which(gobs0>1 & gobs0<=2 )
    MAE_clim = mae(gobs0[i0],gclim0[i0])
    crpsr = round(mean(EnsCrps(ens_vr[i0,],gobs0[i0])),digits=2)
    crpssr_clim[ngrid] = round(1-crpsr/MAE_clim,digits = 2)
    len[ngrid] <- length(i0)
    rm(Pens);rm(gobs0);rm(gclim0);rm(ens_vr)
  }
} #End grid

write(t(crpssr_clim),file=paste(ObsDir,'/drizzle/ver_C_',nmonth,'.txt',sep=''),ncolumns=1)
write(t(len),file=paste(ObsDir,'/drizzle/verL_C_',nmonth,'.txt',sep=''),ncolumns=1)

proc.time() - ptm


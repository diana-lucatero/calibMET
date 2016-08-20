## makes Q-Q plots obs-sim values

ObsDir = 'E:/calibMET/precip'
month = seq(1,12,1);mlead = seq(1,7,1)
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)

library(SpecsVerification);library(Metrics)
library(matrixStats);library(boot)

grid_filename = paste(ObsDir,'/obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
info_filename = paste(ObsDir,'/ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))

nmonth = 1
index     = which(ens_info[1:293,3] == month[nmonth])
#index_out = which(ens_info[index,4] != Nens[nmonth]) ## For Feb and Aug
#index     = index[-(index_out)]
ens_index = ens_info[index,c(31,32)]
no_ens    = Nens[nmonth];dimIndex  = length(index)
index_f   = matrix(nrow=1, ncol=dimIndex*no_ens)
index_f   = sapply(1:dimIndex, function(ix) c(ens_index[ix,1]:ens_index[ix,2]))
Pclim = array(NaN, dim=c(dimIndex[1],7))

ens0 = array(NaN,dim = c(375*30,662))
obs0 = array(NaN,dim = c(25*30,662))
for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  ## Observations
  obs_namefile = paste(ObsDir,'/observation/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_namefile),nrow=216,ncol=293));obs[obs<0]=0;
  obs = obs[index,2:216]; obs = obs[,1:30]
  ## Ensembles Raw
  ens_filename0 = paste(ObsDir,'/ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens = t(matrix(scan(ens_filename0),nrow=217,ncol=8322))
  ens = ens[index_f,3:217];ens = ens[,1:30]
  ens0[,ngrid] = matrix(ens,nrow = (375*30),ncol = 1)
  obs0[,ngrid] = matrix(obs,nrow = (25*30),ncol = 1)
  #;rm(ens0)
} # End grid
ens1 = matrix(ens0,nrow = (375*30*662),ncol = 1)
obs1 = matrix(obs0,nrow = (25*30*662),ncol = 1)
ens_vr = sapply(0:((dimIndex*30*662)-1), function(nind) ens1[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
dataOE <- data.frame(ens=ens_vr,obs=obs1)
## Conditional porbabilities h(y|fk) for precipitation

Nlead = 215 #Lead time in days
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)

ObsDir = 'G:/calibMET/precip/'
grid_filename = paste(ObsDir,'obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
info_filename = paste(ObsDir,'ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))

## Compute min and max , ens , obs for each grid
probs0 = c(0,0.25,0.5,0.75,1)
## Create function
myfunctionEns <- function(ObsDir,ngrid,probs0){
  ens_filename <- paste(ObsDir,'ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens <- t(matrix(scan(ens_filename),nrow=217,ncol=8322)); ens[ens<0]=0; ens = ens[,3:217]
  qEns <- quantile(ens,probs = probs0,na.rm= TRUE) 
  rm(ens)
  return(qEns)
}

myfunctionObs <- function(ObsDir,ngrid,probs0){
  obs_filename <- paste(ObsDir,'observation/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs <- t(matrix(scan(obs_filename),nrow=216,ncol=293));obs[obs<0]=0; obs = obs[,2:216]
  qObs <- quantile(obs,probs = probs0,na.rm= TRUE)
  rm(obs)
  return(qObs)
}


qEnsf = sapply(1:dim(grid)[1], function(x) myfunctionEns(ObsDir,x,probs0))
qObsf = sapply(1:dim(grid)[1], function(x) myfunctionObs(ObsDir,x,probs0))


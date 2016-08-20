##Computes verification for ensembles generated with correct_LS.R (08-06-16)

ObsDir = 'E:/calibMET/precip'
#meth = c('dc','qm')
meth = 'bma'
month = seq(1,12,1);mlead = seq(1,7,1)
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)

library(SpecsVerification);library(Metrics)
library(matrixStats);library(boot)

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
Pclim = array(NaN, dim=c(dimIndex[1],7))
# Verification arrays raw
biasr = array(NaN,dim=c(dim(grid)[1],7));corrr = array(NaN,dim=c(dim(grid)[1],7))
crpsr = array(NaN,dim=c(dim(grid)[1],7));crpssr_clim = array(NaN,dim=c(dim(grid)[1],7))
rh00 = array(NaN,dim=c(dim(grid)[1],7,(no_ens+1))); sharp00 = array(NaN,dim=c(dim(grid)[1],7,dimIndex))
# Verification arrays methods
biasm = array(NaN,dim=c(dim(grid)[1],7,length(meth)));corrm = array(NaN,dim=c(dim(grid)[1],7,length(meth)))
crpsm = array(NaN,dim=c(dim(grid)[1],7,length(meth)));crpssm_raw = array(NaN,dim=c(dim(grid)[1],7,length(meth)))
crpssm_clim = array(NaN,dim=c(dim(grid)[1],7,length(meth)))
rh = array(NaN,dim=c(dim(grid)[1],7,length(meth),(no_ens+1))); sharp = array(NaN,dim=c(dim(grid)[1],7,length(meth),dimIndex)) 
for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  ## Observations
  obs_namefile = paste(ObsDir,'/observation/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_namefile),nrow=216,ncol=293));obs[obs<0]=0;
  obs = obs[index,2:216]
  # Compute accumulated precipitation
  Pobs_m = cbind(rowSums(obs[,1:30]),rowSums(obs[,31:60]),rowSums(obs[,61:90]),rowSums(obs[,91:120]),rowSums(obs[,121:150]),rowSums(obs[,151:180]),rowSums(obs[,181:210]))
  rm(obs)
  ## Climatology
  for (nyear in 1:dimIndex[1]){
    Pclim[nyear,] = colMeans(Pobs_m[-(nyear),])
  }
  ## Ensembles Raw
  ens_filename0 = paste(ObsDir,'/ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens0 = t(matrix(scan(ens_filename0),nrow=217,ncol=8322))
  Pens00 = ens0[index_f,3:217];rm(ens0)
  # Compute accumulated precipitation raw
  Pens_m    = cbind(rowSums(Pens00[,1:30]),rowSums(Pens00[,31:60]),rowSums(Pens00[,61:90]),rowSums(Pens00[,91:120]),rowSums(Pens00[,121:150]),rowSums(Pens00[,151:180]),rowSums(Pens00[,181:210]))
  ## Verification raw
  for (ilead in 1:7){
    Pens0 = Pens_m[,ilead];gobs = Pobs_m[,ilead];gclim = Pclim[,ilead]
    MAE_clim = mae(gobs,gclim)
    biasr[ngrid,ilead] = round((mean(Pens0,na.rm = TRUE)/mean(gobs,na.rm = TRUE)-1)*100,digits=2) #(mean(gens0)/mean(gobs)-1)*100
    ens_vr = sapply(0:(dimIndex-1), function(nind) Pens0[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
    corrr[ngrid,ilead] = round(Corr(ens_vr,gobs)[1],digits=2)
    crpsr[ngrid,ilead] = round(mean(EnsCrps(ens_vr,gobs)),digits=2)
    rh00[ngrid,ilead,]   = Rankhist(ens_vr,gobs)
    crpssr_clim[ngrid,ilead] = round(1-crpsr[ngrid,ilead]/MAE_clim,digits = 2)
    sh = t(apply(ens_vr,1, quantile, probs = c(0.25,0.75),na.rm=TRUE))
    sharp00[ngrid,ilead,] = sh[,2]-sh[,1]
  }
  ## Ensembles methods
  for (nmet in 1:2){#length(meth)
    ens_namefile = paste(ObsDir,'/ensembles_c/ens_',meth[nmet],'_g',grid[ngrid,2],'_m',month[nmonth],'.txt',sep='')
    ens = t(matrix(scan(ens_namefile),nrow=215,ncol=dimIndex[1]*no_ens))
    # Compute accumulated precipitation methods
    Pens_met    = cbind(rowSums(ens[,1:30]),rowSums(ens[,31:60]),rowSums(ens[,61:90]),rowSums(ens[,91:120]),rowSums(ens[,121:150]),rowSums(ens[,151:180]),rowSums(ens[,181:210]))
    rm(ens)
    ## Verification methods
    for (ilead in 1:7){
      gens0 = Pens_met[,ilead];gobs = Pobs_m[,ilead];gclim = Pclim[,ilead]
      MAE_clim = mae(gobs,gclim)
      biasm[ngrid,ilead,nmet] = round((mean(gens0,na.rm = TRUE)/mean(gobs,na.rm = TRUE)-1)*100,digits=2)
      ens_v = sapply(0:(dimIndex-1), function(nind) gens0[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_v = t(ens_v)
      corrm[ngrid,ilead,nmet] = round(Corr(ens_v,gobs)[1],digits=2)
      crpsm[ngrid,ilead,nmet] = round(mean(EnsCrps(ens_v,gobs)),digits=2)
      rh[ngrid,ilead,nmet,] = Rankhist(ens_v,gobs)
      crpssm_raw[ngrid,ilead,nmet] = round(1-crpsm[ngrid,ilead,nmet]/crpsr[ngrid,ilead],digits = 2)
      crpssm_clim[ngrid,ilead,nmet] = round(1-crpsm[ngrid,ilead,nmet]/MAE_clim,digits = 2)
      sh0 = t(apply(ens_v,1, quantile, probs = c(0.25,0.75),na.rm=TRUE))
      sharp[ngrid,ilead,nmet,] = sh0[,2]-sh0[,1]
    } # end lead time methods
  } # End methods
} #End grid

# Save verification 662(ngrid) x 7(nlead)
# Raw
write(t(biasr),file=paste(ObsDir,'/verification/R_bias_m',nmonth,'.txt',sep=''),ncolumns=7)
write(t(corrr),file=paste(ObsDir,'/verification/R_corr_m',nmonth,'.txt',sep=''),ncolumns=7)
write(t(crpsr),file=paste(ObsDir,'/verification/R_crps_m',nmonth,'.txt',sep=''),ncolumns=7)
write(t(crpssr_clim),file=paste(ObsDir,'/verification/R_crpss_clim_m',nmonth,'.txt',sep=''),ncolumns=7)
for (il in 1:7) {
  write(t(rh00[,il,]),file=paste(ObsDir,'/verification/R_RH_m',nmonth,'_lt',il,'.txt',sep=''),ncolumns=(no_ens+1))
  write(t(sharp00[,il,]),file=paste(ObsDir,'/verification/R_sharp_m',nmonth,'_lt',il,'.txt',sep=''),ncolumns=dimIndex)
}
# Methods
for (imeth in 1:2){
  write(t(biasm[,,imeth]),file=paste(ObsDir,'/verification/',meth[imeth],'_bias_m',nmonth,'.txt',sep=''),ncolumns=7)
  write(t(corrm[,,imeth]),file=paste(ObsDir,'/verification/',meth[imeth],'_corr_m',nmonth,'.txt',sep=''),ncolumns=7)
  write(t(crpsm[,,imeth]),file=paste(ObsDir,'/verification/',meth[imeth],':crps_m',nmonth,'.txt',sep=''),ncolumns=7)
  write(t(crpssm_raw[,,imeth]),file=paste(ObsDir,'/verification/',meth[imeth],'_crpss_raw_m',nmonth,'.txt',sep=''),ncolumns=7)
  write(t(crpssm_clim[,,imeth]),file=paste(ObsDir,'/verification/',meth[imeth],'_crpss_clim_m',nmonth,'.txt',sep=''),ncolumns=7)
  for (il in 1:7) {
    write(t(rh[,il,imeth,]),file=paste(ObsDir,'/verification/',meth[imeth],'_RH_m',nmonth,'_lt',il,'.txt',sep=''),ncolumns=(no_ens+1))
    write(t(sharp[,il,imeth,]),file=paste(ObsDir,'/verification/',meth[imeth],'_sharp_m',nmonth,'_lt',il,'.txt',sep=''),ncolumns=dimIndex)
  }
}


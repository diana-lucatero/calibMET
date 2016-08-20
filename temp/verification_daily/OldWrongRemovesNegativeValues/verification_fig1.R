##Computes verification temperature daily

ptm = proc.time()

ObsDir = 'E:/calibMET/temp'
meth = c('dc','qm','bma')
month = seq(1,12,1);mlead = seq(1,7,1)
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
Nlead = 7
library(SpecsVerification);library(Metrics)
library(matrixStats)

s0 <- c(1,31,61,91,121,151,181)
s1 <- c(30,60,90,120,150,180,215)
d  <- c(rep(30,6),35)

grid_filename = paste(ObsDir,'/obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=5,ncol=724))
info_filename = paste(ObsDir,'/ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))

nmonth <- 12
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
biasr = array(NaN,dim=c(dim(grid)[1],Nlead))
crpssr_clim = array(NaN,dim=c(dim(grid)[1],Nlead))
rh00 = array(NaN,dim=c(dim(grid)[1],Nlead,(no_ens+1)))
sharp00 = array(NaN,dim=c(dim(grid)[1],Nlead))
# Verification arrays methods
biasm = array(NaN,dim=c(dim(grid)[1],Nlead,length(meth)))
crpssm_clim = array(NaN,dim=c(dim(grid)[1],Nlead,length(meth)))
rh = array(NaN,dim=c(dim(grid)[1],Nlead,length(meth),(no_ens+1)))
sharp = array(NaN,dim=c(dim(grid)[1],Nlead,length(meth))) 

for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  ## Observations
  obs_namefile = paste(ObsDir,'/observations/obs_station_',grid[ngrid,2],'.txt',sep='')
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
  for (ilead in 1:Nlead){#Nlead
    Pens0 = Pens00[,s0[ilead]:s1[ilead]]
    gobs = obs[,s0[ilead]:s1[ilead]]
    gclim = Pclim[,s0[ilead]:s1[ilead]]
    Pens <- array(Pens0,dim = dim(Pens0)[1]*dim(Pens0)[2])
    gobs0 <- matrix(gobs,nrow  = dimIndex*dim(Pens0)[2],ncol=1)
    gclim0 <- matrix(gclim,nrow = dimIndex*dim(Pens0)[2],ncol=1)
    ens_vr = sapply(0:(dimIndex*dim(Pens0)[2]-1), function(nind) Pens[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
    # Bias
    biasr[ngrid,ilead] = round((mean(Pens0,na.rm = TRUE)/mean(gobs,na.rm = TRUE)-1)*100,digits=2)
    # CRPSS
    MAE_clim = mae(gobs0,gclim0)
    crpsr = round(mean(EnsCrps(ens_vr,gobs0)),digits=2)
    crpssr_clim[ngrid,ilead] = round(1-crpsr/MAE_clim,digits = 2)
    # RH
    rh00[ngrid,ilead,]   = Rankhist(ens_vr,gobs0)
    # Sharpness
    sh = t(apply(ens_vr,1, quantile, probs = c(0.25,0.75),na.rm=TRUE))
    sharp00[ngrid,ilead] = mean(sh[,2]-sh[,1])
    rm(Pens);rm(gobs0);rm(gclim0);rm(ens_vr)
  }
  ## Ensembles methods
  for (nmet in 1:length(meth)){#length(meth)
    ens_namefile = paste(ObsDir,'/ensembles_c/ens_',meth[nmet],'_g',grid[ngrid,2],'_m',month[nmonth],'.txt',sep='')
    ens = t(matrix(scan(ens_namefile),nrow=215,ncol=dimIndex[1]*no_ens));ens[ens<0]=0;
    ## Verification methods
    for (ilead in 1:Nlead){#Nlead
      gens0 = ens[,s0[ilead]:s1[ilead]]
      gobs = obs[,s0[ilead]:s1[ilead]]
      gclim = Pclim[,s0[ilead]:s1[ilead]]
      Pens <- array(gens0,dim = dim(gens0)[1]*dim(gens0)[2])
      gobs0 <- matrix(gobs,nrow  = dimIndex*dim(gens0)[2],ncol = 1)
      gclim0 <- matrix(gclim,nrow  = dimIndex*dim(gens0)[2],ncol = 1)
      ens_v = sapply(0:(dimIndex*dim(gens0)[2]-1), function(nind) Pens[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_v = t(ens_v)
      # Bias
      biasm[ngrid,ilead,nmet] = round((mean(gens0,na.rm = TRUE)/mean(gobs,na.rm = TRUE)-1)*100,digits=2)
      # Rank histogram
      rh[ngrid,ilead,nmet,] = Rankhist(ens_v,gobs0)
      # CRPSS
      crpsm = round(mean(EnsCrps(ens_v,gobs0)),digits=2)
      MAE_clim = mae(gobs0,gclim0)
      crpssm_clim[ngrid,ilead,nmet] = round(1-crpsm/MAE_clim,digits = 2)
      # Sharpness
      sh0 = t(apply(ens_v,1, quantile, probs = c(0.25,0.75),na.rm=TRUE))
      sharp[ngrid,ilead,nmet] = mean(sh0[,2]-sh0[,1])
    } # end lead time methods
  } # End methods
} #End grid

#Save verification 662(ngrid) x 7(nlead)
#Raw
write(t(biasr),file=paste(ObsDir,'/verification_daily/R_bias_m',nmonth,'.txt',sep=''),ncolumns=Nlead)
write(t(crpssr_clim),file=paste(ObsDir,'/verification_daily/R_crpss_clim_m',nmonth,'.txt',sep=''),ncolumns=Nlead)
write(t(sharp00),file=paste(ObsDir,'/verification_daily/R_sharp_m',nmonth,'.txt',sep=''),ncolumns=Nlead)
for (il in 1:Nlead) {
  write(t(rh00[,il,]),file=paste(ObsDir,'/verification_daily/R_RH_m',nmonth,'_lt',il,'.txt',sep=''),ncolumns=(no_ens+1))
}# End lead

# Methods
for (imeth in 1:length(meth)){
  write(t(biasm[,,imeth]),file=paste(ObsDir,'/verification_daily/',meth[imeth],'_bias_m',nmonth,'.txt',sep=''),ncolumns=Nlead)
  write(t(crpssm_clim[,,imeth]),file=paste(ObsDir,'/verification_daily/',meth[imeth],'_crpss_clim_m',nmonth,'.txt',sep=''),ncolumns=Nlead)
  write(t(sharp[,,imeth]),file=paste(ObsDir,'/verification_daily/',meth[imeth],'_sharp_m',nmonth,'.txt',sep=''),ncolumns=Nlead)
  for (il in 1:Nlead) {
    write(t(rh[,il,imeth,]),file=paste(ObsDir,'/verification_daily/',meth[imeth],'_RH_m',nmonth,'_lt',il,'.txt',sep=''),ncolumns=(no_ens+1))
  }# end lead
}# End meth



proc.time() - ptm


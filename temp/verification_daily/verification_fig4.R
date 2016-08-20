##Computes verification temperature daily

ptm = proc.time()

ObsDir = 'G:/calibMET/temp'
meth = c('dc','qm','bma')
month = seq(1,12,1);mlead = seq(1,7,1)
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
Nlead = 7
library(SpecsVerification);library(Metrics)
library(matrixStats)

s0 <- c(1,31,61,91,121,151,181)
s1 <- c(30,60,90,120,150,180,215)
d  <- c(rep(30,6),35)

#Thresholds
nthre <- 5
thre <- matrix(c(-5,5,5,15,15,25),nrow = 3,ncol = 2,byrow = T)

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
crpssr_clim = array(NaN,dim=c(dim(grid)[1],Nlead,nthre))
len = array(NaN,dim=c(dim(grid)[1],Nlead,nthre))
# Verification arrays methods
crpssm_clim = array(NaN,dim=c(dim(grid)[1],Nlead,length(meth),nthre))
lenm = array(NaN,dim=c(dim(grid)[1],Nlead,length(meth),nthre))

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
  for (ilead in 1:Nlead){#Nlead
    Pens0 = Pens00[,s0[ilead]:s1[ilead]]
    gobs = obs[,s0[ilead]:s1[ilead]]
    gclim = Pclim[,s0[ilead]:s1[ilead]]
    Pens <- array(Pens0,dim = dim(Pens0)[1]*dim(Pens0)[2])
    gobs0 <- matrix(gobs,nrow  = dimIndex*dim(Pens0)[2],ncol=1)
    gclim0 <- matrix(gclim,nrow = dimIndex*dim(Pens0)[2],ncol=1)
    ens_vr = sapply(0:(dimIndex*dim(Pens0)[2]-1), function(nind) Pens[((nind*no_ens)+1):((nind*no_ens)+no_ens)]); ens_vr = t(ens_vr)
    
    ## Threshold 1
    i0 <- which(gobs0< -5)
    if (length(i0) == 0) {crpssr_clim[ngrid,ilead,1] = NaN}
    else {
      MAE_clim = mae(gobs0[i0],gclim0[i0])
      crpsr = round(mean(EnsCrps(ens_vr[i0,],gobs0[i0])),digits=2)
      crpssr_clim[ngrid,ilead,1] = round(1-crpsr/MAE_clim,digits = 2)
      len[ngrid,ilead,1] <- length(i0)
    }#end if
    rm(i0)
    
    ## Threshold 2:4
    for (it in 1:3){
      i0 <- which(gobs0>thre[it,1] & gobs0<=thre[it,2] )
      if (length(i0) == 0) {crpssr_clim[ngrid,ilead,(it+1)] = NaN}
      else {
        MAE_clim = mae(gobs0[i0],gclim0[i0])
        crpsr = round(mean(EnsCrps(ens_vr[i0,],gobs0[i0])),digits=2)
        crpssr_clim[ngrid,ilead,(it+1)] = round(1-crpsr/MAE_clim,digits = 2)
        len[ngrid,ilead,(it+1)] <- length(i0)
      }#end if
      rm(i0)
    }
    
    ## Threshold 5
    i0 <- which(gobs0>25) #i0 <- which(gobs0>1 & gobs0<=2 )
    if (length(i0) == 0) {crpssr_clim[ngrid,ilead,5] = NaN}
    else {
      MAE_clim = mae(gobs0[i0],gclim0[i0])
      crpsr = round(mean(EnsCrps(ens_vr[i0,],gobs0[i0])),digits=2)
      crpssr_clim[ngrid,ilead,5] = round(1-crpsr/MAE_clim,digits = 2)
      len[ngrid,ilead,5] <- length(i0)
    }#end if
    rm(i0)
  }# End lead time raw
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
      
      ## Threshold 1
      i0 <- which(gobs0< -5)
      if (length(i0) == 0) {crpssm_clim[ngrid,ilead,nmet,1] = NaN}
      else {
        MAE_clim = mae(gobs0[i0],gclim0[i0])
        crpsm = round(mean(EnsCrps(ens_v[i0,],gobs0[i0])),digits=2)
        crpssm_clim[ngrid,ilead,nmet,1] = round(1-crpsm/MAE_clim,digits = 2)
        lenm[ngrid,ilead,nmet,1] <- length(i0)
      }# End if
      rm(i0)
      
      ## Threshold 2:4
      for (it in 1:3){
        i0 <- which(gobs0>thre[it,1] & gobs0<=thre[it,2] )
        if (length(i0) == 0) {crpssm_clim[ngrid,ilead,nmet,(it+1)] = NaN}
        else {
          MAE_clim = mae(gobs0[i0],gclim0[i0])
          crpsm = round(mean(EnsCrps(ens_v[i0,],gobs0[i0])),digits=2)
          crpssm_clim[ngrid,ilead,nmet,(it+1)] = round(1-crpsm/MAE_clim,digits = 2)
          lenm[ngrid,ilead,nmet,(it+1)] <- length(i0)
        }# End if
        rm(i0)
      }
      
      ## Threshold 5
      i0 <- which(gobs0>25)
      if (length(i0) == 0) {crpssm_clim[ngrid,ilead,nmet,5] = NaN}
      else {
        MAE_clim = mae(gobs0[i0],gclim0[i0])
        crpsm = round(mean(EnsCrps(ens_v[i0,],gobs0[i0])),digits=2)
        crpssm_clim[ngrid,ilead,nmet,5] = round(1-crpsm/MAE_clim,digits = 2)
        lenm[ngrid,ilead,nmet,5] <- length(i0)
      }# End if
      rm(i0)

    } # end lead time methods
  } # End methods
} #End grid

#Save verification 662(ngrid) x 7(nlead)
#Raw
for (it in 1:5){
  write(t(crpssr_clim[,,it]),file=paste(ObsDir,'/verification_daily/R_crpss_clim_m',nmonth,'_thre',it,'.txt',sep=''),ncolumns=Nlead)
  write(t(len[,,it]),file=paste(ObsDir,'/verification_daily/R_len_clim_m',nmonth,'_thre',it,'.txt',sep=''),ncolumns=Nlead)
}

# Methods
for (imeth in 1:length(meth)){
  for (it in 1:5){
    write(t(crpssm_clim[,,imeth,it]),file=paste(ObsDir,'/verification_daily/',meth[imeth],'_crpss_clim_m',nmonth,'_thre',it,'.txt',sep=''),ncolumns=Nlead)
    write(t(lenm[,,imeth,it]),file=paste(ObsDir,'/verification_daily/',meth[imeth],'_len_clim_m',nmonth,'_thre',it,'.txt',sep=''),ncolumns=Nlead)
  }
}# End meth



proc.time() - ptm


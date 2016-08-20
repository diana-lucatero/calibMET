## Plots densities per grid point

ObsDir = 'E:/calibMET/precip/'
grid_filename = paste(ObsDir,'obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
info_filename = paste(ObsDir,'ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))
prob0 = c(0,0.5,1)
## Compute corrected monthly agregated
for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  # Observed
  obs_filename = paste(ObsDir,'observation/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_filename),nrow=216,ncol=293));obs[obs<0]=0;obs = obs[,2:216]
  ## Compute stats obs
  qobs = quantile(obs,probs = prob0,na.rm = TRUE); qobs = round(qobs, digits = 2)
  # Ensembles
  ens_filename = paste(ObsDir,'ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322)); ens[ens<0]=0;ens = ens[,3:217]
  ## Compute stats ens
  qens = quantile(ens,probs = prob0,na.rm = TRUE); qens = round(qens, digits = 2)
  ## Plots ens
  name_plt = paste(ObsDir,'plots/ens_',ngrid,'.png',sep='')
  png(name_plt)
  h = hist(ens,breaks =1000)
  y0 = max(h$counts)/2; x0 = qens[3]/2
  plot(h)
  text0 = paste('min = ',qens[1],'\nmean = ',qens[2],'\nmax = ',qens[3],sep='')
  text(x0,y0,text0)
  ## Plots obs
  name_plt = paste(ObsDir,'plots/obs_',ngrid,'.png',sep='')
  png(name_plt)
  h = hist(obs,breaks =1000)
  y0 = max(h$counts)/2; x0 = qobs[3]/2
  plot(h)
  text0 = paste('min = ',qobs[1],'\nmean = ',qobs[2],'\nmax = ',qobs[3],sep='')
  text(x0,y0,text0)
  graphics.off()
}

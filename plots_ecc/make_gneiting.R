## Reproduces Gneiting (2013) of year 2003/02 with lead time = 30 day
## Precipitation and temperature, nodes 646 and 647

ObsDir = 'E:/calibMET/'
library(psych)
meth <- c('dc','qm','bma')
nmeth <- 3

year = 2003; nmonth = 2; index = 157

grid = c(646,647)
info_filename = paste(ObsDir,'precip/ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))
ens_index = ens_info[index,c(31,32)];
ens_index = c(ens_index[1]:ens_index[2])
no_ens    = length(ens_index)

ensT = array(NaN,dim=c(2,no_ens));ensP = array(NaN,dim=c(2,no_ens))
ensTc = array(NaN,dim=c(2,no_ens));ensPc = array(NaN,dim=c(2,no_ens))
for (ngrid in 1:length(grid)){ # Loop for space :dim(grid)[1]
  print(ngrid)
  ## Raw
  # Temperature
  ens_filename = paste(ObsDir,'temp/ensembles_idw/ens_station_',grid[ngrid],'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322)) #ens[ens<0]=0;
  ens = ens[,3:217];ensT[ngrid,] = ens[ens_index,30]; rm(ens)
  # Precipitation
  ens_filename = paste(ObsDir,'precip/ensembles_idw/ens_station_',grid[ngrid],'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322)); ens[ens<0]=0;
  ens = ens[,3:217];ensP[ngrid,] = ens[ens_index,30]; rm(ens)
  ## Corrected BMA
  # Temperature
  ens_filename = paste(ObsDir,'temp/ensembles_c/ens_',meth[nmeth],'_g',grid[ngrid],'_m2.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=215,ncol=(24*51))); #ens[ens<0]=0;
  ensTc[ngrid,] = ens[613:663,30]; rm(ens)
  # Precipitation
  ens_filename = paste(ObsDir,'precip/ensembles_c/ens_',meth[nmeth],'_g',grid[ngrid],'_m2.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=215,ncol=(24*51))); ens[ens<0]=0;
  ensPc[ngrid,] = ens[613:663,30]; rm(ens)
}

# # Plot
# #m <- matrix(seq(1:3),nrow = 1,ncol = 3,byrow = TRUE) # No Meth + 1 (for legend)
# name_plt = paste(ObsDir,'plots_ecc/plot_raw.png',sep='')
# png(name_plt)
# #layout(mat = m)
# 
# 
# data = t(rbind(ensT,ensP))
# colnames(data) <- c('T646','T647','P646','P647') 
# pairs.panels(data,main = 'raw 2003 02 lt +30 days')
# graphics.off()
# 
# name_plt = paste(ObsDir,'plots_ecc/plot_',meth[nmeth],'.png',sep='')
# png(name_plt)
# data = t(rbind(ensTc,ensPc))
# colnames(data) <- c('T646','T647','P646','P647') 
# main0 = paste0(meth[nmeth],' 2003 02 lt +30 days')
# pairs.panels(data,main = main0)
# graphics.off()
# 
# ## Ranks
# 
# rkP1 = rank(ensP[1,],ties.method = 'first')
# rkP2 = rank(ensP[2,],ties.method = 'first')
# rkT1 = rank(ensT[1,],ties.method = 'first')
# rkT2 = rank(ensT[2,],ties.method = 'first')
# 
# ensTc_r = array(NaN,dim=c(2,no_ens));ensPc_r = array(NaN,dim=c(2,no_ens))
# for (iex in 1:no_ens){
#   ind0 = which(rank(ensPc[1,],ties.method='first') == rkP1[iex])
#   ensPc_r[1,iex] = ensPc[1,ind0]
#   ind0 = which(rank(ensPc[2,],ties.method='first') == rkP2[iex])
#   ensPc_r[2,iex] = ensPc[2,ind0]
# }
# 
# for (iex in 1:no_ens){
#   ind0 = which(rank(ensTc[1,],ties.method='first') == rkT1[iex])
#   ensTc_r[1,iex] = ensTc[1,ind0]
#   ind0 = which(rank(ensTc[2,],ties.method='first') == rkT2[iex])
#   ensTc_r[2,iex] = ensTc[2,ind0]
# }
# 
# 
# name_plt = paste(ObsDir,'plots_ecc/plot_',meth[nmeth],'_ecc.png',sep='')
# png(name_plt)
# data = t(rbind(ensTc_r,ensPc_r))
# colnames(data) <- c('T646','T647','P646','P647') 
# main0 = paste0(meth[nmeth],' ECC 2003 02 lt +30 days')
# pairs.panels(data,main = main0)
# graphics.off()
# 
# 

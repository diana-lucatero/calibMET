## Computes and saves ranks
## Check Feb and Aug -- they might differ
Nlead = 215 #Lead time in days
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)

Dir = 'C:/Users/DianaL/Desktop/NCAR_visit/calibrateMET/'

# grid_filename = paste(Dir,'obs_ens.txt',sep='')
# grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
grid = seq(1,72,1)
info_filename = paste(Dir,'ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))
no_fore = 310

for (ngrid in 1:length(grid)){ # Loop for space :length(grid)
  print(ngrid)
  ens_filename = paste(Dir,'data/precip/ensembles_nc/ens_station_',grid[ngrid],'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322)); ens[ens<0]=0; nams = ens[,1:2];ens = ens[,3:217]
  Pens_r = array(NaN, dim = c(dim(ens)[1],dim(ens)[2]))
  for (ilead in 1:Nlead){#Nlead
    Pens = ens[,ilead]
    ens_vr = sapply(0:(no_fore-1), function(nind) rank(Pens[ens_info[(nind+1),31]:ens_info[(nind+1),32]],ties.method='first'))
    ens_vr = unlist(ens_vr)
    Pens_r[,ilead] = ens_vr
  }
  rkEns = cbind(nams,Pens_r)
  ## Save ranks
  name_file = paste(Dir,'data/precip/ensembles_r/ens_station_',grid[ngrid],'.txt',sep='')
  write.table(rkEns, file=name_file, row.names=FALSE, col.names=FALSE)
}

## Check which files are missing
## 12 months, 662 grid points

Dir = 'E:/calibMET/precip/ensembles_c_bma/'

grid_filename = paste(Dir,'obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))

exist = array(NaN, dim = c(12,662))

for (nmonth in 1:12){#12
  for (ngrid in 1:662){#662
    name_file = paste0(Dir,'ens_bma_g',grid[ngrid,2],'_m',nmonth,'.txt')
    if (find.file(name_file) == ''){ exist[nmonth,ngrid] = 0 } 
    else { exist[nmonth,ngrid] = 1}
  }
}

write(t(exist),file=paste(Dir,'exist.txt',sep=''),ncolumns=662)
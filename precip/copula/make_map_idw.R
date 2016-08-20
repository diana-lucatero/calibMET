## Makes map of year 2003/02 with lead time = 30 days and ensemble 1
## Converts UTM (+zone = 32) to longlat

library(rgdal)

## Coordinates idw
Dir = 'C:/Users/DianaL/Desktop/NCAR_visit/calibrateMET/plots/'
name_f = paste(Dir,'maps/obs_ens_coord_obs.txt',sep="")
coor   = t(matrix(scan(name_f),nrow=3,ncol=662))


year = 2003; nmonth = 2; index = 157
ObsDir = 'E:/calibMET/precip/'
grid_filename = paste(ObsDir,'obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
info_filename = paste(ObsDir,'ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))
ens_index = ens_info[index,c(31,32)];
ens_index = c(ens_index[1]:ens_index[2])
no_ens    = length(ens_index)
## read closer nodes
node_fn = paste(ObsDir,'ranks/nodes_closer.txt',sep='')
nodes = t(matrix(scan(node_fn),nrow=4,ncol=662))
nodes = nodes[,4]


ens01 = array(NaN,dim = 662); obs0 = array(NaN,dim = 662)
ens02 = array(NaN,dim = 662); ens03 = array(NaN,dim = 662);
ens04 = array(NaN,dim = 662); ens05 = array(NaN,dim = 662)
for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  obs_filename = paste(ObsDir,'observation/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_filename),nrow=216,ncol=293));obs[obs<0]=0;
  ens_filename = paste(ObsDir,'ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322)); ens[ens<0]=0;
  # Extract ensembles_idw
  Pens      = ens[ens_index,3:217]; ens01[ngrid] = Pens[1,30]; rm(Pens);rm(ens)
  # Extract observations
  Pobs = obs[index,2:216]; obs0[ngrid] = Pobs[30];rm(Pobs);rm(obs)
  # Extract ranks
  rank_filename = paste(ObsDir,'ensembles_r_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  rank = t(matrix(scan(rank_filename),nrow=217,ncol=8322))
  Prank      = rank[ens_index,3:217]; Prank0 = Prank[1,30] # First ensemble
  # Extract ensembles_c/dc
  ens_filename = paste(ObsDir,'ensembles_c/ens_dc_g',grid[ngrid,2],'_m',nmonth,'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=215,ncol=(24*51))); ens[ens<0]=0;
  ens02[ngrid] = ens[613,30]
  ## ECC_DC
  Pens1 = ens[613:663,30]
  ind0 = which(rank(Pens1,ties.method='first') == Prank0)
  ens03[ngrid] = Pens1[ind0]
  rm(ens);rm(Pens1);rm(ind0)
  # Extract ensembles_c/qm
  ens_filename = paste(ObsDir,'ensembles_c/ens_qm_g',grid[ngrid,2],'_m',nmonth,'.txt',sep='')
  ens = t(matrix(scan(ens_filename),nrow=215,ncol=(24*51))); ens[ens<0]=0;
  ens04[ngrid] = ens[613,30]
  ## ECC_QM
  Pens2 = ens[613:663,30]
  ind0 = which(rank(Pens2,ties.method='first') == Prank0)
  ens05[ngrid] = Pens2[ind0]
  rm(ens);rm(Pens2);rm(ind0)
}

# Plot
m <- matrix(seq(1:7),nrow = 1,ncol = 7,byrow = TRUE) # No Meth + 1 (for legend)
name_plt = paste(ObsDir,'copula/testmap_idw.png',sep='')
png(name_plt,width=750,height=200)
layout(mat = m, widths=c(2,2,2,2,2,2,0.5))
ColorRamp <- terrain.colors(1000, alpha = 1)
values = cbind(obs0,ens01,ens02,ens03,ens04,ens05)
min = min(values,na.rm = TRUE);max=max(values,na.rm = TRUE)
ColorLevels <- seq(min, max, length=length(ColorRamp))

par(mar = c(0,0,2,0))
grd1 <- data.frame(z = obs0, xc = coor[,2], yc = coor[,3])
coordinates(grd1) <- ~ xc + yc
gridded(grd1) <- TRUE
grd1 <- as(grd1, "SpatialGridDataFrame")
image(grd1, col = ColorRamp, main = 'obs')

par(mar = c(0,0,2,0))
grd2 <- data.frame(z = ens01, xc = coor[,2], yc = coor[,3])
coordinates(grd2) <- ~ xc + yc
gridded(grd2) <- TRUE
grd2 <- as(grd2, "SpatialGridDataFrame")
image(grd2, col = ColorRamp, main = 'idw')

par(mar = c(0,0,2,0))
grd3 <- data.frame(z = ens02, xc = coor[,2], yc = coor[,3])
coordinates(grd3) <- ~ xc + yc
gridded(grd3) <- TRUE
grd3 <- as(grd3, "SpatialGridDataFrame")
image(grd3, col =ColorRamp, main = 'dc')

par(mar = c(0,0,2,0))
grd4 <- data.frame(z = ens03, xc = coor[,2], yc = coor[,3])
coordinates(grd4) <- ~ xc + yc
gridded(grd4) <- TRUE
grd4 <- as(grd4, "SpatialGridDataFrame")
image(grd4, col = ColorRamp, main = 'dc_ecc')

par(mar = c(0,0,2,0))
grd5 <- data.frame(z = ens04, xc = coor[,2], yc = coor[,3])
coordinates(grd5) <- ~ xc + yc
gridded(grd5) <- TRUE
grd5 <- as(grd5, "SpatialGridDataFrame")
image(grd5, col = ColorRamp, main = 'qm')

par(mar = c(0,0,2,0))
grd6 <- data.frame(z = ens05, xc = coor[,2], yc = coor[,3])
coordinates(grd6) <- ~ xc + yc
gridded(grd6) <- TRUE
grd6 <- as(grd6, "SpatialGridDataFrame")
image(grd6, col = ColorRamp, main = 'qm_ecc')

#Plot legend
par(mar = c(3,2,2.5,1))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

graphics.off()


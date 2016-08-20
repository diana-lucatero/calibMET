## Converts UTM (+zone = 32) to longlat

library(rgdal)

Dir = 'C:/Users/DianaL/Desktop/NCAR_visit/calibrateMET/plots/'
  
name_f = paste(Dir,'maps/obs_ens_coord_obs.txt',sep="")
coor   = t(matrix(scan(name_f),nrow=3,ncol=662))

utmcoor <- SpatialPoints(cbind(coor[,2],coor[,3]),proj4string = CRS('+proj=utm +zone=32'))
longlatcoor <- spTransform(utmcoor,CRS('+proj=longlat'))
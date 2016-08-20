## Maps of better methods
## http://rspatial.r-forge.r-project.org/gallery/#fig09.R
## http://www.r-bloggers.com/grid2polygons-2/
# Load libraries
#source("http://www.phaget4.org/R/myImagePlot.R")
quality = c('bias','corr','crpss_clim')
meth = c('R','dc','qm','bma')
# Load coord files
ObsDir = 'E:/calibMET/precip/'
Nmonth = 12;Nlead = 30
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
  
iq = 3 # quality
nmonth = 1 # month

values = array(NaN,dim=c(length(meth),662,Nlead))
for (imeth in 1:length(meth)){#length(meth)
      name_file = paste(ObsDir,'verification_weekly/',meth[imeth],'_',quality[iq],'_m',nmonth,'.txt',sep="")
      val = t(matrix(scan(name_file),nrow=Nlead,ncol=662))
      values[imeth,,] = val
} # End method

    #Plotting

## Coordinates idw
library(gdal)
Dir = 'C:/Users/DianaL/Desktop/NCAR_visit/calibrateMET/plots/'
name_f = paste(Dir,'maps/obs_ens_coord_obs.txt',sep="")
coor   = t(matrix(scan(name_f),nrow=3,ncol=662))

# Colors
ColorRamp <- rainbow(1000, alpha = 1)
#min = min(values,na.rm = TRUE);max=max(values,na.rm = TRUE)
min = -100; max = 100
min = -1; max = 1
ColorLevels <- seq(min, max, length=length(ColorRamp))

 
for (lt in 1:Nlead) { # Start lead time
  m <- matrix(seq(1:5),nrow = 1,ncol = 5,byrow = TRUE) # No Meth + 1 (for legend)
  name_plt = paste(ObsDir,'verification_weekly/plots/',quality[iq],'_m',nmonth,'_lt',lt,'.png',sep='')
  png(name_plt,width=750,height=200)
  #par(cex.lab=3)
  layout(mat = m, widths=c(2,2,2,2,0.5), heights=array(1,dim=c(5)))

  # Plot methods
  for (imet in 1:length(meth)){
    grd1 <- data.frame(z = values[imet,,lt], xc = coor[,2], yc = coor[,3])
    coordinates(grd1) <- ~ xc + yc
    gridded(grd1) <- TRUE
    grd1 <- as(grd1, "SpatialGridDataFrame")
    par(mar = c(0,0,2,0))
    main0 = meth[imet] 
    image(grd1, col = ColorRamp, main = main0,zlim=c(min,max))
  }
    #Plot legend
  par(mar = c(3,2,2.5,1))
  image(1, ColorLevels,
    matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
    col=ColorRamp,
    xlab="",ylab="",
    xaxt="n")
  
  graphics.off()
} # End lead time


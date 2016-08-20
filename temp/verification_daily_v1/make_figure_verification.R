## Maps of better methods
## http://rspatial.r-forge.r-project.org/gallery/#fig09.R
## http://www.r-bloggers.com/grid2polygons-2/
# Load libraries
#source("http://www.phaget4.org/R/myImagePlot.R")
quality = c('bias','corr','crpss_clim')
meth = c('R','dc','qm','bma')
# Load coord files
ObsDir = 'F:/ECMWF_Seasonal_data/BC_lt_T/calibMET/'
Nmonth = 12;Nlead = 7
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

for (iq in 1:length(quality)){
  values = array(NaN,dim=c(Nmonth,Nlead,length(meth)))
  for (imeth in 1:length(meth)){#length(meth)
    for (nmonth in 1:Nmonth){
      name_file = paste(ObsDir,'verification_daily/',meth[imeth],'_',quality[iq],'_m',nmonth,'.txt',sep="")
      val = t(matrix(scan(name_file),nrow=7,ncol=724))
      values[nmonth,,imeth] = colMeans(val,na.rm  = TRUE)
    } # end loop month
  }

  #Plotting

  m <- matrix(seq(1:5),nrow = 1,ncol = 5,byrow = TRUE) # No Meth + 1 (for legend)
  name_plt = paste(ObsDir,'verification_daily/',quality[iq],'.png',sep='')
  png(name_plt,width=750,height=200)
  #par(cex.lab=3)
  layout(mat = m, widths=c(2,2,2,2,0.5), heights=array(1,dim=c(5)))
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                  seq(0,1,length=256),  # Green
                  seq(1,0,length=256))  # Blue

  reverse = seq(12,1,-1)
  yLabels <- Mname#Mname[reverse]
  xLabels <- seq(7,1,-1)
  min = min(values,na.rm = TRUE);max=max(values,na.rm = TRUE)
  ColorLevels <- seq(min, max, length=length(ColorRamp))

  # Plot methods
  for (imet in 1:4){
    dat <- values[,,imet]
    #reverse <- nrow(dat) : 1
    #yLabels <- yLabels[reverse]
    dat <- dat[reverse,]
    par(mar = c(3,2,2.5,1))
    image(t(dat),col=ColorRamp , xlab="",
      ylab="", axes=FALSE, zlim=c(min,max))
    #image(t(dat),col=ColorRamp, zlim=c(min,max))
    axis(1,at=1:length(xLabels), labels=xLabels, cex.axis=1)
    axis(2,at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
     cex.axis=1)
  }
  #Plot legend
  par(mar = c(3,2,2.5,1))
  image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")
} # End qualities
graphics.off()
## Makes heat map with crpss_clim and sharpness

## Makes Heat Maps of Bias for PCP, TEMP and ET0

quality <- c('crpss_clim','sharp')
meth = c('R','dc','qm','bma')
nmeth <- c('RAW','LS','QM','BMA')
library(sp);library(RColorBrewer)
# Load coord files
ObsDir = 'G:/calibMET/'
Nmonth = 12;Nlead = 7
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

reverse = seq(12,1,-1)
yLabels <- Mname[reverse]
xLabels <- seq(1,7,1)

valuesP = array(NaN,dim=c(Nmonth,Nlead,length(meth),length(quality)))
valuesT = array(NaN,dim=c(Nmonth,Nlead,length(meth),length(quality)))
valuesE = array(NaN,dim=c(Nmonth,Nlead,length(meth),length(quality)))
for (imeth in 1:length(meth)){#length(meth)
  for (nmonth in 1:Nmonth){
    ## Precipitation 
    name_file = paste(ObsDir,'precip/verification_daily/',meth[imeth],'_crpss_clim_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=662))
    valuesP[nmonth,,imeth,1] = colMeans(val,na.rm  = TRUE)
    name_file = paste(ObsDir,'precip/verification_daily/',meth[imeth],'_sharp_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=662))
    valuesP[nmonth,,imeth,2] = colMeans(val,na.rm  = TRUE)
    ## Temperature
    name_file = paste(ObsDir,'temp/verification_daily/',meth[imeth],'_crpss_clim_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=724))
    valuesT[nmonth,,imeth,1] = colMeans(val,na.rm  = TRUE)
    name_file = paste(ObsDir,'temp/verification_daily/',meth[imeth],'_sharp_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=724))
    valuesT[nmonth,,imeth,2] = colMeans(val,na.rm  = TRUE)
    ## Ref. Evaporanspiration
    name_file = paste(ObsDir,'ref_evap/verification_daily/',meth[imeth],'_crpss_clim_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=724))
    valuesE[nmonth,,imeth,1] = colMeans(val,na.rm  = TRUE)
    name_file = paste(ObsDir,'ref_evap/verification_daily/',meth[imeth],'_sharp_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=724))
    valuesE[nmonth,,imeth,2] = colMeans(val,na.rm  = TRUE)
  } # end loop month
}

# Plotting
x0 <- c(1:7)
y0 <- c(1:12)

c0 <- c(rep(1,5),seq(2,6),rep(7,5),seq(8,12),rep(13,5),seq(14,18))    
m <- matrix(c0,nrow = 6,ncol = 5,byrow = TRUE) # No Meth + 1 (for legend)

#c0 <- c(rep(1,4),seq(2,5),rep(6,4),seq(7,10))    
#m <- matrix(c0,nrow = 4,ncol = 4,byrow = TRUE) # No Meth + 1 (for legend)

name_plt = paste(ObsDir,'PQ_figures/fig2/fig2.png',sep='')
png(name_plt,width=900,height=750)
#par(cex.lab=3)
layout(mat = m, widths=c(2,2,2,2,2.5), heights=array(c(0.3,1,0.3,1,0.3,1)))
#layout(mat = m, widths=c(2,2,2,2),heights=array(c(0.3,1,0.3,1)))
#colors <- c('burlywood','darkgoldenrod1','firebrick1','darkmagenta','darkgreen','brown4')
pal = colorRampPalette(c("pink", "slateblue2", "purple4", "midnightblue"))
colors <-pal(6)

## PRECIPITATION

par(mar=c(1,3,2,1))
# Title precipitation
plot.new()
text(0.5,0.5,"Precipitation",cex=2,font=2)

par(mar=c(2.5,3,2,1))
# Plot methods precip
for (imet in 1:4){#4
  dat <- t(valuesP[,,imet,1]);dat2 <- t(valuesP[,,imet,2])
  dat0 <- matrix(dat,nrow = 12*7, ncol =1 );dat3 <- matrix(dat2,nrow = 12*7,ncol=1)
  dat1 <- data.frame(x = rep(seq(1, 7, 1), 12), 
                     y = rep(seq(12, 1, -1), each = 7),val = dat0)
  coordinates(dat1) <- ~ x + y
  #gridded(dat1) <- TRUE
  col0 <- rep('xx',length(dat0))
  cex0 <- array(NaN,dim = c(length(dat0)))
  
  cex0[(dat3>=0) & (dat3 < 2) ] <- 0.5
  cex0[(dat3>=2) & (dat3 < 4) ] <- 1.5
  cex0[(dat3>=4) & (dat3 < 6) ] <- 2.5
  cex0[(dat3>=6) & (dat3 < 8) ] <- 3.5
  
  col0[(dat0>= 0) & (dat0 < 0.1) ] <-  colors[1]
  col0[(dat0>= 0.1) & (dat0 < 0.2) ] <-  colors[2]
  col0[(dat0>= 0.2) & (dat0 < 0.3) ] <-  colors[3]
  col0[(dat0>= 0.3) & (dat0 < 0.4) ] <-  colors[4]
  col0[(dat0>= 0.4) & (dat0 < 0.5) ] <-  colors[5]
  
  main0 <- nmeth[imet]
  plot(dat1, pch=19,col= col0,cex=cex0,main=main0)
    #at1 = (0:(length(xLabels)-1))/(length(xLabels)-1)
    #at2 = (0:(length(yLabels)-1))/(length(yLabels)-1)
    axis(1,at=xLabels, labels=xLabels, cex.axis=1)
    axis(2,at=c(1:12), labels=yLabels, las= HORIZONTAL<-1,
         cex.axis=1)
  rm(dat);rm(dat2);rm(dat0);rm(dat1);rm(dat3)
} # End methods

#Plot legend
par(mar = c(4,4.1,2.5,1),pty='s')
d<-data.frame(x = rep(seq(1, 5, 1), 4), y = rep(seq(4, 1, -1), each = 5),val = rnorm(20))
coordinates(d) <- ~ x + y
col1<-rep(colors[1:5],4)
cex1<- rep(c(0.5,1.5,2.5,3.5),each = 5)
ylab<-c('[0,2]','[2,4]','[4,6]','[6,8]')
xlab <- c('[0,0.1]','[0.1,0.2]','[0.2,0.3]','[0.3,0.4]','[0.4,0.5]')
plot(d,col=col1,cex=cex1,pch=19,ylab='sharpness [mm/d]',xaxs = 'r',yaxs = 'i')
mtext('crpss',1,4)
axis(1,at=c(1:5),labels=xlab,las= VERTICAL<-2)
axis(2,at=c(1:4),labels=rev(ylab),las= HORIZONTAL<-1)

## TEMPERATURE

par(mar=c(1,3,2,1),pty='m')
# Title temperature
plot.new()
text(0.5,0.5,"Temperature",cex=2,font=2)

par(mar=c(2.5,3,2,1))
# Plot methods temp
for (imet in 1:4){#4
  dat <- t(valuesT[,,imet,1]);dat2 <- t(valuesT[,,imet,2])
  dat0 <- matrix(dat,nrow = 12*7, ncol =1 );dat3 <- matrix(dat2,nrow = 12*7,ncol=1)
  dat1 <- data.frame(x = rep(seq(1, 7, 1), 12), 
                     y = rep(seq(12, 1, -1), each = 7),val = dat0)
  coordinates(dat1) <- ~ x + y
  #gridded(dat1) <- TRUE
  col0 <- rep('xx',length(dat0))
  cex0 <- array(NaN,dim = c(length(dat0)))
  
  cex0[(dat3>=0) & (dat3 < 1) ] <- 0.5
  cex0[(dat3>=1) & (dat3 < 2) ] <- 1.5
  cex0[(dat3>=2) & (dat3 < 3) ] <- 2.5
  cex0[(dat3>=3) & (dat3 < 5) ] <- 3.5
  
  col0[(dat0>= 0) & (dat0 < 0.1) ] <-  colors[1]
  col0[(dat0>= 0.1) & (dat0 < 0.2) ] <-  colors[2]
  col0[(dat0>= 0.2) & (dat0 < 0.3) ] <-  colors[3]
  col0[(dat0>= 0.3) & (dat0 < 0.4) ] <-  colors[4]
  col0[(dat0>= 0.4) & (dat0 < 0.5) ] <-  colors[5]
  col0[(dat0>= 0.5) & (dat0 < 0.6) ] <-  colors[6]
  
  col0[is.na(dat0)] <- 'white'
  
  main0 <- nmeth[imet]
  plot(dat1, pch=19,col= col0,cex=cex0,main= main0)
  axis(1,at=xLabels, labels=xLabels, cex.axis=1)
  axis(2,at=c(1:12), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=1)
  rm(dat);rm(dat2);rm(dat0);rm(dat1);rm(dat3)
} # End methods

#Plot legend
par(mar = c(4.1,4,2.5,1),pty='s')
d<-data.frame(x = rep(seq(1, 6, 1), 4), y = rep(seq(4, 1, -1), each = 6),val = rnorm(24))
coordinates(d) <- ~ x + y
col1<-rep(colors,4)
cex1<- rep(c(0.5,1.5,2.5,3.5),each = 6)
ylab<-c('[0,2]','[2,4]','[4,6]','[6,8]')
xlab <- c('[0,0.1]','[0.1,0.2]','[0.2,0.3]','[0.3,0.4]','[0.4,0.5]','[0.5,0.6]')
plot(d,col=col1,cex=cex1,pch=19,ylab='sharpness [C]',xaxs = 'r',yaxs = 'i')
mtext('crpss',1,4)
axis(1,at=c(1:6),labels=xlab,las= VERTICAL<-2)
axis(2,at=c(1:4),labels=rev(ylab),las= HORIZONTAL<-1)

## REF. EVAPOTRANSPIRATION

par(mar=c(1,3,2,1),pty='m')
# Title ref. evapo
plot.new()
text(0.5,0.5,"Ref. Evapotranspiration",cex=2,font=2)

par(mar=c(2.5,3,2,1))
# Plot methods ref. evap
for (imet in 1:4){#4
  dat <- t(valuesE[,,imet,1]);dat2 <- t(valuesE[,,imet,2])
  dat0 <- matrix(dat,nrow = 12*7, ncol =1 );dat3 <- matrix(dat2,nrow = 12*7,ncol=1)
  dat1 <- data.frame(x = rep(seq(1, 7, 1), 12), 
                     y = rep(seq(12, 1, -1), each = 7),val = dat0)
  coordinates(dat1) <- ~ x + y
  #gridded(dat1) <- TRUE
  col0 <- rep('xx',length(dat0))
  cex0 <- array(NaN,dim = c(length(dat0)))
  
  cex0[(dat3>=0) & (dat3 < 0.5) ] <- 0.5
  cex0[(dat3>=0.5) & (dat3 < 1) ] <- 1.5
  cex0[(dat3>=1) & (dat3 < 1.5) ] <- 2.5
  cex0[(dat3>=1.5) & (dat3 < 2.5) ] <- 3.5
  
  col0[(dat0>= 0) & (dat0 < 0.1) ] <-  colors[1]
  col0[(dat0>= 0.1) & (dat0 < 0.2) ] <-  colors[2]
  col0[(dat0>= 0.2) & (dat0 < 0.3) ] <-  colors[3]
  col0[(dat0>= 0.3) & (dat0 < 0.4) ] <-  colors[4]
  
  col0[is.na(dat0)] <- 'white'
  
  main0 <- nmeth[imet]
  plot(dat1, pch=19,col= col0,cex=cex0,main=main0)
  axis(1,at=xLabels, labels=xLabels, cex.axis=1)
  axis(2,at=c(1:12), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=1)
  rm(dat);rm(dat2);rm(dat0);rm(dat1);rm(dat3)
} # End methods

#Plot legend
par(mar = c(5.5,4,2.5,1),pty='s')
d<-data.frame(x = rep(seq(1, 4, 1), 4), y = rep(seq(4, 1, -1), each = 4),val = rnorm(16))
coordinates(d) <- ~ x + y
col1<-rep(colors,4)
cex1<- rep(c(0.5,1.5,2.5,3.5),each = 4)
ylab<-c('[0,0.5]','[0.5,1]','[1,1.5]','[1.5,2.5]')
xlab <- c('[0,0.1]','[0.1,0.2]','[0.2,0.3]','[0.3,0.4]')
plot(d,col=col1[1:4],cex=cex1,pch=19,ylab='sharpness [mm/day]',xaxs = 'r',yaxs = 'i')
mtext('crpss',1,4)
axis(1,at=c(1:4),labels=xlab,las= VERTICAL<-2)
axis(2,at=c(1:4),labels=rev(ylab),las= HORIZONTAL<-1)

graphics.off()












# # Plotting
# 
# m <- matrix(c(seq(1:4),9,5,6,7,8,9),nrow = 2,ncol = 5,byrow = TRUE) # No Meth + 1 (for legend)
# name_plt = paste(ObsDir,'verification/',quality[iq],'.png',sep='')
# png(name_plt,width=750,height=200)
# #par(cex.lab=3)
# layout(mat = m, widths=c(2,2,2,2,0.5), heights=array(1,dim=c(5)))
# nHalf = 10;Thresh = 0
# 
# ## Make vector of colors for values below threshold
# rc1 = colorRampPalette(colors = c("red", "white"), space="Lab")(17)    
# ## Make vector of colors for values above threshold
# rc2 = colorRampPalette(colors = c("white", "blue"), space="Lab")(nHalf)
# ColorRamp = c(rc1, rc2)
# 
# edg = edg = c(-1.7,1)
# 
# 
# rb1 = seq(edg[1], Thresh, length.out=17+1)#nHalf+1
# rb2 = seq(Thresh, edg[2], length.out=nHalf+1)[-1]
# ColorLevels = c(rb1, rb2)
# 
# # Plot methods precip
# for (imet in 1:4){
#   dat <- valuesP[,,imet]
#   dat <- dat[reverse,]
#   main0 = meth[imet]
#   image(t(dat),col=ColorRamp, axes=FALSE, zlim=c(edg[1],edg[2]),main = main0,xlab = 'lt (months)')
#   at1 = (0:(length(xLabels)-1))/(length(xLabels)-1)
#   at2 = (0:(length(yLabels)-1))/(length(yLabels)-1)
#   axis(1,at=at1, labels=xLabels, cex.axis=1)
#   axis(2,at=at2, labels=yLabels, las= HORIZONTAL<-1,
#        cex.axis=1)
# } # End methods
# 
# # Plot methods temp
# for (imet in 1:4){
#   dat <- valuesT[,,imet]
#   dat <- dat[reverse,]
#   main0 = meth[imet]
#   image(t(dat),col=ColorRamp, axes=FALSE, zlim=c(edg[1],edg[2]),main = main0,xlab = 'lt (months)')
#   at1 = (0:(length(xLabels)-1))/(length(xLabels)-1)
#   at2 = (0:(length(yLabels)-1))/(length(yLabels)-1)
#   axis(1,at=at1, labels=xLabels, cex.axis=1)
#   axis(2,at=at2, labels=yLabels, las= HORIZONTAL<-1,
#        cex.axis=1)
# } # End methods
# 
# #Plot legend
# par(mar = c(3,2,2.5,1))
# image(1, ColorLevels,
#       matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
#       col=ColorRamp,
#       xlab="",ylab="",
#       xaxt="n",axes = FALSE)
# axis(2,at=ColorLevels, labels=ColorLevels, las= HORIZONTAL<-1,
#      cex.axis=1)
# 
# graphics.off()
## Makes heat map with crpss_clim and sharpness

## Makes Heat Maps of Bias for PCP, TEMP and ET0

quality <- c('crpss_clim','sharp')
meth = c('R','dc','qm','bma')
nmeth <- c('RAW','LS','QM','BMA')
# Load coord files
ObsDir = 'G:/calibMET/'
Nmonth = 12;Nlead = 7
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

reverse = seq(12,1,-1)
yLabels <- Mname[reverse]
xLabels <- seq(1,7,1)

valuesP = array(NaN,dim=c(Nmonth,662,Nlead,length(meth),length(quality)))
valuesT = array(NaN,dim=c(Nmonth,724,Nlead,length(meth),length(quality)))
valuesE = array(NaN,dim=c(Nmonth,724,Nlead,length(meth),length(quality)))
for (imeth in 1:length(meth)){#length(meth)
  for (nmonth in 1:Nmonth){
    ## Precipitation 
    name_file = paste(ObsDir,'precip/verification_daily/',meth[imeth],'_crpss_clim_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=662))
    valuesP[nmonth,,,imeth,1] = val
    name_file = paste(ObsDir,'precip/verification_daily/',meth[imeth],'_sharp_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=662))
    valuesP[nmonth,,,imeth,2] = val
    ## Temperature
    name_file = paste(ObsDir,'temp/verification_daily/',meth[imeth],'_crpss_clim_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=724))
    valuesT[nmonth,,,imeth,1] = val
    name_file = paste(ObsDir,'temp/verification_daily/',meth[imeth],'_sharp_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=724))
    valuesT[nmonth,,,imeth,2] = val
    ## Ref. Evapotranspiration
    name_file = paste(ObsDir,'ref_evap/verification_daily/',meth[imeth],'_crpss_clim_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=724))
    valuesE[nmonth,,,imeth,1] = val
    name_file = paste(ObsDir,'ref_evap/verification_daily/',meth[imeth],'_sharp_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=724))
    valuesE[nmonth,,,imeth,2] = val
  } # end loop month
}

# Plotting

nmonth <- 2
nlead <- 1

c0 <- c(rep(1,5),seq(2,6),rep(7,5),seq(8,12),rep(13,5),seq(14,18)) 
m <- matrix(c0,nrow = 6,ncol = 5,byrow = TRUE) # No Meth + 1 (for legend)
name_plt = paste0(ObsDir,'PQ_figures/fig5/fig5_m',nmonth,'_lt',nlead,'.png')
png(name_plt,width=750,height=750)
layout(mat = m,widths = c(2,2,2,2,2.5),heights=array(c(0.3,1,0.3,1,0.3,1)))# 

colors <- c('burlywood','darkgoldenrod1','firebrick1','darkmagenta','darkgreen','brown4')


## Plotting precipitation
# Coordinates
name_f = paste0(ObsDir,'PQ_figures/fig5/coordP.txt')
coorP   = t(matrix(scan(name_f),nrow=3,ncol=662))

# Title
par(mar=c(1,3,2,1))
plot.new()
text(0.5,0.5,"Precipitation",cex=2,font=2)

# Maps
par(mar=c(2.5,3,2,1))
for (imet in 1:length(meth)){
  data1 <- valuesP[nmonth,,nlead,imet,1]
  data2 <- valuesP[nmonth,,nlead,imet,2]

  col0 <- rep('xx',length(data1))
  cex0 <- array(NaN,dim = c(length(data1)))

  ## Sharpness
  cex0[(data2>=0) & (data2 < 2) ] <- 0.5
  cex0[(data2>=2) & (data2 < 4) ] <- 1
  cex0[(data2>=4) & (data2 < 6) ] <- 1.5
  cex0[(data2>=6) & (data2 < 8) ] <- 2
  cex0[(data2>=8) ] <- 2.5
  ## CRPSS
  col0[(data1 < 0) ] <-  colors[1]
  col0[(data1>= 0) & (data1 < 0.1) ] <-  colors[2]
  col0[(data1>= 0.1) & (data1 < 0.2) ] <-  colors[3]
  col0[(data1>= 0.2) & (data1 < 0.3) ] <-  colors[4]
  col0[(data1>= 0.3) & (data1 < 0.4) ] <-  colors[5]
  col0[(data1>= 0.4) & (data1 < 0.5) ] <-  colors[6]

  grd1 <- data.frame(z = data1, xc = coorP[,2], yc = coorP[,3])
  coordinates(grd1) <- ~ xc + yc
  main0 = nmeth[imet]
  plot(grd1, pch=19, col= col0,cex=cex0,main = main0)
  rm(grd1);rm(col0);rm(cex0)
}

# Legend
par(mar = c(5.1,4,2.5,1),pty='s')
d<-data.frame(x = rep(seq(1, 6, 1), 5), y = rep(seq(5, 1, -1), each = 6),val = rnorm(30))
coordinates(d) <- ~ x + y
col1<-rep(colors,5)
cex1<- rep(c(0.5,1,1.5,2,2.5),each = 6)
ylab<-c('[ 0,2 ]','[ 2,4 ]','[ 4,6 ]','[ 6,8 ]','[ >8 ]')
xlab <- c('[ <0 ]','[ 0,0.1 ]','[ 0.1,0.2 ]','[ 0.2,0.3 ]','[ 0.3,0.4 ]','[ 0.4,0.5 ]')
plot(d,col=col1,cex=cex1,pch=19,ylab='sharpness [mm/day]',xaxs = 'r',yaxs = 'i',
     main = 'Legend')
mtext('crpss',1,4)
axis(1,at=c(1:6),labels=xlab,las =2)
axis(2,at=c(1:5),labels=rev(ylab),las= HORIZONTAL<-1)

## Plotting temperature
# Coordinates
name_f = paste0(ObsDir,'PQ_figures/fig5/coordT.txt')
coorT   = t(matrix(scan(name_f),nrow=3,ncol=724))

# Title
par(mar=c(1,3,2,1),pty='m')
plot.new()
text(0.5,0.5,"Temperature",cex=2,font=2)

# Maps
par(mar=c(2.5,3,2,1))
for (imet in 1:length(meth)){
  data1 <- valuesT[nmonth,,nlead,imet,1]
  data2 <- valuesT[nmonth,,nlead,imet,2]
  
  col0 <- rep('xx',length(data1))
  cex0 <- array(NaN,dim = c(length(data1)))
  
  ## Sharpness
  cex0[(data2>=0) & (data2 < 1) ] <- 0.5
  cex0[(data2>=1) & (data2 < 2) ] <- 1
  cex0[(data2>=2) & (data2 < 3) ] <- 1.5
  cex0[(data2>=3) & (data2 < 5) ] <- 2
  cex0[(data2>=5) ] <- 2.5
  ## CRPSS
  col0[(data1 < 0) ] <-  colors[1]
  col0[(data1>= 0) & (data1 < 0.1) ] <-  colors[2]
  col0[(data1>= 0.1) & (data1 < 0.2) ] <-  colors[3]
  col0[(data1>= 0.2) & (data1 < 0.3) ] <-  colors[4]
  col0[(data1>= 0.3) & (data1 < 0.4) ] <-  colors[5]
  col0[(data1>= 0.4) & (data1 < 0.6) ] <-  colors[6]
  
  col0[is.na(data1)]<-'white'
  
  grd1 <- data.frame(z = data1, xc = coorT[,2], yc = coorT[,3])
  coordinates(grd1) <- ~ xc + yc
  main0 = nmeth[imet]
  plot(grd1, pch=19, col= col0,cex=cex0,main = main0)
}

# Legend
par(mar = c(5.1,4,2.5,1),pty='s')
d<-data.frame(x = rep(seq(1, 6, 1), 5), y = rep(seq(5, 1, -1), each = 6),val = rnorm(30))
coordinates(d) <- ~ x + y
col1<-rep(colors,5)
cex1<- rep(c(0.5,1,1.5,2,2.5),each = 6)
ylab<-c('[0,1]','[1,2]','[2,3]','[3,5]','[>5]')
xlab <- c('[<0]','[0,0.1]','[0.1,0.2]','[0.2,0.3]','[0.3,0.4]','[0.4,0.6]')
plot(d,col=col1,cex=cex1,pch=19,ylab='sharpness [C]',xaxs = 'r',yaxs = 'i',
     main = 'Legend')
mtext('crpss',1,4)
axis(1,at=c(1:6),labels=xlab,las =2)
axis(2,at=c(1:5),labels=rev(ylab),las= HORIZONTAL<-1)

## Plotting ref. evapotranspiration
# Title
par(mar=c(1,3,2,1),pty='m')
plot.new()
text(0.5,0.5,"Ref. Evapotranspiration",cex=2,font=2)

# Maps
par(mar=c(2.5,3,2,1))
for (imet in 1:length(meth)){
  data1 <- valuesE[nmonth,,nlead,imet,1]
  data2 <- valuesE[nmonth,,nlead,imet,2]
  
  col0 <- rep('xx',length(data1))
  cex0 <- array(NaN,dim = c(length(data1)))
  
  ## Sharpness
  cex0[(data2>=0) & (data2 < 1) ] <- 0.5
  cex0[(data2>=1) & (data2 < 2) ] <- 1
  cex0[(data2>=2) & (data2 < 3) ] <- 1.5

  ## CRPSS
  col0[(data1 < 0) ] <-  colors[1]
  col0[(data1>= 0) & (data1 < 0.1) ] <-  colors[2]
  col0[(data1>= 0.1) & (data1 < 0.2) ] <-  colors[3]
  col0[(data1>= 0.2) & (data1 < 0.3) ] <-  colors[4]
  col0[(data1>= 0.3) & (data1 < 0.4) ] <-  colors[5]
  col0[(data1>= 0.4) & (data1 < 0.6) ] <-  colors[6]
  
  col0[is.na(data1)]<-'white'
  grd1 <- data.frame(z = data1, xc = coorT[,2], yc = coorT[,3])
  coordinates(grd1) <- ~ xc + yc
  main0 = nmeth[imet]
  plot(grd1, pch=19, col= col0,cex=cex0,main = main0)
}

# Legend
par(mar = c(5.1,4,2.5,1),pty='s')
d<-data.frame(x = rep(seq(1, 6, 1), 3), y = rep(seq(3, 1, -1), each = 6),val = rnorm(18))
coordinates(d) <- ~ x + y
col1<-rep(colors,3)
cex1<- rep(c(0.5,1,1.5),each = 6)
ylab<-c('[0,1]','[1,2]','[2,3]')
xlab <- c('[<0]','[0,0.1]','[0.1,0.2]','[0.2,0.3]','[0.3,0.4]','[0.4,0.6]')
plot(d,col=col1,cex=cex1,pch=19,ylab='sharpness [mm/day]',xaxs = 'r',yaxs = 'i',
     main = 'Legend')
mtext('crpss',1,4)
axis(1,at=c(1:6),labels=xlab,las =2)
axis(2,at=c(1:3),labels=rev(ylab),las= HORIZONTAL<-1)

graphics.off()
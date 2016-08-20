## Makes Heat Maps of Bias for PCP, TEMP and ET0

quality <- 'bias'
meth = c('R','dc','qm','bma')
nmeth <- c('RAW','LS','QM','BMA')
# Load coord files
ObsDir = 'G:/calibMET/'
Nmonth = 12;Nlead = 7
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

reverse = seq(12,1,-1)
yLabels <- Mname[reverse]
xLabels <- seq(1,7,1)

valuesP = array(NaN,dim=c(Nmonth,Nlead,length(meth)))
valuesT = array(NaN,dim=c(Nmonth,Nlead,length(meth)))
valuesE = array(NaN,dim=c(Nmonth,Nlead,length(meth)))
for (imeth in 1:length(meth)){#length(meth)
  for (nmonth in 1:Nmonth){
    ## Precipitation
    name_file = paste(ObsDir,'precip/verification_daily/',meth[imeth],'_',quality,'_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=662))
    valuesP[nmonth,,imeth] = colMeans(val,na.rm  = TRUE)
    ## Temperature
    name_file = paste(ObsDir,'temp/verification_daily/',meth[imeth],'_',quality,'_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=724))
    valuesT[nmonth,,imeth] = colMeans(val,na.rm  = TRUE)
    ## Ref Evap
    name_file = paste(ObsDir,'ref_evap/verification_daily/',meth[imeth],'_',quality,'_m',nmonth,'.txt',sep="")
    val = t(matrix(scan(name_file),nrow=7,ncol=724))
    valuesE[nmonth,,imeth] = colMeans(val,na.rm  = TRUE)
  } # end loop month
}
  
# Plotting

c0 <- c(rep(1,5),seq(2,6),rep(7,5),seq(8,12),rep(13,5),seq(14,18))    
m <- matrix(c0,nrow = 6,ncol = 5,byrow = TRUE) # No Meth + 1 (for legend)
name_plt = paste(ObsDir,'PQ_figures/fig1/',quality,'.png',sep='')
png(name_plt,width=750,height=750)
#par(cex.lab=3)
layout(mat = m, widths=c(2,2,2,2,0.5), heights=array(c(0.3,1,0.3,1,0.3,1)))


## PRECIPITATION
# Colors precipitation
nHalf = 32;Thresh = 0
## Make vector of colors for values below threshold
rc1 = colorRampPalette(colors = c("red", "white"), space="Lab")(12)    
## Make vector of colors for values above threshold
rc2 = colorRampPalette(colors = c("white", "blue"), space="Lab")(nHalf)
ColorRamp = c(rc1, rc2)
edg = c(-60,160)
rb1 = seq(edg[1], Thresh, length.out=12+1)#nHalf+1
rb2 = seq(Thresh, edg[2], length.out=nHalf+1)[-1]
ColorLevels = c(rb1, rb2)

par(mar=c(1,3,2,1))
# Title precipitation
plot.new()
text(0.5,0.5,"Precipitation",cex=2,font=2)

par(mar=c(2.5,3,2,1))
# Plot methods precip
for (imet in 1:4){
  dat <- valuesP[,,imet]
  dat <- dat[reverse,]
  main0 = nmeth[imet]
  image(t(dat),col=ColorRamp, axes=FALSE, zlim=c(edg[1],edg[2]),main = main0)
  at1 = (0:(length(xLabels)-1))/(length(xLabels)-1)
  at2 = (0:(length(yLabels)-1))/(length(yLabels)-1)
  axis(1,at=at1, labels=xLabels, cex.axis=1)
  axis(2,at=at2, labels=yLabels, las= HORIZONTAL<-1,
   cex.axis=1)
} # End methods

#Plot legend precipitation
#par(mar = c(3,2,2.5,1))
image(1, ColorLevels[-45],
      matrix(data=ColorLevels[-45], ncol=length(ColorLevels[-45]),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n",axes = FALSE,main='%')
at0 <- seq(-62.5,157.5,10); labels0 <- seq(-60,160,10)
#at=(ColorLevels-2.5), labels=ColorLevels
axis(2,at=at0, labels=labels0, las= HORIZONTAL<-1,
     cex.axis=0.9)


## TEMPERATURE
# Colors temperature
nHalf = 18;Thresh = 0
## Make vector of colors for values below threshold
rc1 = colorRampPalette(colors = c("red", "white"), space="Lab")(6)    
## Make vector of colors for values above threshold
rc2 = colorRampPalette(colors = c("white", "blue"), space="Lab")(nHalf)
ColorRamp = c(rc1, rc2)
edg = c(-30,90)
rb1 = seq(edg[1], Thresh, length.out=6+1)#nHalf+1
rb2 = seq(Thresh, edg[2], length.out=nHalf+1)[-1]
ColorLevels = c(rb1, rb2)


par(mar=c(1,3,2,1))
# Title temperature
plot.new()
text(0.5,0.5,"Temperature",cex=2,font=2)

par(mar=c(2.5,3,2,1))
# Plot methods temp
for (imet in 1:4){
  dat <- valuesT[,,imet]
  dat <- dat[reverse,]
  main0 = nmeth[imet]
  image(t(dat),col=ColorRamp, axes=FALSE, zlim=c(edg[1],edg[2]),xlab = 'lt (months)',main = main0)
  at1 = (0:(length(xLabels)-1))/(length(xLabels)-1)
  at2 = (0:(length(yLabels)-1))/(length(yLabels)-1)
  axis(1,at=at1, labels=xLabels, cex.axis=1)
  axis(2,at=at2, labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=1)
} # End methods

#Plot legend temperature
#par(mar = c(3,2,2.5,1))
image(1, ColorLevels[-25],
  matrix(data=ColorLevels[-25], ncol=length(ColorLevels[-25]),nrow=1),
  col=ColorRamp,
  xlab="",ylab="",
  xaxt="n",axes = FALSE,main = '%')
at0 <- seq(-32.5,87.5,10); labels0 <- seq(-30,90,10)
#at=(ColorLevels-5), labels=ColorLevels
axis(2,at=at0, labels=labels0, las= HORIZONTAL<-1,
           cex.axis=1)

## REF EVAP
# Colors ref evapotranspiration
nHalf = 20;Thresh = 0
## Make vector of colors for values below threshold
rc1 = colorRampPalette(colors = c("red", "white"), space="Lab")(4)    
## Make vector of colors for values above threshold
rc2 = colorRampPalette(colors = c("white", "blue"), space="Lab")(nHalf)
ColorRamp = c(rc1, rc2)
edg = c(-20,100)
rb1 = seq(edg[1], Thresh, length.out=4+1)#nHalf+1
rb2 = seq(Thresh, edg[2], length.out=nHalf+1)[-1]
ColorLevels = c(rb1, rb2)


par(mar=c(1,3,2,1))
# Title ref evapo
plot.new()
text(0.5,0.5,"Ref. Evapotranspiration",cex=2,font=2)

par(mar=c(2.5,3,2,1))
# Plot methods temp
for (imet in 1:4){
  dat <- valuesE[,,imet]
  dat <- dat[reverse,]
  main0 = nmeth[imet]
  image(t(dat),col=ColorRamp, axes=FALSE, zlim=c(edg[1],edg[2]),xlab = 'lt (months)',main = main0)
  at1 = (0:(length(xLabels)-1))/(length(xLabels)-1)
  at2 = (0:(length(yLabels)-1))/(length(yLabels)-1)
  axis(1,at=at1, labels=xLabels, cex.axis=1)
  axis(2,at=at2, labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=1)
} # End methods

#Plot legend temperature
#par(mar = c(3,2,2.5,1))
image(1, ColorLevels[-25],
      matrix(data=ColorLevels[-25], ncol=length(ColorLevels[-25]),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n",axes = FALSE,main = '%')
at0 <- seq(-22.5,97.5,10); labels0 <- seq(-20,100,10)
#at=(ColorLevels-5), labels=ColorLevels
axis(2,at=at0, labels=labels0, las= HORIZONTAL<-1,
     cex.axis=1)

graphics.off()
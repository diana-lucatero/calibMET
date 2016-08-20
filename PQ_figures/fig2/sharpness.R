## Plot sharpness

quality = 'sharp'
meth = c('R','dc','qm','bma')
# Load coord files
ObsDir = 'E:/calibMET/precip/'
Nmonth = 12;Nlead = 7
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
Nsam = c(25,24,25,25,25,24,24,23,24,24,24,24)

reverse = seq(12,1,-1)
yLabels <- Mname[reverse]
xLabels <- seq(1,7,1)

values = array(NaN,dim=c(Nmonth,Nlead,length(meth)))
for (imeth in 1:length(meth)){#length(meth)
  for (lt in (1:7)){
    for (nmonth in 1:Nmonth){#Nmonth
      name_file = paste(ObsDir,'verification/',meth[imeth],'_',quality,'_m',nmonth,'_lt',lt,'.txt',sep="")
      val = t(matrix(scan(name_file),nrow=Nsam[nmonth],ncol=662))
      values[nmonth,lt,imeth] = mean(val,na.rm  = TRUE)
    }
  } # end loop month
}

#Plotting

m <- matrix(seq(1:5),nrow = 1,ncol = 5,byrow = TRUE) # No Meth + 1 (for legend)
name_plt = paste(ObsDir,'verification/',quality,'.png',sep='')
png(name_plt,width=750,height=200)
#par(cex.lab=3)
layout(mat = m, widths=c(2,2,2,2,0.5), heights=array(1,dim=c(5)))
nHalf = 11

## Make vector of colors for values below threshold
#rc1 = colorRampPalette(colors = c("red", "white"), space="Lab")(nHalf)    
## Make vector of colors for values above threshold
ColorRamp = colorRampPalette(colors = c("white", "blue"), space="Lab")(nHalf)
#ColorRamp = c(rc1, rc2)
## In your example, this line sets the color for values between 49 and 51. 
#ColorRamp[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 

#min = min(values,na.rm = TRUE);max=max(values,na.rm = TRUE)

ColorLevels = seq(0, 110, length=length(ColorRamp))
#rb2 = seq(Thresh, max, length.out=nHalf+1)[-1]
#ColorLevels = c(rb1, rb2)

# Plot methods
for (imet in 1:4){
  dat <- values[,,imet]
  #reverse <- nrow(dat) : 1
  #yLabels <- yLabels[reverse]
  dat <- dat[reverse,]
  main0 = meth[imet]
  #par(mar = c(3,2,2.5,1))
  image(t(dat),col=ColorRamp , 
        ylab="", axes=FALSE, zlim=c(0,110),main = main0,xlab = 'lt (months)')
  #image(t(dat),col=ColorRamp, zlim=c(min,max))
  at1 = (0:(length(xLabels)-1))/(length(xLabels)-1)
  at2 = (0:(length(yLabels)-1))/(length(yLabels)-1)
  axis(1,at=at1, labels=xLabels, cex.axis=1)
  axis(2,at=at2, labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=1)
}
#Plot legend
par(mar = c(3,2,2.5,1))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      ylab="",
      xaxt="n",axes = FALSE)
axis(2,at=ColorLevels, labels=ColorLevels, las= HORIZONTAL<-1,
     cex.axis=1)

graphics.off()
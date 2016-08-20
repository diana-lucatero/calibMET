## Plot sharpness

quality = 'sharp'
meth = c('R','dc','qm')
# Load coord files
ObsDir = 'E:/calibMET/precip/'
Nmonth = 12;Nlead = 7
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
Nsam = c(25,24,25,25,25,24,24,23,24,24,24,24)

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

m <- matrix(seq(1:4),nrow = 1,ncol = 4,byrow = TRUE) # No Meth + 1 (for legend)
name_plt = paste(ObsDir,'verification/',quality,'.png',sep='')
png(name_plt,width=750,height=200)
#par(cex.lab=3)
layout(mat = m, widths=c(2,2,2,0.5), heights=array(1,dim=c(4)))
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
for (imet in 1:3){
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

graphics.off()
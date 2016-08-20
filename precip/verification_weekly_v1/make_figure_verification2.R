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
  
  reverse = seq(12,1,-1)
  yLabels <- Mname[reverse]
  xLabels <- seq(1,Nlead,1)
  
  for (iq in 1:length(quality)){#length(quality)
    values = array(NaN,dim=c(Nmonth,Nlead,length(meth)))
    for (imeth in 1:length(meth)){#length(meth)
      for (nmonth in 1:Nmonth){
        name_file = paste(ObsDir,'verification_weekly/',meth[imeth],'_',quality[iq],'_m',nmonth,'.txt',sep="")
        val = t(matrix(scan(name_file),nrow=Nlead,ncol=662));val[which(val==-100)]=NaN
        values[nmonth,,imeth] = colMeans(val,na.rm  = TRUE)
      } # end loop month
    }
  
    #Plotting
  
    m <- matrix(seq(1:5),nrow = 1,ncol = 5,byrow = TRUE) # No Meth + 1 (for legend)
    name_plt = paste(ObsDir,'verification_weekly/',quality[iq],'.png',sep='')
    png(name_plt,width=750,height=200)
    #par(cex.lab=3)
    layout(mat = m, widths=c(2,2,2,2,0.5), heights=array(1,dim=c(5)))
    nHalf = 10;Thresh = 0
    
    ## Make vector of colors for values below threshold
    rc1 = colorRampPalette(colors = c("red", "white"), space="Lab")(nHalf)    
    ## Make vector of colors for values above threshold
    rc2 = colorRampPalette(colors = c("white", "blue"), space="Lab")(nHalf)
    ColorRamp = c(rc1, rc2)
    ## In your example, this line sets the color for values between 49 and 51. 
    #ColorRamp[c(nHalf, nHalf+1)] = rgb(t(col2rgb("white")), maxColorValue=256) 
    
    if (iq == 1) {edg = c(-100,100)}
      else if (iq == 2) {edg = c(-1,1)}
      else edg = c(-1,1)
    #min = min(values,na.rm = TRUE);max=max(values,na.rm = TRUE)
    #ColorLevels <- seq(min, max, length=length(ColorRamp))
    
    rb1 = seq(edg[1], Thresh, length.out=nHalf+1)
    rb2 = seq(Thresh, edg[2], length.out=nHalf+1)[-1]
    ColorLevels = c(rb1, rb2)
  
    # Plot methods
    for (imet in 1:4){
      dat <- values[,,imet]
      dat <- dat[reverse,]
      main0 = meth[imet]
      image(t(dat),col=ColorRamp, axes=FALSE, zlim=c(edg[1],edg[2]),main = main0,xlab = 'lt (weeks)')
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
        xlab="",ylab="",
        xaxt="n",axes = FALSE)
    at2 = (0:(length(ColorLevels)-1))/(length(ColorLevels)-1)
    axis(2,at=ColorLevels, labels=ColorLevels, las= HORIZONTAL<-1,
         cex.axis=1)
    
  } # End qualities
  graphics.off()
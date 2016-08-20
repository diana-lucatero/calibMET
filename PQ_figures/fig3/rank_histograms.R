## Rank Histograms

quality = 'RH'
meth = c('R','dc','qm','bma')
nmeth <- c('RAW','LS','QM','BMA')
# Load coord files
ObsDir = 'G:/calibMET/'
Nmonth = 12;Nlead = 7
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)

nmonth = 2

for (lt in 1:Nlead){
  arg0 = seq(1,(Nens[nmonth]+1),1)
  valuesP = array(NaN,dim=c(Nens[nmonth]+1,length(meth)))
  valuesT = array(NaN,dim=c(Nens[nmonth]+1,length(meth)))
  valuesE = array(NaN,dim=c(Nens[nmonth]+1,length(meth)))
  for (imeth in 1:length(meth)){#length(meth)
      # PRECIPITATION
      name_file = paste(ObsDir,'precip/verification_daily/',meth[imeth],'_',quality,'_m',nmonth,'_lt',lt,'.txt',sep="")
      val = t(matrix(scan(name_file),nrow=(Nens[nmonth]+1),ncol=662))
      valuesP[,imeth] = colSums(val,na.rm  = TRUE)
      # TEMPERATURE
      name_file = paste(ObsDir,'temp/verification_daily/',meth[imeth],'_',quality,'_m',nmonth,'_lt',lt,'.txt',sep="")
      val = t(matrix(scan(name_file),nrow=(Nens[nmonth]+1),ncol=724))
      valuesT[,imeth] = colSums(val,na.rm  = TRUE)
      # REF EVAPOTRANSOIRATION
      name_file = paste(ObsDir,'ref_evap/verification_daily/',meth[imeth],'_',quality,'_m',nmonth,'_lt',lt,'.txt',sep="")
      val = t(matrix(scan(name_file),nrow=(Nens[nmonth]+1),ncol=724))
      valuesE[,imeth] = colSums(val,na.rm  = TRUE)
  }

  ylim0 <- c(0,300000)
  ylim1 <- c(0,90000)
  ## START PLOTTING
  c0 <- c(rep(1,4),seq(2,5),rep(6,4),seq(7,10),rep(11,4),seq(12,15))    
  m <- matrix(c0,nrow = 6,ncol = 4,byrow = T) # No Meth + 1 (for legend)

  name_plt = paste0(ObsDir,'PQ_figures/fig3/fig3_m',nmonth,'_lt',lt,'.png')
  png(name_plt,width=750,height=750)
  layout(mat = m,heights = c(0.3,1,0.3,1,0.3,1))#, widths=c(2,2,2,0.5), heights=array(1,dim=c(4)))
  #par(cex.axis = 2)

  # Title precipitation
  par(mar=c(1,3,2,1))
  plot.new()
  text(0.5,0.5,"Precipitation",cex=2,font=2)

  par(mar=c(4,5,2,1))
  # PRECIPITATION
  for (imeth in 1:length(meth)){#length(meth)
    main0 = nmeth[imeth]
    barplot(valuesP[,imeth],names.arg = arg0,main=main0,ylim = ylim0,xlab = 'ens #',ylab = 'frequency',cex.axis = 1.5,cex.lab=1.5)#,
  }

  # Title temperature
  par(mar=c(1,3,2,1))
  plot.new()
  text(0.5,0.5,"Temperature",cex=2,font=2)

  par(mar=c(4,5,2,1))
  # TEMPERATURE
  for (imeth in 1:length(meth)){#length(meth)
    main0 = nmeth[imeth]
    barplot(valuesT[,imeth],names.arg = arg0,main=main0,ylim = ylim1,xlab = 'ens #',ylab = 'frequency',cex.axis = 1.5,cex.lab=1.5)#,ylim = ylim0
  }
  
  # Title ref evapotranspiration
  par(mar=c(1,3,2,1))
  plot.new()
  text(0.5,0.5,"Ref. Evapotranspiration",cex=2,font=2)
  
  par(mar=c(4,5,2,1))
  # REF. EVAPOTRANSPIRATION
  for (imeth in 1:length(meth)){#length(meth)
    main0 = nmeth[imeth]
    barplot(valuesE[,imeth],names.arg = arg0,main=main0,ylim = ylim1,xlab = 'ens #',ylab = 'frequency',cex.axis = 1.5,cex.lab=1.5)#,ylim = ylim0
  }
  rm(valuesP);rm(valuesT);rm(valuesE)
} # End lead time
graphics.off()
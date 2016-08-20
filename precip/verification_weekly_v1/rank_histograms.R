## Make plots of rank histograms and sharpness

meth = c('R','dc','qm','bma')
# Load coord files
ObsDir = 'E:/calibMET/precip/'
Nmonth = 12;Nlead = 30
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
Nsam = c(25,24,25,25,25,24,24,23,24,24,24,24)

m <- matrix(seq(1:8),nrow = 2,ncol = 4,byrow = FALSE) # No Meth + 1 (for legend)
breaks0 = seq(0,10,0.5)
for (nmonth in 1:1){#Nmonth
  for (lt in (1:Nlead)){#Nlead
    print(lt)
    values0 = array(NaN,dim=c(length(meth),662,(Nens[nmonth]+1)))
    values1 = array(NaN,dim=c(length(meth),662,Nsam[nmonth]))
    for (imeth in 1:length(meth)){#length(meth)
      ## Rank histogram
      name_file = paste(ObsDir,'verification_weekly/',meth[imeth],'_RH_m',nmonth,'_lt',lt,'.txt',sep="")
      val = t(matrix(scan(name_file),nrow=(Nens[nmonth]+1),ncol=662))
      values0[imeth,,] = val
      ## Sharpness
      name_file = paste(ObsDir,'verification_weekly/',meth[imeth],'_sharp_m',nmonth,'_lt',lt,'.txt',sep="")
      val = t(matrix(scan(name_file),nrow=Nsam[nmonth],ncol=662))
      values1[imeth,,] = val
    } # End methods read
    ## Plot
    for (ngrid in (1:662)){#662
      name_plt = paste(ObsDir,'verification_weekly/plots/g',ngrid,'_m',nmonth,'_lt',lt,'.png',sep="")
      png(name_plt,width=750,height=350)
      layout(mat = m)#, widths=c(2,2,2,0.5), heights=array(1,dim=c(4)))
      # Plot methods
      for (imet in 1:4){
        main0 = meth[imet]; arg0 = seq(1,(Nens[nmonth]+1),1)
        barplot(values0[imet,ngrid,],names.arg = arg0,main=main0,ylim  = c(0,10))
        #par(mar = c(1,0,1,1))
        hist(values1[imet,ngrid,],main = '',xlab = '',ylab = '',ylim  = c(0,10))
      } # End methods plot
      graphics.off()
    } # End grid plot
  } # End lead time
} # end month
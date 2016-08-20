## Rank Histograms

quality = 'RH'
meth = c('R','dc','qm','bma')
# Load coord files
ObsDir = 'E:/calibMET/temp/'
Nmonth = 12;Nlead = 7
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
m <- matrix(seq(1:4),nrow = 1,ncol = 4,byrow = TRUE) # No Meth + 1 (for legend)

lt = 1;nmonth = 1
par(cex.axis = 2)
name_plt = paste(ObsDir,'verification/','rankhist_month',nmonth,'_lt',lt,'.png',sep="")
png(name_plt,width=750,height=250)
layout(mat = m)#, widths=c(2,2,2,0.5), heights=array(1,dim=c(4)))

ylim0 <- c(0,5000)
arg0 = seq(1,(Nens[nmonth]+1),1)    
#values = array(NaN,dim=c(Nmonth,Nlead,length(meth)))
for (imeth in 1:length(meth)){#length(meth)
      name_file = paste(ObsDir,'verification/',meth[imeth],'_',quality,'_m',nmonth,'_lt',lt,'.txt',sep="")
      val = t(matrix(scan(name_file),nrow=(Nens[nmonth]+1),ncol=724))
      values = colSums(val,na.rm  = TRUE)
      main0 = meth[imeth]
      barplot(values,names.arg = arg0,main=main0,ylim = ylim0)
}
graphics.off()
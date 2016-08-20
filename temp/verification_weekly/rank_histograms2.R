## Rank Histograms

quality = 'RH'
meth = c('R','dc','qm','bma')
# Load coord files
ObsDir = 'E:/calibMET/temp/'
Nmonth = 12;Nlead = 30
Mname = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
m <- matrix(seq(1:12),nrow = 3,ncol = 4,byrow = TRUE) # No Meth + 1 (for legend)

#values = array(NaN,dim=c(Nmonth,Nlead,length(meth)))
for (imeth in 1:length(meth)){#length(meth)
  for (lt in (1:Nlead)){
    name_plt = paste(ObsDir,'verification_weekly/',meth[imeth],'_',quality,'_lt',lt,'.png',sep="")
    png(name_plt,width=750,height=750)
    layout(mat = m)#, widths=c(2,2,2,0.5), heights=array(1,dim=c(4)))
    ## Start plotting
    for (nmonth in 1:Nmonth){#Nmonth
      name_file = paste(ObsDir,'verification_weekly/',meth[imeth],'_',quality,'_m',nmonth,'_lt',lt,'.txt',sep="")
      val = t(matrix(scan(name_file),nrow=(Nens[nmonth]+1),ncol=724))
      values = colSums(val,na.rm  = TRUE)
      main0 = paste(Mname[nmonth],' lt ',lt,sep=''); arg0 = seq(1,(Nens[nmonth]+1),1)
      barplot(values,names.arg = arg0,main=main0)
    }
    graphics.off()
  } # end loop month
}

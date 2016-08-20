## Conditional verification figure

Dir <- 'G:/calibMET/'
meth <- c('R','dc','qm','bma')
nmonth <- 2; Nlead = 7
## Precipitation
pcp <- array(NaN,dim = c(662,Nlead,length(meth),5))
temp <- array(NaN,dim = c(724,Nlead,length(meth),5))
evap <- array(NaN,dim = c(724,Nlead,length(meth),5))
for (imet in 1:length(meth)){#length(meth)
  for (ithre in 1:5){#5
    name_file <- paste0(Dir,'precip/verification_daily/',meth[imet],'_crpss_clim_m',nmonth,'_thre',ithre,'.txt')
    pcp[,,imet,ithre] <- t(matrix(scan(name_file),ncol = 662,nrow = 7))
    name_file <- paste0(Dir,'temp/verification_daily/',meth[imet],'_crpss_clim_m',nmonth,'_thre',ithre,'.txt')
    temp[,,imet,ithre] <- t(matrix(scan(name_file),ncol = 724,nrow = 7))
    name_file <- paste0(Dir,'ref_evap/verification_daily/',meth[imet],'_crpss_clim_m',nmonth,'_thre',ithre,'.txt')
    evap[,,imet,ithre] <- t(matrix(scan(name_file),ncol = 724,nrow = 7))
  }
}

## Start plotting

ilead <- 1

# Data
pcp0 <- pcp[,ilead,,]
pcp1 <- data.frame(matrix(pcp0,nrow=662,ncol = 20))
temp0 <- temp[,ilead,,]
temp1 <- data.frame(matrix(temp0,nrow=724,ncol = 20))
evap0 <- evap[,ilead,,]
evap1 <- data.frame(matrix(evap0,nrow=724,ncol = 20))

c1 <- c(1:6)
m <- matrix(c1,nrow = 6,ncol = 1,byrow = T) # No Meth + 1 (for legend)
name_plt = paste0(Dir,'PQ_figures/fig4/fig4.png')
png(name_plt,width=750,height=750)
layout(mat = m,heights  = c(0.3,1,0.3,1,0.3,1))#, widths=c(2,2,2,0.5), heights=array(1,dim=c(4)))

nmeth <- c('RAW','LS','QM','BMA')
method <- rep(nmeth,5)
c0 <- c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24)
col0 <- rep(c('aquamarine','coral1','darkgoldenrod1','darkorchid1'),5)

## PLOT PRECIPITATION
# Title
par(mar=c(0,3,1,2))
plot.new()
text(0.5,0.5,"Precipitation",cex=2,font=2)

# Boxplot
par(mar=c(4,5.1,1,2)+0.1)
#boxplot(pcp1,las = 2, at = c0,par(mar=c(12,5,4,2)+0.1),names = method,col = col0,xaxt='n')
boxplot(pcp1,las = 2, at = c0,names = method,col = col0,xaxt='n',ylim=c(-0.1,1),
        cex.axis = 1.5)
at0 <- c(2.5,7.5,12.5,17.5,22.5)
labels0 <- c('pcp = 0','0 < pcp <= 5','5 < pcp <= 10','10 < pcp <= 50','pcp > 50')
axis(1,at=at0,labels = labels0,lwd = 0,cex.axis=1.5)
mtext('crpss',2,3.5)

## PLOT TEMPERATURE
# Title
par(mar=c(0,3,2,1))
plot.new()
text(0.5,0.5,"Temperature",cex=2,font=2)

# Boxplot
par(mar=c(4,5.1,1,2)+0.1)
#boxplot(temp1,las = 2, at = c0,par(mar=c(12,5,4,2)+0.1),names = method,col = col0,xaxt='n')
boxplot(temp1,las = 2, at = c0,names = method,col = col0,xaxt='n',ylim=c(-0.1,1),
        cex.axis = 1.5)
at0 <- c(2.5,7.5,12.5,17.5,22.5)
labels0 <- c('temp <= -5','-5 < temp <= 5','5 < temp <= 15','15 < temp <= 25','temp > 25')
axis(1,at=at0,labels = labels0,lwd = 0,cex.axis=1.5)
mtext('crpss',2,3.5)

## PLOT REF. EVAPOTRANSPIRATION
# Title
par(mar=c(0,3,2,1))
plot.new()
text(0.5,0.5,"Ref. Evapotranspiration",cex=2,font=2)

# Boxplot
par(mar=c(4,5.1,1,2)+0.1)
#boxplot(evap1,las = 2, at = c0,par(mar=c(12,5,4,2)+0.1),names = method,col = col0,xaxt='n')
boxplot(evap1,las = 2, at = c0,names = method,col = col0,xaxt='n',ylim=c(-0.1,1),
        cex.axis = 1.5)
at0 <- c(2.5,7.5,12.5,17.5,22.5)
labels0 <- c('evap = 0','0 < evap <= 1','1 < evap <= 3','3 < evap <= 4','evap > 4')
axis(1,at=at0,labels = labels0,lwd = 0,cex.axis=1.5)
mtext('crpss',2,3.5)

graphics.off()

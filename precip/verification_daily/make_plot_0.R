## Makes plot of skill of different procedures to deal with drizzle effect
## Explanations in PQ_figures

Dir <- 'E:/calibMET/precip/drizzle/'

proc <- 'B'

data <- array(NaN,dim=c(12,662))
for (month in 1:12){
  name <- paste0(Dir,'verL_',proc,'_',month,'.txt')
  data[month,] <- matrix(scan(name),ncol = 1,nrow = 662)
}

name_plt = paste(Dir,'verifL_proc',proc,'.png',sep='')
png(name_plt,width=200,height=200)

par(mar=c(3,2,2,2))
mean0 <- mean(data)
labels0 <-  c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')
main0 <- paste0('Procedure ',proc)
boxplot(t(data),main = main0,xaxt="n",col = 635)#,ylim=c(0,1) ,yaxt="n"
abline(h=mean0,col='red')
at0 <- seq(1,12,1);at1 <- seq(0,1,0.1)
axis(1,at=at0,labels = labels0,cex.axis=0.8,las= 3)
#axis(2,at=at1,labels = at1,cex.axis=0.8,las= 3)

graphics.off()
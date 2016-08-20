## makes Q-Q plots obs-sim values
library(seas)
ObsDir = 'E:/calibMET/precip'
year = seq(1990,2014,1)
days = year.length(year,'julian')

grid_filename = paste(ObsDir,'/obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))

dry_p = array(NaN,dim = c(25,662))
for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  ## Observations
  obs_namefile = paste(ObsDir,'/drizzle/observation/obs_station_',grid[ngrid,2],'.txt',sep='')
  obs = t(matrix(scan(obs_namefile),nrow=367,ncol=25));obs[obs<0]=0;
  obs = obs[,2:367]#; obs = obs[,1:30]
  dry_d = rowSums(obs==0,na.rm = T)
  dry_p[,ngrid] = (dry_d/days)*100
}

## Make map of dry days per year

## Coordinates idw
name_f = paste(ObsDir,'/drizzle/obs_ens_coord_obs.txt',sep="")
coor   = t(matrix(scan(name_f),nrow=3,ncol=662))

# Colors
ColorRamp<- colorRampPalette(colors = c("blue", "red"), space="Lab")(7)
#ColorRamp <- terrain.colors(31, alpha = 1)
#min = min(dry_p,na.rm = TRUE);max=max(dry_p,na.rm = TRUE)
min = 20; max = 80
ColorLevels <- seq(min, max, length=length(ColorRamp))

num <- c(seq(1:25),rep(26,5),rep(27,4),28); heights0 = c(1,1,1,1,1,0.3,1.5) 
m <- matrix(num,nrow = 7,ncol = 5,byrow = TRUE) # No Meth + 1 (for legend)
name_plt = paste0(ObsDir,'/drizzle/dry_days_year.png')
png(name_plt,width=750,height=750)
#par(cex.lab=3)
layout(mat = m, heights= heights0)
par(mar=c(1,1,1,1))

for (ny in 1:length(year)) { # Start lead time
  grd1 <- data.frame(z = dry_p[ny,], xc = coor[,2], yc = coor[,3])
  coordinates(grd1) <- ~ xc + yc
  gridded(grd1) <- TRUE
  grd1 <- as(grd1, "SpatialGridDataFrame")
  #par(mar = c(0,0,2,0))
  main0 = year[ny] 
  image(grd1, col = ColorRamp, main = main0,zlim=c(min,max))
}
  #Plot legend
par(mar = c(2,1,1,1))
image(ColorLevels, 1,
        matrix(data=ColorLevels, ncol=1,nrow=length(ColorLevels)),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n",axes=F)

axis(1,at=ColorLevels, labels=ColorLevels, las= HORIZONTAL<-1,
              cex.axis=1)

#Plot all
par(mar = c(2,2.5,1,1))
boxplot(t(dry_p),xaxt = 'n',axes = F,ylim = c(0,100))
axis(1,at=c(1:25), labels=year, las= HORIZONTAL<-1,
     cex.axis=1)
axis(2,at=seq(0,100,10), labels=seq(0,100,10), las= HORIZONTAL<-1,
     cex.axis=1)

#Plot all
par(mar = c(2,0,1,1))
dry_a = matrix(dry_p,nrow = 662*25,ncol = 1)
boxplot(dry_a,xaxt = 'n',axes = F,ylim = c(0,100))
axis(1,at=1, labels='all', las= HORIZONTAL<-1,
     cex.axis=1)
#axis(2,at=seq(0,100,10), labels=seq(0,100,10), las= HORIZONTAL<-1,
#     cex.axis=1)
  
graphics.off()

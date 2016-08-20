## makes Q-Q plots obs-sim values

ObsDir = 'E:/calibMET/precip'
month = seq(1,12,1);mlead = seq(1,7,1)
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)
Nyears <- c(25,24,25,25,25,24,24,23,24,24,24,24)

library(lubridate)#;library(ncdf4)
library(RNetCDF)
source(paste0(ObsDir,'/drizzle/dry_thre.R'))

grid_filename = paste(ObsDir,'/obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
info_filename = paste(ObsDir,'/ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))

## Read percentages
dry <- array(NaN,dim = c(662,25,12))
for (im in 1:12){
  name <- paste0(ObsDir,'/drizzle/dry_per_m_',im,'.txt')
  dry[,,im] <- t(matrix(scan(name),nrow=25,ncol=662))
}

nmonth <- 12
index     = which(ens_info[1:293,3] == month[nmonth])
#index_out = which(ens_info[index,4] != Nens[nmonth]) ## For Feb and Aug
#index     = index[-(index_out)]
#dry <- dry[,-(index_out),]
ens_index = ens_info[index,c(31,32)]
no_ens    = Nens[nmonth];dimIndex  = length(index)
index_f   = matrix(nrow=1, ncol=dimIndex*no_ens)
index_f   = sapply(1:dimIndex, function(ix) c(ens_index[ix,1]:ens_index[ix,2]))
nyears <- dimIndex

thre <- array(NaN,dim = c(no_ens,dimIndex,7,dim(grid)[1]))
for (ngrid in 1:dim(grid)[1]){ # Loop for space :dim(grid)[1]
  print(ngrid)
  ## Ensembles Raw
  ens_filename0 <- paste(ObsDir,'/ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
  ens  <- t(matrix(scan(ens_filename0),nrow=217,ncol=8322))
  ens  <- ens[index_f,3:217]
  thre[,,,ngrid] <- dry_thre(ens,nmonth,dry,ngrid,no_ens,nyears)
  rm(ens)
} # End grid

#par(mar=c(1,1,1,1))
#m <- matrix(c(1:7),nrow = 7,ncol = 1,byrow = TRUE) # No Meth + 1 (for legend)
#name_plt = paste0(ObsDir,'/drizzle/thre_m',nmonth,'.png')
#png(name_plt,width=750,height=750)
#layout(mat = m)
# at0 <- seq(0,10,2)
# for (im in 1:7){
#   par(mar=c(2,2.5,1,1),cex.axis=c(1.5))
#   thre1 <- matrix(thre[,,im,],ncol = no_ens,nrow = nyears*662,byrow = T)
#   boxplot(thre1,ylim=c(0,10),xlab="",ylab="",xaxt="n",axes = FALSE)
#   abline(h=1.5,col='red')
#   axis(2,at=at0,labels = at0,las= HORIZONTAL<-1)
#   box()
# }
# axis(1,at=c(1:no_ens),labels = c(1:no_ens),las= HORIZONTAL<-1)
#graphics.off()

# 
# dimX = ncdim_def("no_ens", "no", c(1:no_ens))
# dimY = ncdim_def("no_year", "no", c(1:nyears))
# dimZ = ncdim_def("no_lead", "no", c(1:7))
# dimW = ncdim_def("no_grid", "no", c(1:662))

thre[which(is.na(thre))]=-9999
# fillvalue <- -9999
# 
# dlname <- "percentage dry days"
# tmp_def <- ncvar_def("per","per",list(dimX,dimY,dimZ,dimW),fillvalue,dlname,prec="single")
# name_cdf <- paste0(ObsDir,'/drizzle/percentage_m',nmonth,'.nc')
# ncout <- nc_create(name_cdf,list(tmp_def))
# ncvar_put(ncout, tmp_def, thre, start=c(1, 1, 1, 1))
# nc_close(ncout)


name_cdf <- paste0(ObsDir,'/drizzle/percentage_m',nmonth,'.nc')
ncout <- create.nc(name_cdf)

dim.def.nc(ncout,"no_ens", no_ens)
dim.def.nc(ncout,"no_year", nyears)
dim.def.nc(ncout,"no_lead", 7)
dim.def.nc(ncout,"no_grid", 662)

var.def.nc(ncout, "per", "NC_DOUBLE", c("no_ens","no_year","no_lead","no_grid"))

att.put.nc(ncout, "per", "missing_value", "NC_DOUBLE", -9999)

var.put.nc(ncout, "per", thre, c(1, 1, 1, 1), c(no_ens,nyears,7,662))

sync.nc(ncout)

close.nc(ncout)




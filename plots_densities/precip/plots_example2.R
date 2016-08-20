## Makes plots of examples for 2003 02 lt +30 grid 646

ObsDir = 'E:/calibMET/precip/'
grid_filename = paste(ObsDir,'obs_ens.txt',sep='')
grid = t(matrix(scan(grid_filename),nrow=6,ncol=662))
info_filename = paste(ObsDir,'ensem_info.txt',sep='')
ens_info = t(matrix(scan(info_filename),nrow=33,ncol=311))
Nens = c(15,51,15,15,51,15,15,51,15,15,51,15)

ngrid = 269

obs_filename = paste(ObsDir,'observation/obs_station_',grid[ngrid,2],'.txt',sep='')
obs = t(matrix(scan(obs_filename),nrow=216,ncol=293));obs[obs<0]=0;
ens_filename = paste(ObsDir,'ensembles_idw/ens_station_',grid[ngrid,2],'.txt',sep='')
ens = t(matrix(scan(ens_filename),nrow=217,ncol=8322)); ens[ens<0]=0;

nmonth = 1 # Try a month of 15 members for better visualization?

index     = which(ens_info[1:293,3] == nmonth)
#index_out = which(ens_info[index,4] != Nens[nmonth])
#index     = index[-(index_out)]
ens_index = ens_info[index,c(31,32)];no_ens    = Nens[nmonth]
dimIndex  = length(index)
index_f   = matrix(nrow=1, ncol=dimIndex*no_ens)
index_f   = sapply(1:dimIndex, function(ix) c(ens_index[ix,1]:ens_index[ix,2]))
# Extract ensembles with start date: nmonth [Jan(1) to Dec(12)]
Pens      = ens[index_f,3:217];rm(ens)
# Extract observations with start date: nmonth [Jan(1) to Dec(12)]
Pobs = obs[index,2:216]


nyear0 = 9 # 2003 01 lt + 30 grid 646
nyear = nyear0 - 1 
nlead = 4
ens0 = Pens[-(index_nyear),nlead] # train ens
obs0 = Pobs[-(nyear+1),nlead] # train obs
dat <- c(ens0,obs0)

a0 = seq(1,no_ens,1)
a1 = rep((no_ens+1),(dimIndex[1]-1))
f0 = rep(iiens,(dimIndex[1]-1))
f1 = c(f0,a1)
data <- data.frame(values = dat, no_ens = f1)

index_nyear = ((nyear*no_ens)+1):((nyear*no_ens)+no_ens)
mens = colMeans(Pens[-(index_nyear),],na.rm = TRUE) # Mean of all = Mean of the mean
mobs = colMeans(Pobs[-(nyear+1),],na.rm = TRUE)
## First method - Bias correction linear scaling
BC_factor = round(mobs[nlead]/mens[nlead],digits = 2) 
ens0 = Pens[-(index_nyear),nlead];d=density(ens0) # train ens
obs0 =Pobs[-(nyear+1),nlead];d1=density(obs0) # train obs
final_ens = t(BC_factor * t(Pens[index_nyear,nlead])) # corr ens
x1 = Pens[index_nyear,nlead] # raw ens
col0 = c('blue','red')
mar0 <- c(3,2,3,1)
mar1 <- c(3,1,3,1)
ylim0 <- c(0,35)

## Thrid method BMA
thres = seq(from=0,to=100,by=1)
Xpredictor  = (Pens[-(index_nyear),nlead]) # train ens
Xpredict    = (Pens[index_nyear,nlead]) # raw ens
Ypredictand = Pobs[-(nyear+1),nlead]  # train obs
train = t(matrix(Xpredictor,nrow = no_ens,ncol=dimIndex[1]-1))
ens_data = ensembleData(forecasts=train,observations=Ypredictand,verbose=FALSE,dates=NULL,forecastHour = NULL,initializationTime = NULL)
calib = fitBMAgamma0(ens_data,exchangeable = exc)
fore = ensembleData(t(Xpredict),verbose=FALSE,dates=NULL,forecastHour = NULL,initializationTime = NULL)
Probcurve = 1 - cdf.fitBMAgamma0(calib, fore, values = thres)

# # START plotting
# m <- matrix(seq(1:6),nrow = 1,ncol = 6,byrow = TRUE) # No Meth + 1 (for legend)
# name_plt = 'E:/calibMET/plots_densities/precip/test_examp.png'
# png(name_plt,width=750,height=200)
# layout(mat = m, widths=c(2,1,2,1,2,1), heights=array(1,dim=c(6)))
# 
# 
# 
# ## Plot LS
# par(mar=mar0)
# main0 = 'LS'
# densityplot(~values, data=data,plot.points = FALSE,from = 0)
# plot(d1,main = main0,lty = 5) # observation train
# plot(d,add = TRUE,col = 'azure3') # ensemble train
# x0 = final_ens # corr ens
# y0 = array(0,dim=c(no_ens,1))
# points(x0,y0,col='red') # final ens
# points(x1,y0,col='blue') # before ens
# abline(v=mens[nlead],col='azure3') # mean ens
# abline(v=mobs[nlead],lty = 5) # mean obs
# abline(v = Pobs[(nyear+1),nlead],col='maroon1') # observation validate
# labels0 = paste('BC factor = ',BC_factor,sep='')
# text(x = 20, y = 0.4, labels = labels0)
# par(mar=mar1)
# ## Boxplot
# data <- data.frame(raw = x1, corr = x0)
# boxplot(data,main = main0,col = col0,ylim = ylim0)
# abline(h = Pobs[(nyear+1),nlead],col='maroon1')
# 
# ## Plot QM
# ens_namefile = paste(ObsDir,'ensembles_c/ens_qm_g',grid[ngrid,2],'_m',nmonth,'.txt',sep='')
# ens = t(matrix(scan(ens_namefile),nrow=215,ncol=dimIndex[1]*no_ens))
# x0 = ens[index_nyear,nlead];rm(ens) # corr ens
# 
# ens00 = matrix(ens0, nrow = no_ens,ncol = (dimIndex[1]-1))  # train ens
# df1 = lapply(1:no_ens, function(iens) density(ens00[iens,])) # ecdf for each ens
# par(mar=mar0)
# main0 = 'QM'
# plot(d1,main = main0,lty = 5) # observation train
# for (iens in 1:no_ens){
#   plot(df1[[iens]],add = TRUE,col = 'azure3')
# }
# points(x1,y0,col='blue') # before ens
# points(x0,y0,col='red') # final ens
# abline(v = Pobs[(nyear+1),nlead],col='maroon1') # observation validate
# ## Boxplot
# par(mar=mar1)
# data <- data.frame(raw = x1, corr = x0)
# boxplot(data,main = main0,col = col0,ylim = ylim0)
# abline(h = Pobs[(nyear+1),nlead],col='maroon1') 
# 
# ## Plot BMA
# ens_namefile = paste(ObsDir,'ensembles_c/ens_bma_g',grid[ngrid,2],'_m',nmonth,'.txt',sep='')
# ens = t(matrix(scan(ens_namefile),nrow=215,ncol=dimIndex[1]*no_ens))
# x0 = ens[index_nyear,nlead];rm(ens) # corr ens
# 
# #ens00 = matrix(ens0, nrow = no_ens,ncol = (dimIndex[1]-1))
# #ecdf1 = lapply(1:no_ens, function(iens) ecdf(ens00[iens,]))
# main0 = 'BMA'
# par(mar=mar0)
# plot(d1,main = main0,lty = 5) # observation train
# for (iens in 1:no_ens){
#   plot(df1[[iens]],add = TRUE,col = 'azure3')
# }
# points(x1,y0,col='blue') # before ens
# points(x0,y0,col='red') # final ens
# abline(v = Pobs[(nyear+1),nlead],col='maroon1') # observation validate
# ## Boxplot
# par(mar=mar1)
# data <- data.frame(raw = x1, corr = x0)
# boxplot(data,main = main0,col = col0,ylim = ylim0)
# abline(h = Pobs[(nyear+1),nlead],col='maroon1')
# 
# graphics.off()

## Computes Eucledian distances and saves closer node

library(matrixStats)

Dir = 'C:/Users/DianaL/Desktop/NCAR_visit/'

filename = paste(Dir,'interpolation/nodes.txt',sep='')
nodes = t(matrix(scan(filename),nrow = 5,ncol=30))
filename = paste(Dir,'interpolation/nodes_obs.txt',sep='')
nd = t(matrix(scan(filename),nrow = 3,ncol=662))
filename = paste(Dir,'interpolation/ens_coord.txt',sep='')
coorEns = t(matrix(scan(filename),nrow = 3,ncol=72))
filename = paste(Dir,'interpolation/obs_coord.txt',sep='')
coorObs = t(matrix(scan(filename),nrow = 3,ncol=662))

min_dist = array(data=NA, dim = c(662))

for (isq in 1:dim(nodes)[1]){ # Start squares dim(nodes)[1]
  print(isq)
  ind = nodes[isq,2:5] # Nodes
  iobs = which(nd[,3]==isq) # Points inside nodes
  # Coordinates edges
  cedges = sapply(1:4, function(ie) coorEns[ind[ie],2:3])
  cedges = t(cedges)
  # Coordinates inside edges
  cind = sapply(1:length(iobs), function(ie) coorObs[iobs[ie],2:3])
  cind = t(cind)
  ## Computes eucledian distances
  distance = array(NaN,dim=c(length(iobs),length(ind)))
  for (iind in 1:length(iobs)){
    dist = sapply(1:4, function(ie) sqrt(sum((c(cedges[ie,1]-cind[iind,1],cedges[ie,2]-cind[iind,2]))^2)))
    Imdist = ind[which(dist == min(dist))] ## Gives the num of node with min distance
    min_dist[iobs[iind]] = Imdist
  }
} # End squares

nodesF = cbind(nd,min_dist)
name_file = paste(Dir,'interpolation/nodes_closer.txt',sep='')
write.table(nodesF, file=name_file, row.names=FALSE, col.names=FALSE)
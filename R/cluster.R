

.ThreeByThreeOffset <- rbind(c(1,0,0),
                             c(-1,0,0),
                             c(0,1,0),
                             c(0,-1,0),
                             c(0,0,1),
                             c(0,0,-1))


dist_weighted_similarity <- function(centers, coords, decay=-.05) {
  csim <- cor(centers) + 1
  D <- as.matrix(dist(coords))
  exp(decay * D) * csim
}

slic_gradient <- function(vox, bvec, mask, offset=.ThreeByThreeOffset) {
  vs <- sweep(offset, 2, vox, "+")
  vsr <- apply(vs,2,range)
  
  mdim <- dim(mask)
  if (vsr[1,1] < 1 || vsr[2,1] > mdim[1]) {
    NA
  } else if (vsr[1,2] < 1 || vsr[2,2] > mdim[2]) {
    NA
  } else if (vsr[1,3] < 1 || vsr[2,3] > mdim[3]) {
    NA
  } else {
    smat <- series(bvec, vs)
    gx <- 1 - cor(smat[,1], smat[,2])
    gy <- 1 - cor(smat[,3], smat[,4])
    gz <- 1 - cor(smat[,5], smat[,6])
    as.vector(gx + gy + gz)
  }
}

min_gradient <- function(bvec, vox, surround, mask) {
  voxneigb <- sweep(surround, 2, vox, "+")
  g <- apply(voxneigb, 1, slic_gradient, bvec, mask)
  voxneigb[which.min(g),]
}

slic_iterate <- function(bvec, voxels, decay, featureCenters, nn.index, nn.dist) {

  unlist(parallel::mclapply(1:nrow(voxels), function(i) {
    print(i)
    centers <- featureCenters[, nn.index[i,]]  
    vals <- series(bvec, voxels[i,,drop=FALSE])
    featsim <- cor(vals, centers)
    scores <- exp(decay * nn.dist[i,]) * (featsim+1)^2
    (nn.index[i,])[which.max(scores)]
  }))
  
}

avgcor <- function(bvec, mask.idx, clusters) {
  mean(sapply(sort(unique(clusters)), function(i) {
    m1 <- series(bvec, mask.idx[which(clusters==i)])
    centroid <- rowMeans(m1)
    mean(cor(m1, centroid))
  }))
}

slic_cluster <- function(mask, bvec, K=500, decay=-.05, iterations=10, nn=8) {
  mask.idx <- which(mask > 0)
  
  ## real coordinates
  cgrid <- indexToCoord(mask, mask.idx)
  ## voxel coordinates
  vgrid <- indexToGrid(mask, mask.idx)
  
  kmeansResult <- kmeans(cgrid, centers=K, iter.max=100)
  spatialCenters <- kmeansResult$centers
  spatialVoxCenters <- round(coordToGrid(mask, spatialCenters))
  
  surround <- as.matrix(expand.grid(x=c(-1,0,1), y=c(-1,0,1), z=c(-1,0,1)))
  
  currentVoxCenters <- t(apply(spatialVoxCenters, 1, function(vox) min_gradient(bvec, vox, surround, mask)))
  currentVoxCoords <- neuroim::gridToCoord(mask, currentVoxCenters)
  featureCenters <- series(bvec, currentVoxCenters)
  
  clusterAssignments <- numeric(length(mask.idx))
  knnres <- get.knnx(currentVoxCoords, cgrid, nn)
  
  print(avgcor(bvec, mask.idx, kmeansResult$cluster))
  
 
  for (i in 1:10) {     
    print(i)
    clusterAssignments <- slic_iterate(bvec, vgrid, -.05, featureCenters, knnres$nn.index, knnres$nn.dist)
    print(avgcor(bvec, mask.idx, clusterAssignments))
    
    voxsplit <- split(data.frame(vgrid), clusterAssignments)
    featureCenters <- do.call(cbind, mclapply(voxsplit, function(vox) rowMeans(series(bvec, as.matrix(vox)))))
    
    currentVoxCoords <- do.call(rbind, mclapply(voxsplit, function(v) {     
      colMeans(neuroim::gridToCoord(mask, v))
    }))
    
    knnres <- get.knnx(currentVoxCoords, cgrid, nn)
  }
  
  bv <- BrainVolume(clusterAssignments, space(mask), indices=mask.idx)
  list(clusvol=bv, clusters=clusterAssignments, centers = featureCenters, coordCenters=currentVoxCoords, voxelSets=voxsplit)
  
}

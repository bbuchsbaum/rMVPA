

.ThreeByThreeOffset <- rbind(c(1,0,0),
                             c(-1,0,0),
                             c(0,1,0),
                             c(0,-1,0),
                             c(0,0,1),
                             c(0,0,-1))



.TwoByTwoOffset <- rbind(c(1,0,0),
                         c(-1,0,0),
                         c(0,1,0),
                         c(0,-1,0))
                             

dist_weighted_similarity <- function(centers, coords, decay=.05) {
  csim <- cor(centers) + 1
  D <- as.matrix(dist(coords))
  exp(-decay * D) * csim
}


get_surround <- function(vox, surround, vol) {
  vs <- sweep(surround, 2, vox, "+")
  vol[vs]
}

slic_gradient <- function(vox, bvec, mask, offset=.TwoByTwoOffset) {
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
    s1 <- series(bvec, vox[1], vox[2], vox[3])
    sum(1 - cor(s1,smat))
    #sum(sapply(seq(1, ncol(smat), by=2), function(k) {
    #  1 - cor(smat[,k], smat[,k+1])
    #}))
  }
}

min_gradient <- function(bvec, vox, surround, mask) {
  voxneigb <- sweep(surround, 2, vox, "+")
  g <- apply(voxneigb, 1, slic_gradient, bvec, mask)
  if (all(is.na(g))) {
    vox 
  } else {
    as.vector(voxneigb[which.min(g),])
  }
}


slic_iterate <- function(bvec, voxels, decay, featureCenters, nn.index, nn.dist) {
  message("decay", decay)

  
  unlist(parallel::mclapply(1:nrow(voxels), function(i) {
    if (i %% 1000 == 0) {
      message("iteration at: ", round(i/nrow(voxels) * 100), "%")
    }
    
    centers <- featureCenters[, nn.index[i,]]  
    vals <- series(bvec, voxels[i,,drop=FALSE])
    featsim <- cor(vals, centers)
    scores <-  (exp(-decay * nn.dist[i,])) + featsim
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

shrink_vector <- function(mask, vgrid, bvec, k=5, iter=1, radius=NULL) {
  if (is.null(radius)) {
    radius <- floor(mean(spacing(mask)) * 2)
  }

  ovec <- bvec
  
  for (i in 1:iter) {
    omat <- do.call(cbind, mclapply(1:nrow(vgrid), function(j) {
      if (j %% 1000 == 0) {
        message("shrinkage at: ", round(j/nrow(vgrid) * 100), "%")
      }
      
      vox <- vgrid[j,]
      cvox <- RegionSphere(mask, vox, radius=radius, nonzero=TRUE)@coords
      vals <- series(ovec, vox[1], vox[2], vox[3])
      vmat <- series(ovec, cvox)
      cres <- cor(vals, vmat)
      rowMeans(series(ovec, cvox[order(cres, decreasing=TRUE)[1:k],]))
    }))
    
    ovec <- SparseBrainVector(omat, space(mask), mask)
  }
  
  ovec
  
}

#' @export
#' @import FNN
slic_cluster <- function(mask, bvec, K=500, decay=.05, iterations=10, nn=8, shrink=0) {
  mask.idx <- which(mask > 0)
  
  ## real coordinates
  cgrid <- indexToCoord(mask, mask.idx)
  ## voxel coordinates
  vgrid <- indexToGrid(mask, mask.idx)
  
  kmeansResult <- kmeans(cgrid, centers=K, iter.max=100)
  spatialCenters <- kmeansResult$centers
  spatialVoxCenters <- round(coordToGrid(mask, spatialCenters))
  
  #surround <- as.matrix(expand.grid(x=c(-1,0,1), y=c(-1,0,1), z=c(-1,0,1)))
  surround <- as.matrix(cbind(as.matrix(expand.grid(x=c(-1,0,1), y=c(-1,0,1))), z=rep(0, 9)))
  
  if (shrink > 0) {
    message("running ", shrink, "shrinkage iterations with k = ", 5)
    bvec <- shrink_vector(mask, vgrid, bvec, k=5, iter=shrink)
  }
  
  currentVoxCenters <- t(apply(spatialVoxCenters, 1, function(vox) min_gradient(bvec, vox, surround, mask)))
  currentVoxCoords <- neuroim:::gridToCoord(mask, currentVoxCenters)
  
 
  
  featureCenters <- series(bvec, currentVoxCenters)
  
  clusterAssignments <- numeric(length(mask.idx))
  knnres <- get.knnx(currentVoxCoords, cgrid, nn)
  
  message("cluster average cor: ", avgcor(bvec, mask.idx, kmeansResult$cluster))
  
 
  for (i in 1:iterations) {     
    print(i)
    clusterAssignments <- slic_iterate(bvec, vgrid, decay, featureCenters, knnres$nn.index, knnres$nn.dist)
    message("cluster average cor: ", avgcor(bvec, mask.idx, clusterAssignments))
    
    voxsplit <- split(data.frame(vgrid), clusterAssignments)
    featureCenters <- do.call(cbind, parallel::mclapply(voxsplit, function(vox) rowMeans(series(bvec, as.matrix(vox)))))
    
    currentVoxCoords <- do.call(rbind, parallel::mclapply(voxsplit, function(v) {     
      colMeans(neuroim:::gridToCoord(mask, as.matrix(v)))
    }))
    
    knnres <- get.knnx(currentVoxCoords, cgrid, nn)
  }
  
  bv <- BrainVolume(clusterAssignments, space(mask), indices=mask.idx)
  list(clusvol=bv, clusters=clusterAssignments, centers = featureCenters, coordCenters=currentVoxCoords, voxelSets=voxsplit)
  
}

#' @export
#' @import FNN
#' 
#' 
# turbo_cluster <- function(mask, bvec, K=500, decay=-.05, iterations=10, nn=8, shrink=0) {
#
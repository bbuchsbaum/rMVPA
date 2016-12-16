
create_sparse_similarity <- function(mask, bvec, dist_thresh=max(spacing(mask)), simfun=cov) {
  library(FNN)
  library(Matrix)
  library(igraph)
  
  mask.idx <- which(mask > 0)
  G <-  indexToCoord(mask, mask.idx)
  
  G_knn <- FNN::get.knn(G,200)
  
  res <- mclapply(1:nrow(G), function(i) {
    print(i)
    nidx  <-  G_knn$nn.index[i,]
    ndist <- G_knn$nn.dist[i,]
    
    #ovals <- mat[,nidx]
    #vals <- mat[,i]
    
    keep <- ndist < dist_thresh
    
    vals <- series(bvec, mask.idx[i] )
    other <- series(bvec,mask.idx[nidx[keep] ])
    cv <- as.vector(simfun(vals, other))
    cbind(i, nidx[keep], cv)
    #sim <- as.vector(cor(vals,ovals))
    #sim[sim < 0] <- 0
    #wsim <- sim * exp(decay * ndist)
    #cbind(i, nidx, wsim)
    
  })
  
  wmat <- do.call(rbind, res)
  wmat <- Matrix::sparseMatrix(dims=c(nrow(G), nrow(G)), i=wmat[,1], j=wmat[,2], x=wmat[,3])
  
  gg <- igraph::graph.adjacency(wmat, mode="undirected", weighted=TRUE)
  
}
create_sparse_adjacency <- function(mask, dist_thresh=max(spacing(mask)), weightfun="binary") {
  library(FNN)
  library(Matrix)
  library(igraph)
  
  mask.idx <- which(mask > 0)
  G <-  indexToCoord(mask, mask.idx)
  G_knn <- FNN::get.knn(G,27)
  
  res <- mclapply(1:nrow(G), function(i) {
    print(i)
    nidx  <-  G_knn$nn.index[i,]
    ndist <- G_knn$nn.dist[i,]
    
    #ovals <- mat[,nidx]
    #vals <- mat[,i]
    
    keep <- ndist < dist_thresh
    cbind(i, nidx[ndist < dist_thresh], rep(1, sum(keep)))
    #sim <- as.vector(cor(vals,ovals))
    #sim[sim < 0] <- 0
    #wsim <- sim * exp(decay * ndist)
    #cbind(i, nidx, wsim)
    
  })
  
  wmat <- do.call(rbind, res)
  wmat <- Matrix::sparseMatrix(dims=c(nrow(G), nrow(G)), i=wmat[,1], j=wmat[,2], x=wmat[,3])
  
  wmat <- (wmat + t(wmat))/2
  gg <- igraph::graph.adjacency(wmat, mode="undirected", weighted=TRUE)
  
}
  
  
  
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
  list(clusvol=bv, 
       clusters=clusterAssignments, 
       centers = featureCenters, 
       coordCenters=currentVoxCoords, 
       voxelSets=voxsplit)
  
}

computeCentroids <- function(bvec, mask.idx, grid, clusters, assignment, medoid=FALSE) {
  if (!medoid) {
    lapply(clusters, function(id) {
      idx <- which(assignment == id)
      mat <- series(bvec, mask.idx[idx])
      coords <- grid[idx,,drop=FALSE]
      list(center=rowMeans(mat), centroid=colMeans(coords))
    })
  } else {
    lapply(clusters, function(id) {
      idx <- which(assignment == id)
      mat <- series(bvec, mask.idx[idx])
      coords <- grid[idx,,drop=FALSE]
      coords_dist <- as.matrix(dist(coords))
      coords_medoid_ind <- which.min(rowSums(coords_dist))
      Dmat <- 1-cor(mat)
      mat_medoid_ind <- which.min(rowSums(Dmat))
      list(center=mat[,mat_medoid_ind], centroid=coords[coords_medoid_ind,])
    })
    
  }
}

## try multiple kmeans initializations, choose one with best intra-cluster correlation.


#' @export
#' @import FNN
#' @import assertthat
#' 
turbo_cluster <- function(mask, bvec, K=500, lambda=.5, iterations=25, connectivity=27, shrink=0, use_medoid=FALSE) {
  assert_that(lambda >= 0 && lambda <= 1)
  assert_that(connectivity > 1 & connectivity <= 27)
  
  mask.idx <- which(mask > 0)
  grid <- indexToCoord(mask, mask.idx)
  vgrid <- indexToGrid(mask, mask.idx)
  
  if (shrink > 0) {
    message("running ", shrink, "shrinkage iterations with k = ", 5)
    bvec <- shrink_vector(mask, vgrid, bvec, k=5, iter=shrink)
  }
  
  kres <- kmeans(grid, K, iter.max=500)
  kvol <- BrainVolume(kres$cluster, space(mask), indices=mask.idx)
  
  clusid <- sort(unique(kres$cluster))
  neib <- get.knn(grid, k=connectivity)
  dthresh <- min(neib$nn.dist[,connectivity])
  
  iter <- 1
  switches <- 1
  iter.max <- iterations
  
  centroids <- computeCentroids(bvec, mask.idx, grid, sort(unique(kres$cluster)), kres$cluster, medoid=use_medoid)
  sp_centroids <- do.call(rbind, lapply(centroids, "[[", "centroid"))
  
  denom <- max(get.knn(sp_centroids, k=1)$nn.dist[,1])
  valmat <- series(bvec, mask.idx)
  
  curclus <- kvol[mask.idx]
  
  
  ## port iteration to rcpp
  while (iter < iter.max && switches > 0) {
    message("turboclust, iteration: ", iter)
    
    newclus <- sapply(1:nrow(vgrid), function(i) {
      
      ind <- neib$nn.index[i,]
      D <- neib$nn.dist[i,]
      keep <- which(D < dthresh)
      ind <- ind[keep]
      
      oclus <- unique(kvol[mask.idx[ind]])
      diffclus <- which(oclus != curclus[i])
      
      if (length(diffclus) > 0) {
      
        candclus <- c(curclus[i], oclus[diffclus])
        
        cval <- sapply(candclus, function(cc) {
           (1 - cor(valmat[,i], centroids[[cc]]$center))/2
        })
        
        dvals <- sapply(candclus, function(cc) {
          sqrt(sum((grid[i,] - centroids[[cc]]$centroid)^2))
        })/denom
        
        cost <- lambda*cval + (1-lambda)*dvals
        #cost <- cvals
        newc <- candclus[which.min(cost)]
      } else {
        curclus[i]
      }
    })
    
    
    
    centroids <- computeCentroids(bvec, mask.idx, grid, sort(unique(newclus)), newclus, medoid=use_medoid)
    switches <- sum(newclus != curclus)
    
    curclus <- newclus
    message("tuboclust: nswitches ", switches)
    
    iter <- iter + 1
  }
  
  kvol <- BrainVolume(newclus, space(mask), indices=mask.idx)
  voxsplit <- split(data.frame(vgrid), newclus)
  
  list(clusvol=kvol, 
       clusters=newclus, 
       centers = do.call(rbind, lapply(centroids, "[[", "center")), 
       coordCenters=do.call(rbind, lapply(centroids, "[[", "centroid")), 
       voxelSets=voxsplit)
  
  
}

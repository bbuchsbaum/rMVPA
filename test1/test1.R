library(neuroim)

PATH <- 
bvec <- loadVector("s3_rbetas_all.nii.gz")
mask <- loadVolume("global_mask.nii")
des <- read.table("mvpa_design.txt", header=TRUE)

dset <- MVPADataset$new(trainVec=bvec,
                        Y=des$syllable,
                        mask=mask,
                        blockVar=des$block)

mask.idx <- which(mask>0)
igrid <- indexToGrid(mask, mask.idx)
voxgrid <- indexToCoord(mask, mask.idx)

kres <- kmeans(voxgrid, 1000, iter.max=100)
clusmap <- split(as.data.frame(igrid), kres$cluster)
hres <- hclust(dist(kres$centers), members=table(kres$cluster))
cset <- lapply(c(8,16,32,64,128,256,500,1000),function(i) cutree(hres, i))

coordsets <- lapply(cset, function(set) {
  s1 <- split(1:length(set), set)
  lapply(s1, function(ind) {
    do.call(rbind, clusmap[ind])
  })
})

model <- loadModel("corclass", config=list(tuneGrid=data.frame(method="spearman", robust=TRUE)))
crossVal <- BlockedCrossValidation(des$block)

fres <- lapply(coordsets, function(clist) {
  print("yes")
  iret <- lapply(clist, function(vox) {
    res <- model$run(dset, as.matrix(vox), crossVal)
    list(vox=vox, perf=performance(res))
  })
  
  out <- array(0, dim(mask))
  for (i in 1:length(iret)) {
    out[as.matrix(iret[[i]]$vox)] <- iret[[i]]$perf[3]
  }
  
  BrainVolume(out, space(mask))
})

for (i in 1:length(fres)) {
  writeVolume(fres[[i]], paste0("hclus_robcor_", i, ".nii"))
}


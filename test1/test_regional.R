library(neuroim)

PATH <- "/Users/bbuchsbaum/analysis/mvpa/NC/"

bvec <- loadVector(paste0(PATH, "/all_betas_enc.nii.gz"))

mask <- loadVolume(paste0(PATH, "/epi_aparc_mask.nii"))
des <- read.table(paste0(PATH, "design_enc.txt"), header=TRUE)

dset <- mvpa_dataset(bvec, mask=mask)
mdes <- mvpa_design(des, des$Video, block_var="run")

cval <- blocked_cross_validation(mdes$block_var)

mod <- mvpa_model(load_model("sda_notune"), dset, mdes, model_type="classification",crossval=cval)

res <- run_regional(mod, mask)
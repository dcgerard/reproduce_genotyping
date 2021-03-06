## Fit updog on all 1000 SNPs

library(updog)
suppressMessages(library(tidyverse))
counts_mat <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)
ploidy <- 6

one_rep <- function(index, counts_mat, size_mat, ploidy) {
  ocounts <- counts_mat[-1, index]
  pcounts <- counts_mat[1, index]
  osize   <- size_mat[-1, index]
  psize   <- size_mat[1, index]
  ocounts <- ocounts[!is.na(ocounts)]
  osize   <- osize[!is.na(osize)]
  utime   <- system.time({
    uout <- updog::flexdog(refvec   = ocounts, 
                           sizevec  = osize,
                           ploidy   = ploidy,
                           p1ref    = pcounts, 
                           p1size   = psize, 
                           model    = "s1",
                           outliers = TRUE,
                           verbose  = FALSE)
  }
  )
  uout$time <- utime
  saveRDS(object = uout, 
          file   = paste0("./Output/updog_fits/uout", index, ".RDS"))
}

library(parallel)
cl <- makeCluster(detectCores() - 2)
parallel::parSapply(cl         = cl, 
                    X          = 1:ncol(counts_mat), 
                    FUN        = one_rep,
                    counts_mat = counts_mat, 
                    size_mat   = size_mat,
                    ploidy     = ploidy)
stopCluster(cl)

# for (index in 1:ncol(counts_mat)) {
#   one_rep(index, counts_mat, size_mat, ploidy)
# }

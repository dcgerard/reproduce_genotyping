## Apply fitPoly to Shirasawa SNPs

library(fitPoly)
suppressMessages(library(updog))
suppressMessages(library(tidyverse))
counts_mat <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)
ploidy <- 6

one_rep <- function(index, counts_mat, size_mat, ploidy) {
  refvec  <- counts_mat[, index]
  sizevec <- size_mat[, index]

  ## Run fitPoly on default settings.
  ## I augment the parent to get fitPoly to sort of fit an S1 population
  fitpoly_df <- data.frame(MarkerName = "SNP",
                           SampleName = 1:(length(refvec) + 1),
                           ratio = c(refvec[1] / sizevec[1], refvec / sizevec))

  pop.parents <- data.frame(popID   = c("pop", "par1", "par2"),
                            parent1 = c("par1", NA, NA),
                            parent2 = c("par2", NA, NA))
  population  <- data.frame(SampleName = 1:(length(refvec) + 1),
                            popID      = "pop")
  population$popID[1] <- "par1"
  population$popID[2] <- "par2"

  ftime <- system.time({
    fp_out <- fitPoly::fitOneMarker(ploidy         = ploidy,
                                    marker         = "SNP",
                                    data           = fitpoly_df,
                                    pop.parents    = pop.parents,
                                    population     = population)
  })
  if (any(!is.na(fp_out$scores))) {
    fp_out$scores    <- fp_out$scores[-1, ]
    fp_out$time      <- ftime
    fp_out$prop_miss <- 1 - mean(fp_out$scores$maxP, na.rm = TRUE)
  }

  saveRDS(object = fp_out,
          file   = paste0("./Output/fp_fits/fpout", index, ".RDS"))
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


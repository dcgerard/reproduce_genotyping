### Simulation Study
set.seed(25346)

## Simulation Function ---------------------------------------------------------
one_rep <- function(unew, usame) {

  set.seed(unew$seed)
  usim             <- usame
  usim$bias        <- unew$bias
  usim$allele_freq <- unew$allele_freq
  usim$od          <- unew$od
  usim$sizevec     <- unew$sizevec

  ## Simulate New Data ----------------------------------------------
  true_geno  <- updog::rgeno(n           = length(usim$sizevec), 
                             ploidy      = usim$ploidy, 
                             model       = "hw",
                             allele_freq = usim$allele_freq)
  refvec     <- updog::rflexdog(sizevec = usim$sizevec, 
                                geno    = true_geno, 
                                ploidy  = usim$ploidy, 
                                seq     = usim$seq_error, 
                                bias    = usim$bias, 
                                od      = usim$od)
  sizevec    <- usim$sizevec

  ## Run updog -------------------------------------------------
  uout <- updog::flexdog(refvec  = refvec, 
                         sizevec = sizevec, 
                         ploidy  = usim$ploidy, 
                         model   = "hw", 
                         verbose = FALSE)

  ## Run Blischak -----------------------------------------------------------------------
  osize_file <- paste0("./Output/blischak_formatted_data_sims/size", unew$seed, ".txt")
  ocounts_file <- paste0("./Output/blischak_formatted_data_sims/counts", unew$seed, ".txt")
  seq_file <- paste0("./Output/blischak_formatted_data_sims/seq_error", unew$seed, ".txt")
  prefix <- paste0("hwe", unew$seed)
  write.table(matrix(sizevec), osize_file, row.names = FALSE, col.names = FALSE)
  write.table(matrix(refvec), ocounts_file, row.names = FALSE, col.names = FALSE)
  write.table(matrix(usim$seq_error), seq_file, row.names = FALSE, col.names = FALSE)
  command_text <- paste0("ebg hwe -p ", usim$ploidy, " -n ", length(refvec), " -l ", 1,
                         " -r ", ocounts_file,
                         " -t ", osize_file,
                         " -e ", seq_file,
                         " 2> trash.txt",
                         " --prefix ", prefix)
  system(command_text)
  system(paste0("mv ", prefix, "* ./Output/blischak_formatted_data_sims/"))
  bgeno <- read.table(paste0("./Output/blischak_formatted_data_sims/", prefix, "-genos.txt"))[, 1]
  ballele_freq <- read.table(paste0("./Output/blischak_formatted_data_sims/", prefix, "-freqs.txt"))[1, 1]

  # ## Run SuperMASSA ----------------------------------------------------------------------
  # ## Does not return estimated genotypes (only returns ploidy and variance estimate) so SKIP
  # super_input_file <- paste0("./Output/blischak_formatted_data_sims/super", unew$seed, ".txt")
  # write.table(x = cbind(paste0("Ind", 1:length(ocounts)), ocounts, osize - ocounts), file = super_input_file,
  #             col.names = FALSE, row.names = FALSE)
  # command_supermass <- paste0("python2 ./supermassa/src/SuperMASSA.py --inference hw --file ", super_input_file,
  #                             " --ploidy_range 6:6")
  # system(command_supermass)

  ## Remove files so don't get bloated
  system(paste0("rm ./Output/blischak_formatted_data_sims/", prefix, "-genos.txt"))
  system(paste0("rm ./Output/blischak_formatted_data_sims/", prefix, "-freqs.txt"))
  system(paste0("rm ", osize_file))
  system(paste0("rm ", ocounts_file))
  system(paste0("rm ", seq_file))

  ## Run fitPoly on default settings.
  fitpoly_df <- data.frame(MarkerName = "SNP",
                           SampleName = 1:length(refvec),
                           ratio = refvec / sizevec)
  fp_out <- fitPoly::fitOneMarker(ploidy = usim$ploidy,
                                  marker = "SNP",
                                  data = fitpoly_df)

  ## Summary Quantities
  parout     <- rep(NA, length = 25)
  parout[1]  <- mean(bgeno == true_geno, na.rm = TRUE)
  parout[2]  <- mean(uout$geno == true_geno, na.rm = TRUE)
  parout[3]  <- sum((bgeno - true_geno) ^ 2, na.rm = TRUE) / sum(!is.na(true_geno) & !is.na(bgeno))
  parout[4]  <- sum((uout$geno - true_geno) ^ 2, na.rm = TRUE) / sum(!is.na(true_geno) & !is.na(uout$geno))
  parout[5]  <- mean(abs(uout$postmean - true_geno))
  parout[6]  <- usim$allele_freq
  parout[7]  <- uout$par$alpha
  parout[8]  <- ballele_freq
  parout[9]  <- usim$od
  parout[10] <- uout$od
  parout[11] <- usim$bias
  parout[12] <- uout$bias
  parout[13] <- usim$seq_error
  parout[14] <- uout$seq
  if (any(!is.na(fp_out$scores))) {
    parout[15] <- mean(fp_out$scores$maxgeno == true_geno, na.rm = TRUE)
    parout[16] <- sum((fp_out$scores$maxgeno - true_geno) ^ 2, na.rm = TRUE) / sum(!is.na(true_geno) & !is.na(fp_out$scores$maxgeno))
    parout[17] <- 1 - mean(fp_out$scores$maxP, na.rm = TRUE)
    parout[18] <- stats::cor(x = fp_out$scores$maxgeno, y = true_geno, use = "complete.obs")
    
    fp_out$scores$P1 +
      fp_out$scores$P2 * 2 +
      fp_out$scores$P3 * 3 +
      fp_out$scores$P4 * 4 +
      fp_out$scores$P5 * 5 +
      fp_out$scores$P6 * 6 ->
      fp_pm
    
    parout[19] <- stats::cor(x = fp_pm, y = true_geno, use = "complete.obs")
    
  } else {
    parout[15:19] <- NA
  }
  parout[20] <- uout$prop_mis
  parout[21] <- stats::cor(x = uout$geno, y = true_geno, use = "complete.obs")
  parout[22] <- stats::cor(x = uout$postmean, y = true_geno, use = "complete.obs")
  parout[23] <- stats::cor(x = bgeno, y = true_geno, use = "complete.obs")
  parout[24] <- stats::cor(x = refvec / sizevec, y = true_geno, use = "complete.obs")
  parout[25] <- unew$seed
  names(parout) <- c("bham", "uham", "bmse", "umse",
                     "umean_mse", "allele_freq", "uallele_freq",
                     "ballele_freq", "od_param", "uod_param", "bias_val",
                     "ubias_val", "seq_error", "useq_error",
                     "fpham", "fpmse", "fpepm", "fpcor", "fpcor_pm",
                     "uepm", "ucor", "ucor_pm", "bcor", "naive_cor",
                     "seed")
  return(parout)
}

## Read in size data to get realistic size distribution --------------------------------
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)


## Parameters to explore --------------------------------------------
bias_seq  <- c(1, 0.75, 0.5)
seq_error <- 0.005 ## Constant throughout
out_prop  <- 0
od_seq    <- c(0, 0.01, 0.02)
ploidy    <- 6
itermax   <- 1000

## Set up `updog` object parameters that don't vary ------------------
usame               <- list()
usame$ploidy  <- ploidy
usame$model   <- "hw"
usame$seq_error     <- seq_error


## Run Simulations ------------------------------------------------
parvals <- expand.grid(allele_freq = seq(0.05, 0.95, length = itermax), bias = bias_seq, od = od_seq)
parvals$seed <- 1:nrow(parvals)
parvals$sizevec <- sapply(size_mat[, 1:itermax], function(x) x[!is.na(x)])

set.seed(1) ## reorder to evenly distribution cluster computation
parvals_order <- sample(1:nrow(parvals))
parvals <- parvals[parvals_order, ]

par_list <- list()
for (list_index in 1:nrow(parvals)) {
  par_list[[list_index]] <- list()
  for (inner_list_index in 1:(ncol(parvals) - 1)) {
    par_list[[list_index]][[inner_list_index]] <- parvals[list_index, inner_list_index]
    names(par_list[[list_index]])[inner_list_index] <- colnames(parvals)[inner_list_index]
  }
}

for (list_index in 1:nrow(parvals)) {
  par_list[[list_index]]$sizevec <- parvals$sizevec[[list_index]]
}

cl <- parallel::makeCluster(parallel::detectCores() - 1)
simout <- t(parallel::parSapply(cl = cl, X = par_list, FUN = one_rep, usame = usame))
parallel::stopCluster(cl)

# simout <- matrix(NA, nrow = nrow(parvals), ncol = 14)
# colnames(simout) <- c("bham", "uham", "bmse", "umse", "umean_mse", "allele_freq", "uallele_freq", "ballele_freq",
#                       "od_param", "uod_param", "bias_val", "ubias_val", "seq_error", "useq_error")
# for (index in 1:nrow(parvals)) {
#   cat(index, "\n")
#   simout[index, ] <- one_rep(unew = par_list[[index]], usame = usame)
# }

write.csv(simout, "./Output/sims_out/sims_out.csv", row.names = FALSE)

### pp sims

## Simulation Function ---------------------------------------------------------
one_rep <- function(unew, usame) {

  set.seed(unew$seed)
  usim             <- usame
  usim$bias        <- unew$bias
  usim$weight_vec  <- c(unew$first_weight, 1 - unew$first_weight)
  usim$od          <- unew$od
  usim$sizevec     <- unew$sizevec

  ## Simulate New Data ----------------------------------------------
  true_geno  <- updog::rgeno(n               = length(usim$sizevec),
                             ploidy          = usim$ploidy,
                             model           = "s1pp",
                             p1geno          = 2,
                             p1_pair_weights = usim$weight_vec)
  refvec     <- updog::rflexdog(sizevec = usim$sizevec,
                                geno    = true_geno,
                                ploidy  = usim$ploidy,
                                seq     = usim$seq_error,
                                bias    = usim$bias,
                                od      = usim$od)
  sizevec    <- usim$sizevec

  ## Run updog -------------------------------------------------
  uouts1 <- updog::flexdog(refvec  = refvec,
                           sizevec = sizevec,
                           ploidy  = usim$ploidy,
                           model   = "s1",
                           verbose = FALSE)
  uouts1pp <- updog::flexdog(refvec  = refvec,
                             sizevec = sizevec,
                             ploidy  = usim$ploidy,
                             model   = "s1pp",
                             verbose = FALSE)

  ## Run Li using updog
  liout <- updog::flexdog(refvec      = refvec,
                          sizevec     = sizevec,
                          ploidy      = usim$ploidy,
                          model       = "hw",
                          bias_init   = 1,
                          seq         = usim$seq_error,
                          od          = 0,
                          update_bias = FALSE,
                          update_od   = FALSE,
                          update_seq  = FALSE,
                          verbose     = FALSE)

  ## Run GATK method using updog
  suppressWarnings({
    gatk_out <- updog::flexdog(refvec      = refvec,
                               sizevec     = sizevec,
                               ploidy      = usim$ploidy,
                               model       = "uniform",
                               bias_init   = 1,
                               od          = 0,
                               seq         = usim$seq_error,
                               update_bias = FALSE,
                               update_od   = FALSE,
                               update_seq  = FALSE,
                               verbose     = FALSE)
  })

  ## Run fitPoly on default settings.
  fitpoly_df <- data.frame(MarkerName = "SNP",
                           SampleName = 1:length(refvec),
                           ratio = refvec / sizevec)
  fp_out <- fitPoly::fitOneMarker(ploidy = usim$ploidy,
                                  marker = "SNP",
                                  data = fitpoly_df)

  parout <- rep(NA, length = 28)
  parout[1] <- mean(uouts1$geno != true_geno, na.rm = TRUE)
  parout[2] <- mean(uouts1pp$geno != true_geno, na.rm = TRUE)
  parout[3] <- mean(liout$geno != true_geno, na.rm = TRUE)
  parout[4] <- mean(gatk_out$geno != true_geno, na.rm = TRUE)
  suppressWarnings({
    parout[5] <- stats::cor(uouts1$postmean, true_geno, use = "complete.obs")
    parout[6] <- stats::cor(uouts1pp$postmean, true_geno, use = "complete.obs")
    parout[7] <- stats::cor(liout$postmean, true_geno, use = "complete.obs")
    parout[8] <- stats::cor(gatk_out$postmean, true_geno, use = "complete.obs")
  })
  parout[9] <- uouts1$prop_mis
  parout[10] <- uouts1pp$prop_mis
  parout[11] <- liout$prop_mis
  parout[12] <- gatk_out$prop_mis
  if (any(!is.na(fp_out$scores))) {
    parout[13] <- mean(fp_out$scores$maxgeno != true_geno, na.rm = TRUE)
    parout[14] <- 1 - mean(fp_out$scores$maxP, na.rm = TRUE)

    fp_out$scores$P1 +
      fp_out$scores$P2 * 2 +
      fp_out$scores$P3 * 3 +
      fp_out$scores$P4 * 4 ->
      fp_pm

    suppressWarnings({
      parout[15] <- stats::cor(x = fp_pm, y = true_geno, use = "complete.obs")
    })

  } else {
    parout[13:15] <- NA
  }
  parout[16] <- uouts1$par$pgeno
  parout[17] <- uouts1pp$par$p1geno
  parout[18] <- uouts1pp$par$p1_pair_weights[1]
  parout[19] <- uouts1$bias
  parout[20] <- uouts1pp$bias
  parout[21] <- usim$bias
  parout[22] <- uouts1$seq
  parout[23] <- uouts1pp$seq
  parout[24] <- usim$seq_error
  parout[25] <- uouts1$od
  parout[26] <- uouts1pp$od
  parout[27] <- usim$od
  parout[28] <- unew$first_weight

  names(parout) <- c("s1_ham",
                     "s1pp_ham",
                     "hw_ham",
                     "u_ham",
                     "s1_cor",
                     "s1pp_cor",
                     "hw_cor",
                     "u_cor",
                     "s1_epm",
                     "s1pp_epm",
                     "hw_epm",
                     "u_epm",
                     "fp_ham",
                     "fp_epm",
                     "fp_cor",
                     "s1_pgeno",
                     "s1pp_pgeno",
                     "s1pp_firstweight",
                     "s1_bias",
                     "s1pp_bias",
                     "bias",
                     "s1_seq",
                     "s1pp_seq",
                     "seq",
                     "s1_od",
                     "s1pp_od",
                     "od",
                     "firstweight")

  return(parout)
}

## Read in size data to get realistic size distribution --------------------------------
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)


## Parameters to explore --------------------------------------------
bias_seq  <- c(1, 0.75, 0.5)
seq_error <- 0.005 ## Constant throughout
out_prop  <- 0
od_seq    <- c(0, 0.01, 0.02)
ploidy    <- 4
itermax   <- 500

## Set up `updog` object parameters that don't vary ------------------
usame           <- list()
usame$ploidy    <- ploidy
usame$model     <- "hw"
usame$seq_error <- seq_error


## Run Simulations ------------------------------------------------
parvals <- expand.grid(seed         = 1:itermax,
                       first_weight = c(0, 1),
                       bias         = bias_seq,
                       od           = od_seq)
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

write.csv(simout, "./Output/sims_out/sims_out_pp.csv", row.names = FALSE)

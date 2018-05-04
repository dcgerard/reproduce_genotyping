## Simulate under F1 but fit both F1 and HWE
one_rep <- function(usame, unew) {
  set.seed(unew$seed)
  usim         <- usame
  usim$bias    <- unew$bias
  usim$od      <- unew$od
  usim$sizevec <- unew$sizevec
  usim$p1geno  <- unew$pgeno

  ## Simulate new data
  true_geno  <- updog::rgeno(n      = length(usim$sizevec), 
                             ploidy = usim$ploidy,
                             model  = "s1", 
                             p1geno = usim$p1geno)
  refvec <- updog::rflexdog(sizevec = usim$sizevec, 
                            geno    = true_geno, 
                            ploidy  = usim$ploidy, 
                            seq     = usim$seq_error, 
                            bias    = usim$bias,
                            od      = usim$od)
  sizevec <- usim$sizevec

  ## run updog assuming HWE
  uout <- updog::flexdog(refvec  = refvec,
                         sizevec = sizevec, 
                         ploidy  = usim$ploidy, 
                         model   = "hw",
                         verbose = FALSE)

  ## run updog assuming S1
  us1_out  <- updog::flexdog(refvec  = refvec,
                             sizevec = sizevec,
                             ploidy  = usim$ploidy,
                             model   = "s1",
                             verbose = FALSE)

  ## Summaries
  parout <- rep(NA, length = 16)

  parout[1] <- mean(us1_out$geno == true_geno, na.rm = TRUE)
  parout[2] <- mean(uout$geno == true_geno, na.rm = TRUE)

  parout[3] <- mean((us1_out$postmean - true_geno)^2, na.rm = TRUE)
  parout[4] <- mean((uout$postmean - true_geno)^2, na.rm = TRUE)

  parout[5] <- usim$p1geno
  parout[6] <- us1_out$par$pgeno
  parout[7] <- uout$par$alpha

  parout[8] <- usim$od
  parout[9] <- us1_out$od
  parout[10] <- uout$od

  parout[11] <- usim$bias
  parout[12] <- us1_out$bias
  parout[13] <- uout$bias

  parout[14] <- usim$seq_error
  parout[15] <- us1_out$seq
  parout[16] <- uout$seq
  
  parout[17] <- unew$seed

  names(parout) <- c("sham", "hham", "smse", "hmse", 
                     "pgeno", "spgeno", "hallele_freq",
                     "od", "sod", "hod", "bias", "sbias", 
                     "hbias", "seq", "sseq", "hseq", "seed")

  return(parout)
}

library(updog)
ploidy <- 6
pvals <- c(5, 4, 3)
qarray <- get_q_array(ploidy)
bias_seq  <- c(1, 0.75, 0.5)
seq_error <- 0.005 ## Constant throughout
od_seq    <- c(0, 0.01, 0.02)
itermax   <- 1000

## Read in size data to get realistic size distribution --------------------------------
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)

parvals       <- expand.grid(pgeno = pvals, bias = bias_seq, seq = seq_error, od = od_seq, seed = 1:itermax)
parvals$sizevec <- sapply(size_mat[, 1:itermax], function(x) x[!is.na(x)])

## Reorder so that computation isn't heavy loaded on some cores
set.seed(1)
par_order <- sample(1:nrow(parvals))
parvals <- parvals[par_order, ]

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

## These parameters do not change
usame           <- list()
usame$seq_error <- seq_error
usame$ploidy    <- ploidy

## Run method

cl <- parallel::makeCluster(parallel::detectCores() - 2)
simout <- t(parallel::parSapply(cl = cl, X = par_list, FUN = one_rep, usame = usame))
parallel::stopCluster(cl)

write.csv(simout, "./Output/hwe_f1_sims/hwe_f1_sims_out.csv", row.names = FALSE)

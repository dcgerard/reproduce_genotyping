## Simulate under F1 but fit both F1 and HWE
one_rep <- function(usame, unew) {
  set.seed(unew$seed)
  usim             <- usame
  usim$bias_val    <- unew$bias
  usim$od_param    <- unew$od
  usim$input$osize <- unew$input$osize
  usim$p1geno      <- unew$pgeno
  usim$p2geno      <- unew$pgeno

  ## Simulate new data
  rout <- updog::rupdog(usim)
  ocounts    <- rout$input$ocounts
  osize      <- rout$input$osize
  true_ogeno <- rout$ogeno

  ## run updog assuming HWE
  bias_start <- exp(-2:2 * 0.7) ## plus to minus three sd
  llike_old <- -Inf
  for (index in 1:length(bias_start)) {
    utemp <- updog::updog_vanilla(ocounts = ocounts, osize = osize, ploidy = usim$input$ploidy, model = "hw",
                                  out_prop = 0, update_outprop = FALSE, bias_val = bias_start[index], non_mono_max = Inf)
    if (utemp$llike > llike_old) {
      uout <- utemp
      llike_old <- uout$llike
    }
  }

  ## run updog assuming S1
  us1_out  <- updog::updog(ocounts = ocounts, osize = osize, ploidy = usim$input$ploidy, model = "s1", out_prop = 0,
                           update_outprop = FALSE, non_mono_max = Inf)

  ## Summaries
  parout <- rep(NA, length = 16)

  parout[1] <- mean(us1_out$ogeno == rout$ogeno, na.rm = TRUE)
  parout[2] <- mean(uout$ogeno == rout$ogeno, na.rm = TRUE)

  parout[3] <- mean((us1_out$postmean - rout$ogeno)^2, na.rm = TRUE)
  parout[4] <- mean((uout$postmean - rout$ogeno)^2, na.rm = TRUE)

  parout[5] <- rout$p1geno
  parout[6] <- us1_out$p1geno
  parout[7] <- uout$allele_freq

  parout[8] <- rout$od_param
  parout[9] <- us1_out$od_param
  parout[10] <- uout$od_param

  parout[11] <- rout$bias_val
  parout[12] <- us1_out$bias_val
  parout[13] <- uout$bias_val

  parout[14] <- rout$seq_error
  parout[15] <- us1_out$seq_error
  parout[16] <- uout$seq_error

  names(parout) <- c("sham", "hham", "smse", "hmse", "pgeno", "spgeno", "hallele_freq",
                     "od", "sod", "hod", "bias", "sbias", "hbias", "seq", "sseq", "hseq")

  return(parout)
}

library(updog)
ploidy <- 6
pvals <- c(5, 4, 3)
qarray <- get_q_array(ploidy)
bias_seq  <- c(1, 0.75, 0.5, 0.25)
seq_error <- 0.005 ## Constant throughout
out_prop  <- 0
od_seq    <- c(0, 0.01, 0.05)
itermax   <- 1000

## Read in size data to get realistic size distribution --------------------------------
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)

parvals       <- expand.grid(pgeno = pvals, bias = bias_seq, seq = seq_error, od = od_seq, seed = 1:itermax)
parvals$osize <- sapply(size_mat[, 1:itermax], function(x) x[!is.na(x)])

par_list <- list()
for (list_index in 1:nrow(parvals)) {
  par_list[[list_index]] <- list()
  for (inner_list_index in 1:(ncol(parvals) - 1)) {
    par_list[[list_index]][[inner_list_index]] <- parvals[list_index, inner_list_index]
    names(par_list[[list_index]])[inner_list_index] <- colnames(parvals)[inner_list_index]
  }
}

for (list_index in 1:nrow(parvals)) {
  par_list[[list_index]]$input$osize <- parvals$osize[[list_index]]
}

## These parameters do not change
usame <- list()
usame$seq_error    <- seq_error
usame$input$ploidy <- ploidy
usame$out_prop     <- out_prop
class(usame)       <- "updog"
usame$input$model  <- "s1"
usame$out_mean     <- 1/2
usame$out_disp     <- 1/3
usame$allele_freq  <- -1

## Run method

cl <- parallel::makeCluster(parallel::detectCores() - 2)
simout <- t(parallel::parSapply(cl = cl, X = par_list, FUN = one_rep, usame = usame))
parallel::stopCluster(cl)

write.csv(simout, "./Output/hwe_f1_sims/hwe_f1_sims_out.csv", row.names = FALSE)

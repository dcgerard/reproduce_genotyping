### Simulation Study
set.seed(25346)

## Simulation Function ---------------------------------------------------------
one_rep <- function(unew, usame) {

  err_file_name <- paste0("./Output/blischak_formatted_data_sims/err", unew$seed, ".txt")
  zz <- file(err_file_name, open = "wt")
  sink(zz, type = "message")
  try({
  set.seed(unew$seed)
  usim             <- usame
  usim$bias_val    <- unew$bias_val
  usim$allele_freq <- unew$allele_freq
  usim$od_param    <- unew$od_param
  usim$input$osize <- unew$input$osize

  ## Simulate New Data ----------------------------------------------
  rout       <- updog::rupdog(usim)
  ocounts    <- rout$input$ocounts
  osize      <- rout$input$osize
  true_ogeno <- rout$ogeno

  ## Run updog -------------------------------------------------
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

  ## Run Blischak -----------------------------------------------------------------------
  osize_file <- paste0("./Output/blischak_formatted_data_sims/size", unew$seed, ".txt")
  ocounts_file <- paste0("./Output/blischak_formatted_data_sims/counts", unew$seed, ".txt")
  seq_file <- paste0("./Output/blischak_formatted_data_sims/seq_error", unew$seed, ".txt")
  prefix <- paste0("hwe", unew$seed)
  write.table(matrix(osize), osize_file, row.names = FALSE, col.names = FALSE)
  write.table(matrix(ocounts), ocounts_file, row.names = FALSE, col.names = FALSE)
  write.table(matrix(usame$seq_error), seq_file, row.names = FALSE, col.names = FALSE)
  command_text <- paste0("ebg hwe -p ", usim$input$ploidy, " -n ", length(ocounts), " -l ", 1,
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

  ## Summary Quantities
  parout     <- rep(NA, length = 14)
  parout[1]  <- mean(bgeno == rout$ogeno, na.rm = TRUE)
  parout[2]  <- mean(uout$ogeno == rout$ogeno, na.rm = TRUE)
  parout[3]  <- sum((bgeno - rout$ogeno) ^ 2, na.rm = TRUE) / sum(!is.na(rout$ogeno) & !is.na(bgeno))
  parout[4]  <- sum((uout$ogeno - rout$ogeno) ^ 2, na.rm = TRUE) / sum(!is.na(rout$ogeno) & !is.na(uout$ogeno))
  parout[5]  <- mean(abs(uout$postmean - uout$ogeno))
  parout[6]  <- rout$allele_freq
  parout[7]  <- uout$allele_freq
  parout[8]  <- ballele_freq
  parout[9]  <- rout$od_param
  parout[10] <- uout$od_param
  parout[11] <- rout$bias_val
  parout[12] <- uout$bias_val
  parout[13] <- rout$seq_error
  parout[14] <- uout$seq_error
  names(parout) <- c("bham", "uham", "bmse", "umse", "umean_mse", "allele_freq", "uallele_freq", "ballele_freq",
                     "od_param", "uod_param", "bias_val", "ubias_val", "seq_error", "useq_error")

  ## only errors below this are not recorded
  sink(type="message")
  close(zz)
  system(paste0("rm ", err_file_name))

  return(parout)
  })
}

## Read in size data to get realistic size distribution --------------------------------
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)


## Parameters to explore --------------------------------------------
bias_seq  <- c(1, 0.75, 0.5, 0.25)
seq_error <- 0.005 ## Constant throughout
out_prop  <- 0
od_seq    <- c(0, 0.01, 0.05)
ploidy    <- 6
itermax   <- 1000

## Set up `updog` object parameters that don't vary ------------------
usame               <- list()
usame$input$ploidy  <- ploidy
usame$input$model   <- "hw"
usame$seq_error     <- seq_error
usame$out_prop      <- out_prop
usame$p1geno        <- -1
usame$p2geno        <- -1
usame$out_mean      <- 1/2
usame$out_disp      <- 1/3
class(usame)        <- "updog"


## Run Simulations ------------------------------------------------
parvals <- expand.grid(allele_freq = seq(0.05, 0.95, length = itermax), bias_val = bias_seq, od_param = od_seq)
parvals$seed <- 1:nrow(parvals)
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

cl <- parallel::makeCluster(parallel::detectCores() - 2)
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




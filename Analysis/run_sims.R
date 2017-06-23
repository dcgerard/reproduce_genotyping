### Simulation Study
set.seed(25346)

## Simulation Function ---------------------------------------------------------
one_rep <- function(unew, usame) {
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
  uout <- updog::updog_vanilla(ocounts = ocounts, osize = osize, ploidy = usim$input$ploidy, model = "hw",
                               out_prop = 0, update_outprop = FALSE, non_mono_max = Inf)

  ## Run Blischak -----------------------------------------------------------------------
  write.table(matrix(osize), "./Output/blischak_formatted_data_sims/size.txt", row.names = FALSE, col.names = FALSE)
  write.table(matrix(ocounts), "./Output/blischak_formatted_data_sims/counts.txt", row.names = FALSE, col.names = FALSE)
  write.table(matrix(seq_error), "./Output/blischak_formatted_data_sims/seq_error.txt", row.names = FALSE, col.names = FALSE)
  command_text <- paste0("ebg hwe -p ", ploidy, " -n ", length(ocounts), " -l ", 1,
                         " -r ./Output/blischak_formatted_data_sims/counts.txt ",
                         "-t ./Output/blischak_formatted_data_sims/size.txt ",
                         "-e ./Output/blischak_formatted_data_sims/seq_error.txt",
                         " 2> trash.txt")
  system(command_text)
  system("mv hwe* ./Output/blischak_formatted_data_sims/")
  bgeno <- read.table("./Output/blischak_formatted_data_sims/hwe-genos.txt")[, 1]
  ballele_freq <- read.table("./Output/blischak_formatted_data_sims/hwe-freqs.txt")[1, 1]

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
  return(parout)
}

## Read in size data to get realistic size distribution --------------------------------
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)


## Parameters to explore --------------------------------------------
bias_seq  <- c(1, 0.75, 0.5, 0.25)
seq_error <- 0.005 ## Constant throughout
out_prop  <- 0
od_seq    <- c(0, 0.01, 0.1)
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
parvals$osize <- sapply(size_mat, function(x) x[!is.na(x)])

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

simout <- matrix(NA, nrow = nrow(parvals), ncol = 14)
colnames(simout) <- c("bham", "uham", "bmse", "umse", "umean_mse", "allele_freq", "uallele_freq", "ballele_freq",
                      "od_param", "uod_param", "bias_val", "ubias_val", "seq_error", "useq_error")
for (index in 1:nrow(parvals)) {
  cat(index, "\n")
  simout[index, ] <- one_rep(unew = par_list[[index]], usame = usame)
}

write.csv(simout, "./Output/sims_out/sims_out.csv", row.names = FALSE)




### Simulation Study
set.seed(25346)


library(updog)

size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)
counts_mat <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)


bias_seq  <- c(1, 0.75, 0.5, 0.25)
seq_error <- 0.005 ## Constant throughout
out_prop  <- 0
od_seq    <- c(0, 0.01, 0.1)
ploidy    <- 6



## These vary ----
bias_val <- bias_seq[2]
allele_freq <- runif(1)
od_param <- od_seq[3]
## ---------------
index <- 1
usim               <- list()
usim$input$osize   <- size_mat[, index]
usim$input$ploidy  <- ploidy
usim$input$model   <- "hw"
usim$seq_error     <- seq_error
usim$out_prop      <- out_prop
usim$bias_val      <- bias_val
usim$p1geno        <- -1
usim$p2geno        <- -1
usim$allele_freq   <- allele_freq
usim$od_param      <- od_param
usim$out_mean      <- 1/2
usim$out_disp      <- 1/3
class(usim)        <- "updog"

rout <- rupdog(usim)
ocounts <- rout$input$ocounts
osize   <- rout$input$osize
true_ogeno <- rout$ogeno
true_ogeno[is.na(true_ogeno)] <- -1

uout <- updog_vanilla(ocounts = ocounts, osize = osize, ploidy = ploidy, model = "hw", out_prop = 0, update_outprop = FALSE, non_mono_max = Inf)
uout$ogeno
rout$ogeno

mean(uout$ogeno == rout$ogeno, na.rm = TRUE)
boxplot(uout$postmean ~ rout$ogeno)

write.table(matrix(osize), "./Output/blischak_formatted_data_sims/size.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(ocounts), "./Output/blischak_formatted_data_sims/counts.txt", row.names = FALSE, col.names = FALSE)
write.table(matrix(seq_error), "./Output/blischak_formatted_data_sims/seq_error.txt", row.names = FALSE, col.names = FALSE)
command_text <- paste0("ebg hwe -p ", ploidy, " -n ", length(ocounts), " -l ", 1,
                       " -r ./Output/blischak_formatted_data_sims/counts.txt ",
                       "-t ./Output/blischak_formatted_data_sims/size.txt ",
                       "-e ./Output/blischak_formatted_data_sims/seq_error.txt")
system(command_text)
system("mv hwe* ./Output/blischak_formatted_data_sims/")
bgeno <- read.table("./Output/blischak_formatted_data_sims/hwe-genos.txt")[, 1]
ballele_freq <- read.table("./Output/blischak_formatted_data_sims/hwe-freqs.txt")[1, 1]

mean(bgeno == rout$ogeno, na.rm = TRUE)
mean(uout$ogeno == rout$ogeno, na.rm = TRUE)

mean(abs(uout$postmean - uout$ogeno))

plot_geno(ocounts, osize, ploidy, ogeno = uout$ogeno, bias_val = uout$bias_val, seq_error = uout$seq_error)
plot_geno(ocounts, osize, ploidy, ogeno = bgeno)

allele_freq
uout$allele_freq
ballele_freq


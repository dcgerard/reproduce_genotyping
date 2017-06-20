## Argument for outliers

library(updog)
library(tidyverse)
osize <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)[, 3]
ocounts <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)[, 3]
ploidy <- 6

uout <- updog_vanilla(ocounts = ocounts, osize = osize, ploidy = ploidy)
plot(uout, plot_beta = FALSE)

pvec <- updog::get_pvec(ploidy = ploidy, bias_val = uout$bias_val,
                        seq_error = uout$seq_error)

outcount <- ocounts[which.min(uout$prob_ok)]
outsize <- osize[which.min(uout$prob_ok)]
pvalue <- rmutil::pbetabinom(q = outcount, size = outsize, m = pvec[5],
                             s = (1 - uout$od_param) / uout$od_param)
pvalue
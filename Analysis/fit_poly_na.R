## Look at distribution of prop_mis

suppressMessages(library(updog))
suppressMessages(library(tidyverse))
nsnps <- 1000

fp_prop_mis <- rep(NA, length = nsnps)
u_prop_mis  <- rep(NA, length = nsnps)
for (index in seq_len(nsnps)) {
  fout <- readRDS(paste0("./Output/fp_fits/fpout", index, ".RDS"))
  uout <- readRDS(paste0("./Output/updog_fits/uout", index, ".RDS"))
  if (!is.null(fout$prop_miss)) {
    fp_prop_mis[index] <- fout$prop_miss
  }
  u_prop_mis[index] <- uout$prop_mis
}

pl <- qplot(u_prop_mis, fp_prop_mis, size = I(0.3)) +
  theme_bw() +
  xlab("updog") +
  ylab("fitPoly") +
  geom_abline(color = "red", lty = 2)

pdf(file = "./Output/fig/prop_mis_shir.pdf", height = 2.4, width = 2.4)
print(pl)
dev.off()

weird_snp <- which(u_prop_mis < 0.05 & fp_prop_mis > 0.2)

fout <- readRDS(paste0("./Output/fp_fits/fpout", weird_snp, ".RDS"))
uout <- readRDS(paste0("./Output/updog_fits/uout", weird_snp, ".RDS"))

uout$input$p1ref <- NULL
pl <- plot(uout)
pl <- pl + guides(color = guide_legend(title="Genotype"),
                  alpha = guide_legend(title="Maximum\nPosterior\nProbability"))
pdf(file = "./Output/fig/ufit_weird.pdf", height = 2.4, width = 3.7)
print(pl)
dev.off()

# plot_geno(refvec = uout$input$refvec,
#           sizevec = uout$input$sizevec,
#           ploidy = uout$input$ploidy,
#           geno = fout$scores$maxgeno[-1])
# plot(fout$scores$maxP[-1], uout$maxpostprob)
#
# weird_snp2 <- which(u_prop_mis > 0.145 & fp_prop_mis < 0.55)
#
# fout <- readRDS(paste0("./Output/fp_fits/fpout", weird_snp2, ".RDS"))
# uout <- readRDS(paste0("./Output/updog_fits/uout", weird_snp2, ".RDS"))
#
# pl <- plot_geno(refvec = uout$input$refvec,
#                 sizevec = uout$input$sizevec,
#                 ploidy = uout$input$ploidy)
#
# plot_geno(refvec = uout$input$refvec,
#           sizevec = uout$input$sizevec,
#           ploidy = uout$input$ploidy,
#           geno = fout$scores$maxgeno[-1])
#
# plot_geno(refvec = uout$input$refvec,
#           sizevec = uout$input$sizevec,
#           ploidy = uout$input$ploidy,
#           geno = uout$geno,
#           maxpostprob = 1 - uout$prob_outlier)
#
# plot(fout$scores$maxP[-1], uout$maxpostprob)


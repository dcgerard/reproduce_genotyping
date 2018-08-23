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

setEPS()
postscript(file = "./Output/fig/prop_mis_shir.eps",
           height = 2.4,
           width = 2.4,
           colormodel = "cmyk",
           paper = "special",
           horizontal = FALSE)
print(pl)
dev.off()

weird_snp <- which(u_prop_mis < 0.05 & fp_prop_mis > 0.2)

fout <- readRDS(paste0("./Output/fp_fits/fpout", weird_snp, ".RDS"))
uout <- readRDS(paste0("./Output/updog_fits/uout", weird_snp, ".RDS"))


maxcount <- max(max(uout$input$refvec, na.rm = TRUE), max(uout$input$sizevec - uout$input$refvec, na.rm = TRUE))
dfdat <- data_frame(A = uout$input$refvec, a = uout$input$sizevec - uout$input$refvec,
                    ogeno = factor(uout$geno, levels = 0:uout$input$ploidy),
                    prob_ok = 1 - uout$prob_outlier)
dfdat$snp <- paste0("SNP", index)
pk <- updog:::xi_fun(p = (0:uout$input$ploidy) / uout$input$ploidy, h = uout$bias, eps = uout$seq)
slopevec <- pk/(1 - pk)
xend <- pmin(rep(maxcount, uout$input$ploidy + 1), maxcount/slopevec)
yend <- pmin(rep(maxcount, uout$input$ploidy + 1), maxcount * slopevec)
df_lines <- data_frame(x = rep(0, uout$input$ploidy + 1), y = rep(0, uout$input$ploidy + 1), xend = xend, yend = yend)

uout$input$p1ref <- NULL
pl <- plot_geno(uout,
                refvec = uout$input$refvec,
                sizevec = uout$input$sizevec,
                ploidy = uout$input$ploidy,
                geno = uout$geno,
                seq = uout$seq,
                bias = uout$bias,
                use_colorblind = TRUE)
pl <- pl + guides(color = guide_legend(title="Genotype")) +
  geom_segment(data = df_lines,
               mapping = aes(x = x, y = y, xend = xend, yend = yend),
               lty = 2,
               color = "grey50",
               size = 0.5)
setEPS()
postscript(file = "./Output/fig/ufit_weird.eps",
           height = 2.4,
           width = 3.7,
           colormodel = "cmyk",
           paper = "special",
           horizontal = FALSE)
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


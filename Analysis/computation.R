## Computation time plots of updog
suppressMessages(library(tidyverse))
nsnps <- 1000
time_vec <- rep(NA, length = nsnps)
for (index in 1:nsnps) {
  uout <- readRDS(paste0("./Output/updog_fits/uout", index, ".RDS"))
  time_vec[index] <- uout$time[3]
}

pl <- ggplot(data = data_frame(time = time_vec), mapping = aes(x = time)) +
  geom_histogram(bins = 20, fill = "white", color = "black") +
  theme_bw() +
  xlab("Time (s)")

pdf(file = "./Output/fig/comp_time.pdf", family = "Times", colormodel = "cmyk", width = 6, height = 2)
print(pl)
dev.off()

fptimevec <- rep(NA, length = nsnps)
for (index in 1:nsnps) {
  fpout <- readRDS(paste0("./Output/fp_fits/fpout", index, ".RDS"))
  if (any(!is.na(fpout$scores))) {
    fptimevec[index] <- fpout$time[3]
  }
}

cat(paste0("updog 95% Interval of Computation Time: (", paste(round(quantile(time_vec, probs = c(0.025, 0.975)), digits = 2), collapse = ", "), ")\n\n"),
    paste0("fitPoly 95% Interval of Computation Time: (", paste(round(quantile(fptimevec, probs = c(0.025, 0.975), na.rm = TRUE), digits = 2), collapse = ", "), ")"),
    file = "./Output/text/computation_95.txt")

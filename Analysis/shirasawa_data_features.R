library(updog)
suppressMessages(library(tidyverse))
numfiles <- length(list.files("./Output/updog_fits/"))
parmat <- matrix(NA, nrow = numfiles, ncol = 5)
colnames(parmat) <- c("bias", "seq", "od", "pgeno", "out_prop")
for (index in 1:numfiles) {
  uout <- readRDS(paste0("./Output/updog_fits/uout", index, ".RDS"))
  parmat[index, 1] <- uout$bias
  parmat[index, 2] <- uout$seq
  parmat[index, 3] <- uout$od
  parmat[index, 4] <- uout$par$pgeno
  parmat[index, 5] <- uout$out_prop
}
pardat <- as_data_frame(parmat)
saveRDS(object = pardat, file = "./Output/shirasawa_snps/shir_features.RDS")

longdat <- pardat %>% mutate(logbias = log2(bias)) %>%
  select("Overdispersion" = od, "Log2-Bias" = logbias, "Sequencing Error" = seq, "Outlier Proportion" = out_prop) %>%
  gather(key = "Parameter", value = "Value")
pl <- ggplot(data = longdat, mapping = aes(x = Value)) +
  geom_histogram(bins = 20, fill = "white", color = "black") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Estimates") +
  ylab("Counts") +
  facet_wrap(~ Parameter, scales = "free")

pdf(file = "./Output/fig/ufit_features.pdf", colormodel = "cmyk", family = "Times",
    height = 6.5, width = 6.5)
print(pl)
dev.off()

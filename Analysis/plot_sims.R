## Plot the simulation results -------------------------------
suppressMessages(library(tidyverse))
dat <- as_data_frame(read.csv("./Output/sims_out/sims_out.csv", row.names = NULL))

longdat <- dat %>% select(bham, uham, od_param, bias_val) %>%
  gather(key = "Method", value = "PropCorrect", bham:uham)

pl <- ggplot(data = longdat, mapping = aes(y = PropCorrect, x = Method)) +
  facet_grid(bias_val ~ od_param) +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))

pdf(file = "./Output/fig/prop_correct_box.pdf", family = "Times", colormodel = "cmyk")
print(pl)
dev.off()

longdat <- dat %>% select(bham, uham, allele_freq, od_param, bias_val) %>%
  gather(key = "Method", value = "PropCorrect", bham:uham)

pl <- ggplot(data = longdat, mapping = aes(y = PropCorrect, x = allele_freq, color = Method)) +
  facet_grid(bias_val ~ od_param) +
  geom_smooth() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 1/2) +
  xlab("Allele Frequency") +
  ylab("Proportion Correct")
pdf(file = "./Output/fig/prop_correct_lines.pdf", family = "Times", colormodel = "cmyk")
print(pl)
dev.off()

longdat <- dat %>% select(ballele_freq, uallele_freq, allele_freq, od_param, bias_val) %>%
  gather(key = "Method", value = "est_allele_freq", ballele_freq:uallele_freq)

pl <- ggplot(data = longdat, mapping = aes(y = est_allele_freq, x = allele_freq, color = Method)) +
  facet_grid(bias_val ~ od_param) +
  geom_smooth() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 1/2) +
  xlab("Allele Frequency") +
  ylab("Estimated Allele Frequency")

pdf(file = "./Output/fig/allele_freq_est.pdf", family = "Times", colormodel = "cmyk")
print(pl)
dev.off()

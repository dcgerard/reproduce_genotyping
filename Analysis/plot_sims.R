## Plot the simulation results -------------------------------
suppressMessages(library(tidyverse))
dat <- as_data_frame(read.csv("./Output/sims_out/sims_out.csv", row.names = NULL))

longdat <- dat %>% transmute(Blischak = bham, updog = uham, od_param = od_param, bias_val = bias_val) %>%
  gather(key = "Method", value = "PropCorrect", Blischak:updog)

pl <- ggplot(data = longdat, mapping = aes(y = PropCorrect, x = Method)) +
  facet_grid(bias_val ~ od_param) +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  ylab("Proportion Correct")

pdf(file = "./Output/fig/prop_correct_box.pdf", family = "Times", colormodel = "cmyk")
print(pl)
dev.off()

longdat <- dat %>% transmute(updog = uham, Blischak = bham, allele_freq = allele_freq, od_param = od_param, bias_val = bias_val) %>%
  gather(key = "Method", value = "PropCorrect", updog:Blischak)


pl <- ggplot(data = longdat, mapping = aes(y = PropCorrect, x = allele_freq, color = Method, group = Method)) +
  facet_grid(bias_val ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Proportion Correct") +
  geom_point(size = 0.1) +
  geom_smooth(color = "black")
pdf(file = "./Output/fig/prop_correct_lines.pdf", family = "Times", colormodel = "cmyk")
print(pl)
dev.off()

longdat <- dat %>% transmute(Blischak = ballele_freq, updog = uallele_freq, allele_freq = allele_freq, od_param = od_param, bias_val = bias_val) %>%
  gather(key = "Method", value = "est_allele_freq", Blischak:updog)

pl <- ggplot(data = longdat, mapping = aes(y = est_allele_freq, x = allele_freq, color = Method, group = Method)) +
  facet_grid(bias_val ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Estimated Allele Frequency") +
  geom_point(size = 0.1) +
  geom_smooth(color = "black") +
  geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 1/2)
pdf(file = "./Output/fig/allele_freq_est.pdf", family = "Times", colormodel = "cmyk")
print(pl)
dev.off()


longdat <- dat %>% select(od_param, uod_param, bias_val, ubias_val, seq_error, useq_error, allele_freq)

## very little difference in estimate of OD at different biases and allele_frequencies

pl <- ggplot(data = longdat, mapping = aes(x = allele_freq, y = uod_param)) +
  facet_grid(bias_val ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Estimated OD Parameter") +
  geom_point(size = 0.1) +
  geom_hline(mapping = aes(yintercept = od_param), lty = 2, color = "red")
pdf(file = "./Output/fig/od_v_af.pdf", family = "Times", colormodel = "cmyk", height = 8, width = 6.5)
print(pl)
dev.off()

pl <- ggplot(data = longdat, mapping = aes(x = as.factor(bias_val), y = uod_param)) +
  facet_grid(.~od_param) +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Bias Parameter") +
  ylab("Estimated OD Parameter") +
  geom_hline(mapping = aes(yintercept = od_param), lty = 2, color = "red")
pdf(file = "./Output/fig/od_v_bias.pdf", family = "Times", colormodel = "cmyk", height = 3, width = 6.5)
print(pl)
dev.off()

pl_od <- ggplot(data = longdat, mapping = aes(x = as.factor(od_param), y = uod_param)) +
  geom_boxplot(outlier.size = 0.1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(expression(tau)) +
  ylab(expression(hat(tau))) +
  geom_hline(yintercept = 0, lty = 2, alpha = 1 / 2, color = "red") +
  geom_hline(yintercept = 0.01, lty = 2, alpha = 1 / 2, color = "red") +
  geom_hline(yintercept = 0.1, lty = 2, alpha = 1 / 2, color = "red")



## Bias estimates
pl <- ggplot(data = longdat, mapping = aes(x = allele_freq, y = log2(ubias_val))) +
  facet_grid(bias_val ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Log2 Estimated Bias Parameter") +
  geom_point(size = 0.1) +
  geom_hline(mapping = aes(yintercept = log2(bias_val)), lty = 2, color = "red")
pdf(file = "./Output/fig/bias_v_af.pdf", family = "Times", colormodel = "cmyk", height = 8, width = 6.5)
print(pl)
dev.off()

pl_bias <- ggplot(data = longdat, mapping = aes(x = as.factor(od_param), y = log2(ubias_val))) +
  facet_grid(.~bias_val) +
  geom_boxplot(outlier.size = 0.1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(expression(tau)) +
  ylab(expression(log[2](hat(h)))) +
  geom_hline(mapping = aes(yintercept = log2(bias_val)), lty = 2, color = "red")


## Sequencing error rate
pl <- ggplot(data = longdat, mapping = aes(x = allele_freq, y = useq_error)) +
  facet_grid(bias_val ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Estimated Sequencing Error Rate") +
  geom_point(size = 0.1) +
  geom_hline(mapping = aes(yintercept = seq_error), lty = 2, color = "red")
pdf(file = "./Output/fig/seq_v_af.pdf", family = "Times", colormodel = "cmyk", height = 8, width = 6.5)
print(pl)
dev.off()

pl_seq <- ggplot(data = longdat, mapping = aes(x = as.factor(od_param), y = useq_error)) +
  geom_boxplot(outlier.size = 0.1) +
  geom_hline(yintercept = 0.005, col = "red", lty = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(expression(tau)) +
  ylab(expression(hat(epsilon)))


## Plots for main part of paper.
library(gridExtra)
pdf(file = "./Output/fig/param_ests.pdf", family = "Times", colormodel = "cmyk", height = 2.3, width = 6.5)
grid.arrange(pl_od, pl_bias, pl_seq, ncol = 3)
dev.off()







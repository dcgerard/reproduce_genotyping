## Plots for output of hwe_f1_sims.R
suppressMessages(library(tidyverse))
library(ggthemes)
suppressMessages(simout <- read_csv("./Output/hwe_f1_sims/hwe_f1_sims_out.csv"))

simout %>%
  select(sham, hham, pgeno, od, bias) %>%
  gather(key = "method", value = "prop_correct", sham, hham) ->
  smalldat

pl <- ggplot(data = smalldat, mapping = aes(x = as.factor(pgeno), y = prop_correct, color = method)) +
  geom_boxplot() +
  facet_grid(bias ~ od) +
  xlab("Parent Genotype") +
  ylab("Proportion Correct") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_colorblind(name = "Prior",
                         labels = c("HWE", "S1"))

pdf(file = "./Output/fig/hwe_s1_prop_correct.pdf", colormodel = "cmyk")
print(pl)
dev.off()

simout %>%
  select(pgeno, spgeno, od, bias) %>%
  group_by(pgeno, od, bias) %>%
  summarise(prop_correct = mean(pgeno == spgeno)) ->
  dfcorrect

pl <- dfcorrect %>%
  ggplot(mapping = aes(x = pgeno, y = prop_correct)) +
  geom_line() +
  geom_point() +
  xlab("Parent Genotype") +
  ylab("Proportion of Simulations Parent\nGenotype Estimted Correctly") +
  facet_grid(bias ~ od) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = c(3, 4, 5)) +
  ylim(0, 1)

pdf(file = "./Output/fig/hwe_s1_parent_correct.pdf", colormodel = "cmyk")
print(pl)
dev.off()


## MSE
simout %>%
  select(smse, hmse, pgeno, od, bias) %>%
  gather(key = "method", value = "mse", smse, hmse) ->
  smalldat

pl <- ggplot(data = smalldat, mapping = aes(x = as.factor(pgeno), y = mse, color = method)) +
  geom_boxplot() +
  facet_grid(bias ~ od) +
  xlab("Parent Genotype") +
  ylab("MSE") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_colorblind(name = "Prior",
                         labels = c("HWE", "S1"))

pdf(file = "./Output/fig/hwe_s1_mse.pdf", colormodel = "cmyk")
print(pl)
dev.off()

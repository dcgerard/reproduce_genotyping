## Plots for run_pp_sims.R

suppressMessages(library(tidyverse))
suppressMessages(library(ggthemes))
sdat <- as_data_frame(read.csv("./Output/sims_out/sims_out_pp.csv"))

## Ability to estimate preferential pairing

sdat %>%
  select(s1pp_firstweight, firstweight, bias, od) ->
  smalldat

smalldat %>%
  ggplot(mapping = aes(x = factor(firstweight), y = s1pp_firstweight)) +
  geom_boxplot() +
  facet_grid(bias ~ od) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("True Weight") +
  ylab("Estimated Weight") ->
  pl

pdf(file = "./Output/fig/est_weight_pp.pdf",
    width = 6.5, height = 7.5,
    family = "Times")
print(pl)
dev.off()

## Proportion correct of top 3 methods

sdat %>%
  transmute(S1pp = 1 - s1pp_ham, S1 = 1 - s1_ham, fitPoly = 1 - fp_ham,
            od = od, bias = bias, weight = factor(firstweight)) %>%
  gather(key = "Method", value = "prop_correct", S1pp:fitPoly) ->
  subdat

subdat$Method <- factor(subdat$Method,levels = c("S1pp", "S1", "fitPoly"))

pl <- ggplot(data = subdat,
             mapping = aes(x = weight, y = prop_correct, color = Method)) +
  facet_grid(bias ~ od) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggthemes::scale_color_colorblind() +
  xlab("Weight on 1") +
  ylab("Proportion Correct")

pdf(file = "./Output/fig/prop_correct_pp.pdf",
    width = 6.5, height = 7.5,
    family = "Times")
print(pl)
dev.off()

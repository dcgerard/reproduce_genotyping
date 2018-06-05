## Plots output of run_sims.R

suppressMessages(library(tidyverse))
suppressMessages(library(ggthemes))
simsout <- as_data_frame(read.csv("./Output/sims_out/sims_out.csv"))

## Proportion Correct
simsout %>%
  select(bham, uham, fpham, gham, od_param, bias_val, allele_freq) %>%
  mutate(af = cut(allele_freq, breaks = seq(0, 1, by = 0.25))) %>%
  transmute(Li = bham, updog = uham, fitPoly = fpham, GATK = gham,
            od = od_param, bias = bias_val, af = af) %>%
  gather(key = "Method", value = "prop_correct", Li:GATK) ->
  subdat

subdat$Method <- factor(subdat$Method,levels = c("updog", "fitPoly", "GATK", "Li"))

pl <- ggplot(data = subdat,
             mapping = aes(x = af, y = prop_correct, color = Method)) +
  facet_grid(bias ~ od) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggthemes::scale_color_colorblind() +
  xlab("Allele Frequency") +
  ylab("Proportion Correct")

pdf(file = "./Output/fig/prop_correct.pdf",
    width = 6.5, height = 7.5,
    family = "Times")
print(pl)
dev.off()

## Correlation

simsout %>%
  select(ucor_pm, bcor, fpcor_pm, gcor_pm, naive_cor, od_param, bias_val, allele_freq) %>%
  mutate(af = cut(allele_freq, breaks = seq(0, 1, by = 0.25))) %>%
  transmute(Li = bcor, updog = ucor_pm, fitPoly = fpcor_pm,
            GATK = gcor_pm, Naive = naive_cor,
            od = od_param, bias = bias_val, af = af) %>%
  gather(key = "Method", value = "Correlation", Li:Naive) ->
  subdat

subdat$Method <- factor(subdat$Method,levels = c("updog", "fitPoly", "Naive", "GATK", "Li"))


pl <- ggplot(data = subdat,
             mapping = aes(x = af, y = Correlation, color = Method)) +
  facet_grid(bias ~ od, scales = "free_y") +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggthemes::scale_color_colorblind() +
  xlab("Allele Frequency") +
  ylab("Correlation")

pdf(file = "./Output/fig/corr.pdf",
    width = 6.5, height = 7.5,
    family = "Times")
print(pl)
dev.off()

## Estimating Proportion Misclassified
simsout %>%
  select(uepm, fpepm, od_param, bias_val, allele_freq) %>%
  mutate(af = cut(allele_freq, breaks = seq(0, 1, by = 0.25))) %>%
  transmute(updog = uepm, fitPoly = fpepm,
            od = od_param, bias = bias_val, af = af) %>%
  gather(key = "Method", value = "epm", updog, fitPoly) ->
  epmdat

simsout %>%
  select(uham, fpham, od_param, bias_val, allele_freq) %>%
  mutate(af = cut(allele_freq, breaks = seq(0, 1, by = 0.25))) %>%
  transmute(updog = 1 - uham, fitPoly = 1 - fpham,
            od = od_param, bias = bias_val, af = af) %>%
  gather(key = "Method", value = "prop_miss", updog:fitPoly) ->
  pcdat

pcdat$epm <- epmdat$epm

pcdat %>%
  mutate(diff = epm - prop_miss) ->
  pcdat

pcdat$Method <- factor(pcdat$Method, levels = c("updog", "fitPoly"))

pl <- ggplot(data = pcdat,
             mapping = aes(x = af, y = diff, color = Method)) +
  facet_grid(bias ~ od) +
  geom_boxplot(outlier.size = 0.5) +
  geom_hline(yintercept = 0, lty = 2, color = ggthemes::colorblind_pal()(3)[3]) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggthemes::scale_color_colorblind() +
  xlab("Allele Frequency") +
  ylab("Difference")


pdf(file = "./Output/fig/diff.pdf",
    width = 6.5, height = 7.5,
    family = "Times")
print(pl)
dev.off()


## Allele Frequency

simsout %>%
  transmute(Li          = ballele_freq,
            updog       = uallele_freq,
            allele_freq = allele_freq,
            od          = od_param,
            bias        = bias_val) %>%
  gather(key = "Method", value = "est_allele_freq", Li:updog) ->
  longdat



longdat %>% ggplot(mapping = aes(y     = est_allele_freq,
                                 x     = allele_freq,
                                 color = Method,
                                 group = Method)) +
  facet_grid(bias ~ od) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Estimated Allele Frequency") +
  geom_point(size = 0.1) +
  geom_smooth(color = ggthemes::colorblind_pal()(3)[3]) +
  geom_abline(intercept = 0, slope = 1, lty = 2, alpha = 1/2) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  ggthemes::scale_color_colorblind() ->
  pl

pdf(file = "./Output/fig/allele_freq_est.pdf",
    family = "Times",
    height = 7.5,
    width = 6.5)
print(pl)
dev.off()

## Rest of parameter estimates -------------------------------------


simsout %>%
  select(od_param,
         uod_param,
         bias_val,
         ubias_val,
         seq_error,
         useq_error,
         allele_freq) ->
  longdat

## very little difference in estimate of OD at different biases and allele_frequencies

pl <- ggplot(data = longdat,
             mapping = aes(x = allele_freq,
                           y = uod_param)) +
  facet_grid(bias_val ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Estimated OD Parameter") +
  geom_point(size = 0.1) +
  geom_hline(mapping = aes(yintercept = od_param), lty = 2, color = "red")
pdf(file = "./Output/fig/od_v_af.pdf",
    family = "Times",
    colormodel = "cmyk",
    height = 7.5, width = 6.5)
print(pl)
dev.off()

pl <- ggplot(data = longdat,
             mapping = aes(x = as.factor(bias_val),
                           y = uod_param)) +
  facet_grid(.~od_param) +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Bias Parameter") +
  ylab("Estimated OD Parameter") +
  geom_hline(mapping = aes(yintercept = od_param), lty = 2, color = "red")
pdf(file = "./Output/fig/od_v_bias.pdf",
    family = "Times",
    colormodel = "cmyk",
    height = 3,
    width = 6.5)
print(pl)
dev.off()

pl_od <- ggplot(data = longdat,
                mapping = aes(x = as.factor(od_param),
                              y = uod_param)) +
  geom_boxplot(outlier.size = 0.1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(expression(tau)) +
  ylab(expression(hat(tau))) +
  geom_hline(yintercept = 0, lty = 2, alpha = 1 / 2, color = "red") +
  geom_hline(yintercept = 0.01, lty = 2, alpha = 1 / 2, color = "red") +
  geom_hline(yintercept = 0.05, lty = 2, alpha = 1 / 2, color = "red")



## Bias estimates
pl <- ggplot(data = longdat,
             mapping = aes(x = allele_freq,
                           y = log2(ubias_val))) +
  facet_grid(bias_val ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Log2 Estimated Bias Parameter") +
  geom_point(size = 0.1) +
  geom_hline(mapping = aes(yintercept = log2(bias_val)), lty = 2, color = "red")
pdf(file = "./Output/fig/bias_v_af.pdf",
    family = "Times",
    colormodel = "cmyk",
    height = 7.5,
    width = 6.5)
print(pl)
dev.off()

pl_bias <- ggplot(data = longdat,
                  mapping = aes(x = as.factor(od_param),
                                y = log2(ubias_val))) +
  facet_grid(.~bias_val) +
  geom_boxplot(outlier.size = 0.1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(expression(tau)) +
  ylab(expression(log[2](hat(h)))) +
  geom_hline(mapping = aes(yintercept = log2(bias_val)), lty = 2, color = "red")


## Sequencing error rate
pl <- ggplot(data = longdat,
             mapping = aes(x = allele_freq,
                           y = useq_error)) +
  facet_grid(bias_val ~ od_param) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Allele Frequency") +
  ylab("Estimated Sequencing Error Rate") +
  geom_point(size = 0.1) +
  geom_hline(mapping = aes(yintercept = seq_error), lty = 2, color = "red")
pdf(file = "./Output/fig/seq_v_af.pdf",
    family = "Times",
    colormodel = "cmyk",
    height = 7.5,
    width = 6.5)
print(pl)
dev.off()

pl_seq <- ggplot(data = longdat,
                 mapping = aes(x = as.factor(od_param),
                               y = useq_error)) +
  geom_boxplot(outlier.size = 0.1) +
  geom_hline(yintercept = 0.005, col = "red", lty = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(expression(tau)) +
  ylab(expression(hat(epsilon)))


## Plots for main part of paper.
suppressMessages(library(gridExtra))
pdf(file = "./Output/fig/param_ests.pdf",
    family = "Times",
    colormodel = "cmyk",
    height = 2.3,
    width = 6.5)
grid.arrange(pl_od, pl_bias, pl_seq, ncol = 3)
dev.off()

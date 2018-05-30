## Plots for oracle estimator

suppressMessages(library(updog))
suppressMessages(library(tidyverse))
odat <- readRDS("./Output/oracle/odat.RDS")

## Misclassification error rate --------------------------------
## Function to extract n that has at max pmiss <= err
get_cutoff <- function(n, pmiss, err = 0.05) {
  if (all(pmiss > err)) {
    return(NA)
  } else {
    return(min(n[pmiss <= err]))
  }
}

maxerr <- 0.05
odat %>%
  group_by(ploidy, seq, bias, od, alpha) %>%
  summarize(min_n = get_cutoff(n = n, pmiss = pmiss, err = maxerr)) %>%
  ungroup() ->
  sumdat

ploidy_vec <- unique(sumdat$ploidy)
alpha_vec <- unique(sumdat$alpha)

for(ploidy_index in 1:length(ploidy_vec)) {
  for(alpha_index in 1:length(alpha_vec)) {
    alpha_current <- alpha_vec[alpha_index]
    ploidy_current <- ploidy_vec[ploidy_index]
    sumdat %>%
      filter(ploidy == ploidy_current, alpha == alpha_current) %>%
      ggplot(mapping = aes(x = bias, y = od, fill = log(min_n))) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "steelblue") +
      theme_bw() +
      geom_text(aes(label = min_n)) +
      ylab("Overdispersion") +
      xlab("Bias") +
      ggtitle(paste0("Ploidy: ", ploidy_current, ", Allele Freq: ", alpha_current)) +
      theme(legend.position = "none") ->
      pl
    pdf(file = paste0("./Output/oracle/fig/n_for_05_alpha", alpha_current * 100, "_ploidy", ploidy_current, ".pdf"), height = 4, width = 5, family = "Times")
    print(pl)
    dev.off()

  }
}

## Correlation -------------------------------------------------
## Function to extract n that has at max pmiss <= err
get_cutoff <- function(n, cor, mincor = 0.9) {
  if (all(cor < mincor)) {
    return(NA)
  } else {
    return(min(n[cor >= mincor]))
  }
}

mincor <- 0.9
odat %>%
  group_by(ploidy, seq, bias, od, alpha) %>%
  summarize(min_n = get_cutoff(n = n, cor = cor, mincor = mincor)) %>%
  ungroup() ->
  sumdat

ploidy_vec <- unique(sumdat$ploidy)
alpha_vec <- unique(sumdat$alpha)

for(ploidy_index in 1:length(ploidy_vec)) {
  for(alpha_index in 1:length(alpha_vec)) {
    alpha_current <- alpha_vec[alpha_index]
    ploidy_current <- ploidy_vec[ploidy_index]
    sumdat %>%
      filter(ploidy == ploidy_current, alpha == alpha_current) %>%
      ggplot(mapping = aes(x = bias, y = od, fill = log(min_n))) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "steelblue") +
      theme_bw() +
      geom_text(aes(label = min_n)) +
      ylab("Overdispersion") +
      xlab("Bias") +
      ggtitle(paste0("Ploidy: ", ploidy_current, ", Allele Freq: ", alpha_current)) +
      theme(legend.position = "none") ->
      pl
    pdf(file = paste0("./Output/oracle/fig/cor_n_for_05_alpha", alpha_current * 100, "_ploidy", ploidy_current, ".pdf"), height = 4, width = 5, family = "Times")
    print(pl)
    dev.off()

  }
}

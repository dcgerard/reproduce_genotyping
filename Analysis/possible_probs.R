## Plot the possible probabilities given the sequencing error rate and
## read-mapping bias.

suppressMessages(library(tidyverse))
library(updog)
ploidy <- 6
bias_vec <- c(0.1, 0.5, 1, 2, 10)
seq_error_vec <- c(0.01, 0.05, 0.1)
porig <- as.factor(paste0(0:ploidy, "/", ploidy))

for (bindex in 1:length(bias_vec)) {
    bias <- bias_vec[bindex]
    for (sindex in 1:length(seq_error_vec)) {
        seq_error <- seq_error_vec[sindex]
        pvec <- updog:::xi_fun(p   = (0:ploidy) / ploidy,
                               h   = bias,
                               eps = seq_error)
        slopevec <- pvec / (1 - pvec)
        xend <- pmin(rep(1, ploidy + 1), 1 / slopevec)
        yend <- pmin(rep(1, ploidy + 1), 1 * slopevec)
        df_lines <- data.frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1),
                               xend = xend, yend = yend, porig = porig)
        df_lines$bias <- bias
        df_lines$seq_error <- seq_error
        if (bindex == 1 & sindex == 1) {
            df_tot <- df_lines
        } else {
            df_tot <- bind_rows(df_tot, df_lines)
        }

    }
}

pl <- ggplot(data = df_tot, mapping = aes(x = x, y = y,
                                            xend = xend, yend = yend,
                                            color = porig)) +
    geom_segment() +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_rect(fill = "white")) +
    facet_grid(seq_error ~ bias) +
    xlab("Bias") +
    ylab("Sequencing Error") +
    ggthemes::scale_color_colorblind(name = "Original\nProbabilities")

setEPS()
postscript(file = "./Output/fig/prob_plots.eps",
           family = "Times",
           colormodel = "cmyk",
           width = 6,
           height = 3,
           paper = "special",
           horizontal = FALSE)
print(pl)
dev.off()


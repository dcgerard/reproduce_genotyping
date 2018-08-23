## Plot 0.025, 0.5, and 0.975 quantiles of prior on sequencing error
## rate and bias parameter.

library(updog)
suppressMessages(library(tidyverse))
ploidy <- 6

bias_lower <- exp(stats::qnorm(0.025, mean = 0, sd = 0.7))
bias_upper <- exp(stats::qnorm(0.975, mean = 0, sd = 0.7))

seq_error_lower  <- updog:::expit(stats::qnorm(0.025, mean = -4.7, sd = 1))
seq_error_median <- updog:::expit(-4.7)
seq_error_upper <- updog:::expit(stats::qnorm(0.975, mean = -4.7, sd = 1))

get_df_lines <- function(seq_error, bias_val, ploidy) {
    porig <- as.factor(paste0(0:ploidy, "/", ploidy))
    pvec <- updog:::xi_fun(p = (0:ploidy) / ploidy,
                           h = bias_val,
                           eps = seq_error)
    slopevec <- pvec / (1 - pvec)
    xend <- pmin(rep(1, ploidy + 1), 1 / slopevec)
    yend <- pmin(rep(1, ploidy + 1), 1 * slopevec)
    df_lines <- data.frame(x = rep(0, ploidy + 1), y = rep(0,
        ploidy + 1), xend = xend, yend = yend, porig = porig)
    df_lines$bias <- format(bias_val, digits = 2)
    df_lines$seq_error <- format(seq_error, digits = 2)
    return(df_lines)
}

df_tot <- bind_rows(get_df_lines(seq_error = seq_error_lower, bias_val = 1, ploidy = 6),
                    get_df_lines(seq_error = seq_error_upper, bias_val = 1, ploidy = 6),
                    get_df_lines(seq_error = 0.01, bias_val = bias_lower, ploidy = 6),
                    get_df_lines(seq_error = 0.01, bias_val = bias_upper, ploidy = 6))
df_tot$bias_seq <- paste0("(", df_tot$bias, ", ", df_tot$seq_error, ")")
df_tot$description <- c(rep("Seq Lower", ploidy + 1),
                        rep("Seq Upper", ploidy + 1),
                        rep("Bias Lower", ploidy + 1),
                        rep("Bias Upper", ploidy + 1))

pl <- ggplot(data = df_tot, mapping = aes(x = x, y = y, xend = xend, yend = yend, col = porig)) +
    facet_wrap(~ description) +
    geom_segment() +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
        ggthemes::scale_color_colorblind(name = "Original\nProbabilities")

setEPS()
postscript(file = "./Output/fig/prior_quantiles.eps",
           colormodel = "cmyk",
           family = "Times",
           height = 3,
           width = 4,
           paper = "special",
           horizontal = FALSE)
print(pl)
dev.off()

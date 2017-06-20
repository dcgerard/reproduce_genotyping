## Arguement for read-mapping bias using mendelian segregation

library(updog)
library(tidyverse)
osize <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)[, 2]
ocounts <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)[, 2]
col <- ocounts / osize > 0.95
ploidy <- 6

possible_probs <- 0:ploidy / ploidy
slopevec <- possible_probs / (1 - possible_probs)
maxcount <- max(c(ocounts, osize - ocounts), na.rm = TRUE)
xend <- pmin(rep(maxcount, ploidy + 1), maxcount / slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
df_lines <- data.frame(x = rep(0, ploidy + 1), y = rep(0,
        ploidy + 1), xend = xend, yend = yend)

dat <- data_frame(A = ocounts, a = osize - ocounts, col = col)
pl <- ggplot(data = dat, mapping = aes(x = a, y = A, color = col)) +
    geom_point(size = 0.5) +
    theme_bw() +
    xlim(0, maxcount) +
    ylim(0, maxcount) +
    ylab("Counts A") +
    xlab("Counts a") +
    geom_segment(data = df_lines,
        mapping = ggplot2::aes_string(x = "x", y = "y", xend = "xend",
        yend = "yend"), lty = 2, alpha = 1 / 2, color = "black",
        size = 0.5) +
    ggthemes::scale_color_colorblind(guide = FALSE)

pdf(file = "./Output/fig/bias_arg.pdf", family = "Times",
    colormodel = "cmyk", height = 3, width = 3)
print(pl)
dev.off()

mean(col, na.rm = TRUE)

pvalue <- pbinom(sum(col, na.rm = TRUE), size = sum(!is.na(col)), prob = 0.25)
pvalue
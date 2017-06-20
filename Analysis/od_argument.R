## Argument for overdispersion
library(updog)
library(tidyverse)
osize <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)[, 1]
ocounts <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)[, 1]

## Remove the counts close to 1 -----------------------------
which_keep <- ocounts / osize < .99
ocounts <- ocounts[which_keep]
osize <- osize[which_keep]

## Calculate intervals --------------------------------------
maxn <- max(ocounts) * 2
ploidy <- 6

possible_probs <- 0:ploidy / ploidy

## intervals at genotype = 5 ---------------------------------
lower_seq5 <- stats::qbinom(p = 0.025, size = 1:maxn,
    prob = possible_probs[ploidy])
upper_seq5 <- stats::qbinom(p = 0.975, size = 1:maxn,
    prob = possible_probs[ploidy])
xseq_lower5 <- 1:maxn - lower_seq5
xseq_upper5 <- 1:maxn - upper_seq5
dat5 <- data_frame(Alower = lower_seq5, alower = xseq_lower5,
    Aupper = upper_seq5, aupper = xseq_upper5)

## intervals at genotype = 4 ----------------------------------
lower_seq4 <- stats::qbinom(p = 0.025, size = 1:maxn,
    prob = possible_probs[ploidy - 1])
upper_seq4 <- stats::qbinom(p = 0.975, size = 1:maxn,
    prob = possible_probs[ploidy - 1])
xseq_lower4 <- 1:maxn - lower_seq4
xseq_upper4 <- 1:maxn - upper_seq4
dat4 <- data_frame(Alower = lower_seq4, alower = xseq_lower4,
    Aupper = upper_seq4, aupper = xseq_upper4)

## Find bad counts -------------------------------------------
bad_ocounts <- (ocounts > dat5$Aupper[osize]) |
               (ocounts < dat5$Alower[osize] &
                ocounts > dat4$Aupper[osize]) |
               (ocounts < dat4$Alower[osize])

## make plot
slopevec <- possible_probs / (1 - possible_probs)
maxcount <- max(c(ocounts, osize - ocounts))
xend <- pmin(rep(maxcount, ploidy + 1), maxcount / slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
df_lines <- data.frame(x = rep(0, ploidy + 1), y = rep(0,
        ploidy + 1), xend = xend, yend = yend)

dat_counts <- data_frame(A = ocounts, a = osize - ocounts, bad = bad_ocounts)

pl <- ggplot(data = dat_counts, mapping = aes(x = a, y = A, col = bad_ocounts)) +
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
print(pl)    

colvec <- ggthemes::colorblind_pal()(4)

pl <- pl + geom_line(data = dat5, mapping = aes(x = alower, y = Alower),
    color = colvec[3], lty = 1, alpha = 1 / 2)
pl <- pl + geom_line(data = dat5, mapping = aes(x = aupper, y = Aupper),
    color = colvec[3], lty = 1, alpha = 1 / 2)
pl <- pl + geom_line(data = dat4, mapping = aes(x = alower, y = Alower),
    color = colvec[4], lty = 1, alpha = 1 / 2)
pl <- pl + geom_line(data = dat4, mapping = aes(x = aupper, y = Aupper),
    color = colvec[4], lty = 1, alpha = 1 / 2)

pdf(file = "./Output/fig/od_arg.pdf", colormodel = "cmyk",
    family = "Times", height = 3, width = 3)
print(pl)
dev.off()

## Calculate p-value of seeing so many points outside the 95 percent intervals
mean(bad_ocounts)
pvalue <- stats::pbinom(q = sum(bad_ocounts), size = length(bad_ocounts), prob = 0.05, lower.tail = FALSE)

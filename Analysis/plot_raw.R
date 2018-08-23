## Plot the raw data ----------------------------------
library(updog)
suppressMessages(library(tidyverse))

ploidy      <- 6
osize_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)[, 1:3]
ocounts_mat <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)[, 1:3]

l1 <- gather(osize_mat, key = "SNP", value = "Size")
l2 <- gather(ocounts_mat, key = "SNP", value = "Counts")
longdat <- l1
longdat$Counts <- l2$Counts
longdat$A <- longdat$Counts
longdat$a <- longdat$Size - longdat$Counts
possible_probs <- 0:ploidy / ploidy
maxcount <- max(max(longdat$A, na.rm = TRUE), max(longdat$a, na.rm = TRUE), na.rm = TRUE)
slopevec <- possible_probs / (1 - possible_probs)
xend <- pmin(rep(maxcount, ploidy + 1), maxcount / slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
smalldat <- data_frame(SNP = rep(c("SNP1", "SNP2", "SNP3"), each = ploidy + 1),
    xend = rep(xend, times = 3), yend = rep(yend, times = 3))
smalldat$xstart <- 0
smalldat$ystart <- 0
smalldat$out <- FALSE
smalldat$ok_size <- 0.5

longdat$out <- FALSE
longdat$ok_size <- 0.5
pos <- which.max(filter(longdat, SNP == "SNP3")$a)
longdat$out[longdat$SNP == "SNP3"][pos] <- TRUE
longdat$ok_size[longdat$SNP == "SNP3"][pos] <- 1

pl <- ggplot(data = longdat, mapping = aes(x = a, y = A, pch = out, size = ok_size)) +
    geom_point() +
    theme_bw() +
    facet_grid(. ~ SNP) +
    xlim(0, maxcount) +
    ylim(0, maxcount) +
    geom_segment(data = smalldat, mapping = aes(x = xstart, y = ystart, xend = xend, yend = yend), lty = 2, color = "grey50") +
    theme(strip.background = element_rect(fill = "white"),
          legend.position = "none") +
    xlab("Counts a") +
    ylab("Counts A") +
    scale_size_continuous(range = c(0.4, 1))

setEPS()
postscript(file = "./Output/fig/snp_examples.eps",
           colormodel = "cmyk",
           family = "Times",
           height = 2.2,
           width = 6,
           paper = "special",
           horizontal = FALSE)
print(pl)
dev.off()

## Single well-behaved SNP -------------------------------------------------------------------

osize_vec   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)[, 7]
ocounts_vec <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)[, 7]

dfdat <- data_frame(a = osize_vec - ocounts_vec, A = ocounts_vec)
maxcount <- max(c(dfdat$a, dfdat$A))
xend <- pmin(rep(maxcount, ploidy + 1), maxcount / slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
smalldat <- data_frame(SNP = rep(c("SNP1", "SNP2", "SNP3"), each = ploidy + 1),
                       xend = rep(xend, times = 3), yend = rep(yend, times = 3))
smalldat$xstart <- 0
smalldat$ystart <- 0
smalldat$out <- FALSE
smalldat$ok_size <- 0.5

pl <- ggplot(data = dfdat, mapping = aes(x = a, y = A)) +
  geom_point(size = 0.2) +
  theme_bw() +
  xlim(0, maxcount) +
  ylim(0, maxcount) +
  theme(strip.background = element_rect(fill = "white"),
        legend.position = "none") +
  geom_segment(data = smalldat, mapping = aes(x = xstart, y = ystart, xend = xend, yend = yend), lty = 2, color = "grey50") +
  xlab("Counts a") +
  ylab("Counts A") +
  scale_size_continuous(range = c(0.4, 1))

setEPS()
postscript(file = "./Output/fig/single_snp.eps",
           colormodel = "cmyk",
           family = "Times",
           height = 2.2,
           width = 2.4,
           paper = "special",
           horizontal = FALSE)
print(pl)
dev.off()

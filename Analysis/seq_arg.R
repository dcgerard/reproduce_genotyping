## Argument for modeling sequencing error
library(updogAlpha)
suppressMessages(library(tidyverse))
library(ggthemes)

ploidy      <- 6
osize_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)
ocounts_mat <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)

index <- 13
plot_geno(ocounts = ocounts_mat[, index], osize = osize_mat[, index], ploidy = ploidy)
uout1 <- updog_vanilla(ocounts = ocounts_mat[, index], osize = osize_mat[, index],
                       ploidy = ploidy, update_outmean = FALSE, update_bias_val = FALSE,
                       update_outdisp = FALSE, update_od_param = FALSE, update_outprop = FALSE,
                       update_pgeno = TRUE, update_seq_error = FALSE, od_param = 10^-6,
                       out_prop = 0,
                       seq_error = 0, seq_error_mean = -Inf)
uout2 <- updog_vanilla(ocounts = ocounts_mat[, index], osize = osize_mat[, index],
                       ploidy = ploidy, update_outmean = FALSE, update_bias_val = FALSE,
                       update_outdisp = FALSE, update_od_param = FALSE, update_outprop = FALSE,
                       update_pgeno = TRUE, update_seq_error = TRUE, od_param = 10 ^ -6, out_prop = 0)


## Make the plot
ploidy <- 6
dfdat <- data_frame(A = c(uout1$input$ocounts, uout2$input$ocounts),
                    a = c(uout1$input$osize - uout1$input$ocounts, uout2$input$osize - uout2$input$ocounts),
                    genotype = factor(c(uout1$ogeno, uout2$ogeno)),
                    fit = c(rep("Seq Error Set to 0", length(uout2$input$ocounts)), rep("Seq Error Estimated", length(uout1$input$ocounts))))

possible_probs <- 0:ploidy / ploidy
slopevec <- possible_probs / (1 - possible_probs)
maxcount <- max(c(dfdat$a, dfdat$A))
xend <- pmin(rep(maxcount, ploidy + 1), maxcount / slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
smalldat <- data_frame(SNP = rep(c("SNP1", "SNP2", "SNP3"), each = ploidy + 1),
                       xend = rep(xend, times = 3), yend = rep(yend, times = 3))
smalldat$xstart <- 0
smalldat$ystart <- 0
smalldat$out <- FALSE
smalldat$ok_size <- 0.5

dfdat$bad_points <- (dfdat$A / (dfdat$A + dfdat$a) > 0.98 & dfdat$genotype == 5) * 1


pal <- colorblind_pal()(3)[c(2,1,3)]
pl <- ggplot(data = dfdat, mapping = aes(x = a, y = A)) +
  geom_point(aes(colour = genotype, size = bad_points, pch = as.factor(bad_points))) +
  theme_bw() +
  facet_grid(. ~ fit) +
  xlim(0, maxcount) +
  ylim(0, maxcount) +
  geom_segment(data = smalldat, mapping = aes(x = xstart, y = ystart, xend = xend, yend = yend), lty = 2, alpha = 1 / 10) +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Counts a") +
  ylab("Counts A") +
  scale_color_manual(values = pal) +
  scale_size_continuous(range = c(0.4, 1), guide = FALSE) +
  scale_shape_discrete(guide = FALSE)

pdf(file = "./Output/fig/seq_error_example.pdf", colormodel = "cmyk",
    family = "Times", height = 2.4, width = 5)
print(pl)
dev.off()

sum(dfdat$bad_points)

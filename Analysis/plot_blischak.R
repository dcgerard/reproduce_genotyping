## Plot results from Blischak fits
library(updog)
suppressMessages(library(tidyverse))
geno_mat <- read.table("./Output/blischak_formatted_data/diseq-genos.txt")
size_mat <- read.table("./Output/blischak_formatted_data/size_mat.txt")
count_mat <- read.table("./Output/blischak_formatted_data/counts_mat.txt")
seq_error <- read.table("./Output/blischak_formatted_data/seq_error.txt")
ploidy <- 6


## SNP 58 is interesting to compare
for (index in 1:3) {
  ocounts <- count_mat[, index]
  osize   <- size_mat[, index]
  ogeno   <- geno_mat[, index]
  ocounts <- ocounts[ocounts != -9]
  osize <- osize[osize != -9]
  ogeno <- ogeno[ogeno != -9]

  dfdat <- data_frame(A = ocounts, a = osize - ocounts, geno = factor(ogeno, levels = 0:ploidy))
  dfdat$snp <- paste0("SNP", index)


  maxcount <- max(max(count_mat[, 1:3], na.rm = TRUE), max(size_mat[, 1:3] - count_mat[, 1:3], na.rm = TRUE))
  pk <- updog:::xi_fun(p = (0:ploidy) / ploidy, h = 1, eps = seq_error[index, 1])
  slopevec <- pk/(1 - pk)
  xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
  yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
  df_lines <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
  df_lines$snp <- paste0("SNP", index)

  if (index == 1) {
    dfdat_tot <- dfdat
    df_lines_tot <- df_lines
  } else {
    dfdat_tot <- bind_rows(dfdat_tot, dfdat)
    df_lines_tot <- bind_rows(df_lines_tot, df_lines)
  }
}

pl <- ggplot(data = dfdat_tot, mapping = aes(x = a, y = A, col = geno)) +
  facet_grid(. ~ snp) +
  geom_point(size = 0.4) +
  xlim(0, maxcount) +
  ylim(0, maxcount) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  guides(color = guide_legend(title = "Genotype", keyheight = 0.15, default.unit = "inch")) +
  xlab("Counts a") +
  ylab("Counts A") +
  geom_segment(data = df_lines_tot, mapping = aes(x = x, y = y, xend = xend, yend = yend),
               lty = 2, color = "grey50", size = 0.5) +
  ggthemes::scale_color_colorblind()

setEPS()
postscript(file = "./Output/fig/blischak_fits.eps",
           family = "Times",
           colormodel = "cmyk",
           width = 6.5,
           height = 2.1,
           paper = "special",
           horizontal = FALSE)
print(pl)
dev.off()


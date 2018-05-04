## Plot supermassa results

## The fits were generated using http://statgen.esalq.usp.br/SuperMASSA/

library(updog)
suppressMessages(library(tidyverse))
ploidy <- 6
for (index in 1:3) {
  snpdat <- read.table(paste0("./Output/supermassa_formatted_data/osnp", index, ".txt"))
  snpgeno <- read.table(paste0("./Output/supermassa_formatted_data/supermassa_out", index, ".txt"))
  snpgeno[, 2] <- stringr::str_replace(snpgeno[, 2], "\\(", "")
  snpgeno[, 2] <- stringr::str_replace(snpgeno[, 2], ",", "")
  snpgeno[, 2] <- as.numeric(snpgeno[, 2])
  snpgeno[, 3] <- stringr::str_replace(snpgeno[, 3], "\\)", "")
  snpgeno[, 3] <- as.numeric(snpgeno[, 3])
  names(snpdat) <- c("ID", "A", "a")
  snpdat$geno <- factor(snpgeno[, 2], levels = 0:ploidy)
  snpdat$snp  <- paste0("SNP", index)

  maxcount <- 373
  pk <- updog:::xi_fun(p = (0:ploidy) / ploidy, h = 1, eps = 0)
  slopevec <- pk/(1 - pk)
  xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
  yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
  df_lines <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
  df_lines$snp <- paste0("SNP", index)

  if (index == 1) {
    snptot <- snpdat
    df_lines_tot <- df_lines
  } else {
    snptot <- bind_rows(snptot, snpdat)
    df_lines_tot <- bind_rows(df_lines_tot, df_lines)
  }
}

pl <- ggplot(data = snptot, mapping = aes(x = a, y = A, col = geno)) +
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
               lty = 2, alpha = 1 / 2, color = "black", size = 0.5) +
  ggthemes::scale_color_colorblind()

pdf(file = "./Output/fig/supermassa_fits.pdf", colormodel = "cmyk",
    family = "Times", height = 2.1, width = 6.5)
print(pl)
dev.off()

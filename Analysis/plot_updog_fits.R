## Plot three SNP's

library(updog)
suppressMessages(library(tidyverse))
counts_mat <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)[, 1:3]
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)[, 1:3]
ploidy <- 6

maxcount <- max(max(counts_mat, na.rm = TRUE), max(size_mat - counts_mat, na.rm = TRUE))

summaries_mat <- matrix(NA, nrow = ncol(counts_mat), ncol = 9) ## contains parent info and parameter estimates
colnames(summaries_mat) <- c("A", "a", "seq_error", "bias_val",
                             "od_param", "out_prop", "ogeno", "prob_ok", "snp")

for (index in 1:3) {
  uout <- readRDS(paste0("./Output/updog_fits/uout", index, ".RDS"))
  dfdat <- data_frame(A = uout$input$ocounts, a = uout$input$osize - uout$input$ocounts,
                      ogeno = factor(uout$ogeno, levels = 0:ploidy),
                      prob_ok = uout$prob_ok)
  dfdat$snp <- paste0("SNP", index)

  pk <- get_pvec(ploidy = ploidy, bias_val = uout$bias_val, seq_error = uout$seq_error)
  slopevec <- pk/(1 - pk)
  xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
  yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
  df_lines <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
  df_lines$snp <- paste0("SNP", index)

  summaries_mat[index, ] <- c(uout$input$p1counts, uout$input$p1size - uout$input$p1counts,
                              uout$seq_error, uout$bias_val,
                              uout$od_param, uout$out_prop, uout$p1geno, 1 - uout$p1_prob_out,
                              index)

  if(index == 1) {
    dfdat_tot <- dfdat
    df_lines_tot <- df_lines
  } else {
    dfdat_tot <- bind_rows(dfdat_tot, dfdat)
    df_lines_tot <- bind_rows(df_lines_tot, df_lines)
  }
}


summaries_df <- as_data_frame(summaries_mat)
summaries_df$snp <- paste0("SNP", summaries_df$snp)
summaries_df$ogeno <- factor(summaries_df$ogeno, levels = 0:ploidy)

pl <- ggplot(data = dfdat_tot, mapping = aes(x = a, y = A, col = ogeno, alpha = prob_ok)) +
  facet_grid(. ~ snp) +
  geom_point(size = 0.4) +
  xlim(0, maxcount) +
  ylim(0, maxcount) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  guides(color = guide_legend(title = "Genotype", keyheight = 0.15, default.unit = "inch"),
         alpha = guide_legend(title = "Probability\nNon-outlier", keyheight = 0.15, default.unit = "inch")) +
  xlab("Counts a") +
  ylab("Counts A") +
  geom_point(data = summaries_df, pch = 17, size = 0.4) +
  geom_segment(data = df_lines_tot, mapping = aes(x = x, y = y, xend = xend, yend = yend),
               lty = 2, alpha = 1 / 2, color = "black", size = 0.5) +
  ggthemes::scale_color_colorblind()

pdf(file = "./Output/fig/updog_fits.pdf", family = "Times", colormodel = "cmyk",
    width = 6.5, height = 2.1)
print(pl)
dev.off()

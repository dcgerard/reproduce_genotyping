## Combine plots from updog, blischak, and supermassa
suppressMessages(library(updog))
suppressMessages(library(tidyverse))
ploidy <- 6


## Start out with supermassa ------------------------------------------------------------------
for (index in 1:3) {
  snpdat <- read.table(paste0("./Output/supermassa_formatted_data/osnp", index, ".txt"))
  snpgeno <- read.table(paste0("./Output/supermassa_formatted_data/supermassa_out", index, ".txt"))
  snpgeno[, 2] <- stringr::str_replace(snpgeno[, 2], "\\(", "")
  snpgeno[, 2] <- stringr::str_replace(snpgeno[, 2], ",", "")
  snpgeno[, 2] <- as.numeric(snpgeno[, 2])
  snpgeno[, 3] <- stringr::str_replace(snpgeno[, 3], "\\)", "")
  snpgeno[, 3] <- as.numeric(snpgeno[, 3])
  names(snpdat) <- c("ID", "A", "a")
  snpdat$geno <- snpgeno[, 2]
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

df_lines_tot$Method <- "SuperMASSA"
snptot$Method       <- "SuperMASSA"
dflines_combo <- df_lines_tot
snp_combo     <- select(snptot, -ID)


## Now do Blischak --------------------------------------------------------------------------------
geno_mat <- read.table("./Output/blischak_formatted_data/diseq-genos.txt")
size_mat <- read.table("./Output/blischak_formatted_data/size_mat.txt")
count_mat <- read.table("./Output/blischak_formatted_data/counts_mat.txt")
seq_error <- read.table("./Output/blischak_formatted_data/seq_error.txt")
for (index in 1:3) {
  ocounts <- count_mat[, index]
  osize   <- size_mat[, index]
  ogeno   <- geno_mat[, index]
  ocounts <- ocounts[ocounts != -9]
  osize <- osize[osize != -9]
  ogeno <- ogeno[ogeno != -9]

  dfdat <- data_frame(A = ocounts, a = osize - ocounts, geno = ogeno)
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

dfdat_tot$Method    <- "Blischak"
df_lines_tot$Method <- "Blischak"

dflines_combo <- bind_rows(dflines_combo, df_lines_tot)
snp_combo     <- bind_rows(snp_combo, dfdat_tot)

snp_combo$prob_ok <- 1

## Updog method ------------------------------------------------------------------------------
counts_mat <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)[, 1:3]
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)[, 1:3]
maxcount <- max(max(counts_mat, na.rm = TRUE), max(size_mat - counts_mat, na.rm = TRUE))
summaries_mat <- matrix(NA, nrow = ncol(counts_mat), ncol = 9) ## contains parent info and parameter estimates
colnames(summaries_mat) <- c("A", "a", "seq_error", "bias_val",
                             "od_param", "out_prop", "ogeno", "prob_ok", "snp")
for (index in 1:3) {
  uout <- readRDS(paste0("./Output/updog_fits/uout", index, ".RDS"))
  dfdat <- data.frame(A       = uout$input$refvec,
                      a       = uout$input$sizevec- uout$input$refvec,
                      geno    = uout$geno,
                      prob_ok = 1 - uout$prob_outlier)
  dfdat$snp <- paste0("SNP", index)

  pk <- updog:::xi_fun(p = (0:ploidy) / ploidy, eps = uout$seq, h = uout$bias)
  slopevec <- pk/(1 - pk)
  xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
  yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
  df_lines <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
  df_lines$snp <- paste0("SNP", index)

  summaries_mat[index, ] <- c(uout$input$p1ref,
                              uout$input$p1size - uout$input$p1ref,
                              uout$seq,
                              uout$bias,
                              uout$od,
                              uout$out_prop,
                              uout$par$pgeno,
                              1,
                              index)

  if(index == 1) {
    dfdat_tot <- dfdat
    df_lines_tot <- df_lines
  } else {
    dfdat_tot <- bind_rows(dfdat_tot, dfdat)
    df_lines_tot <- bind_rows(df_lines_tot, df_lines)
  }
}

dfdat_tot$Method <- "updog"
df_lines_tot$Method <- "updog"
snp_combo <- bind_rows(snp_combo, select(dfdat_tot, A, a, geno, snp, Method, prob_ok))
dflines_combo <- bind_rows(dflines_combo, df_lines_tot)


## Finally, add fitPoly --------------------------------------------------
for (index in 1:3) {
  fpout <- readRDS(paste0("./Output/fp_fits/fpout", index, ".RDS"))
  dfdat <- data.frame(A       = counts_mat[, index],
                      a       = size_mat[, index] - counts_mat[, index],
                      geno    = fpout$scores$maxgeno,
                      prob_ok = 1)
  dfdat$snp <- paste0("SNP", index)

  pk <- unlist(c(select(fpout$modeldata, mu0:mu6)))
  slopevec <- pk/(1 - pk)
  xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
  yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
  df_lines <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
  df_lines$snp <- paste0("SNP", index)

  if(index == 1) {
    dfdat_tot <- dfdat
    df_lines_tot <- df_lines
  } else {
    dfdat_tot <- bind_rows(dfdat_tot, dfdat)
    df_lines_tot <- bind_rows(df_lines_tot, df_lines)
  }
}

dfdat_tot$Method <- "fitPoly"
df_lines_tot$Method <- "fitPoly"
snp_combo <- bind_rows(snp_combo, select(dfdat_tot, A, a, geno, snp, Method, prob_ok))
dflines_combo <- bind_rows(dflines_combo, df_lines_tot)
snp_combo$geno <- factor(snp_combo$geno, levels = 0:ploidy)

## Now plot all of them together --------------------------------------------------------------
snp_combo$geno <- factor(snp_combo$geno, levels = 0:ploidy)
possible_colors <- ggthemes::colorblind_pal()(5)[c(5, 3:1)]

pl <- ggplot(data = snp_combo, mapping = aes(x = a, y = A, color = geno, alpha = prob_ok)) +
  facet_grid(Method ~ snp) +
  geom_point(size = 0.3) +
  xlim(0, maxcount) +
  ylim(0, maxcount) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_segment(data = dflines_combo, mapping = aes(x = x, y = y, xend = xend, yend = yend),
               lty = 2, alpha = 1 / 2, color = "black", size = 0.5) +
  xlab("Counts a") +
  ylab("Counts A") +
  scale_alpha_continuous(name = "Probability\nNon-outlier") +
  scale_color_manual(name = "Estimated\nGenotype", values = possible_colors) +
  guides(colour = guide_legend(override.aes = list(size=1.5)),
         alpha = guide_legend(override.aes = list(size=1.5)))

pdf(file = "./Output/fig/real_data_plots.pdf", colormodel = "cmyk", family = "Times",
    height = 7.5, width = 6.5)
print(pl)
dev.off()



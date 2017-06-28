## Hypothetical problem.
set.seed(20)
library(updog)
suppressMessages(library(tidyverse))
uout <- readRDS("./Output/updog_fits/uout1.RDS")

uout$input$ploidy <- 4
uout$p1geno       <- 4
uout$p2geno       <- 0
uout$bias_val     <- 0.5
uout$seq_error    <- 0.005
uout$od_param     <- 0.001
uout$input$model  <- "f1"
ploidy            <- 4

rout <- rupdog(uout)

maxcount <- max(c(rout$input$ocounts, rout$input$osize - rout$input$ocounts))

unew1 <- updog_vanilla(ocounts = rout$input$ocounts, osize = rout$input$osize, ploidy = rout$input$ploidy, seq_error_sd = Inf, bias_val_sd = Inf)

unew2 <- updog_vanilla(ocounts = rout$input$ocounts, osize = rout$input$osize, ploidy = rout$input$ploidy, bias_val_sd = Inf)

unew3 <- updog_vanilla(ocounts = rout$input$ocounts, osize = rout$input$osize, ploidy = rout$input$ploidy)

## No labeling -----------------------------------------------------------------------
dfdat <- data_frame(A = rout$input$ocounts, a = rout$input$osize - rout$input$ocounts)
dfdat$geno <- NA
dfdat$type <- "A"

pk <- get_pvec(ploidy = ploidy, bias_val = 1, seq_error = 0)
slopevec <- pk/(1 - pk)
xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
df_lines <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
df_lines$type <- "A"

## No penalty on either --------------------------------------------------------------
dftemp <- data_frame(A = unew1$input$ocounts, a = unew1$input$osize - unew1$input$ocounts, geno = unew1$ogeno)
dftemp$type <- "B"
dfdat <- bind_rows(dfdat, dftemp)

pk <- get_pvec(ploidy = ploidy, bias_val = unew1$bias_val, seq_error = unew1$seq_error)
slopevec <- pk/(1 - pk)
xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
df_lines_temp <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
df_lines_temp$type <- "B"
df_lines <- bind_rows(df_lines, df_lines_temp)

## Penalty only on sequencing error rate ---------------------------------------------
dftemp <- data_frame(A = unew2$input$ocounts, a = unew2$input$osize - unew2$input$ocounts, geno = unew2$ogeno)
dftemp$type <- "C"
dfdat <- bind_rows(dfdat, dftemp)

pk <- get_pvec(ploidy = ploidy, bias_val = unew2$bias_val, seq_error = unew2$seq_error)
slopevec <- pk/(1 - pk)
xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
df_lines_temp <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
df_lines_temp$type <- "C"
df_lines <- bind_rows(df_lines, df_lines_temp)

## Penalty on both bias and sequencing error rate -----------------------------------
dftemp <- data_frame(A = unew3$input$ocounts, a = unew3$input$osize - unew3$input$ocounts, geno = unew3$ogeno)
dftemp$type <- "D"
dfdat <- bind_rows(dfdat, dftemp)

pk <- get_pvec(ploidy = ploidy, bias_val = unew3$bias_val, seq_error = unew3$seq_error)
slopevec <- pk/(1 - pk)
xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
df_lines_temp <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
df_lines_temp$type <- "D"
df_lines <- bind_rows(df_lines, df_lines_temp)

## Plot Results ----------------------------------------------------------------------
dfdat$geno <- factor(dfdat$geno, levels = 0:4)

pl <- ggplot(data = dfdat, mapping = aes(x = a, y = A, col = geno)) +
  facet_wrap(~ type) +
  geom_point(size = 0.2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Counts a") +
  ylab("Counts A") +
  guides(color = guide_legend(title = "Genotype")) +
  scale_color_hue(drop = FALSE) +
  geom_segment(data = df_lines, mapping = aes(x = x, y = y, xend = xend, yend = yend),
               lty = 2, alpha = 1 / 2, color = "black", size = 0.5)

pdf(file = "./Output/fig/ident_prob.pdf", family = "Times", colormodel = "cmyk",
    width = 4, height = 3)
print(pl)
dev.off()


cat(file = "./Output/text/hypo_text.txt",
    "Unpenalized Sequencing Error Rate:", unew1$seq_error, "\n",
    "Unpenalized Bias:", unew1$bias_val, "\n",
    "Seq error rate when only seq error rate is penalized:", unew2$seq_error, "\n",
    "Bias when only sequencing rate is penalized:", unew2$bias_val, "\n",
    "Penalized sequencing error rate:", unew3$seq_error, "\n",
    "Penalized bias:", unew3$bias_val)

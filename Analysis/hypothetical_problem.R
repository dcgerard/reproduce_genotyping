## Hypothetical problem.


suppressMessages(library(updog))
suppressMessages(library(tidyverse))
library(ggthemes)
uout <- readRDS("./Output/updog_fits/uout1.RDS")
sizevec <- uout$input$sizevec
nind <- length(sizevec)

ploidy    <- 4
p1geno    <- 4
p2geno    <- 0
bias      <- 0.8
seq_error <- 0.005
od        <- 0.01
model     <- "f1"
ploidy    <- 4

set.seed(1)
geno_vec <- rgeno(n      = nind,
                  ploidy = ploidy,
                  model  = model,
                  p1geno = p1geno,
                  p2geno = p2geno)

refvec <- rflexdog(sizevec = sizevec,
                   geno    = geno_vec,
                   ploidy  = ploidy,
                   seq     = seq_error,
                   bias    = bias,
                   od      = od)

maxcount <- max(c(refvec, sizevec - refvec))

## start at more bias values to make sure
## reach maximum likelihood in this weird situation
unew1 <- flexdog(refvec    = refvec,
                 sizevec   = sizevec,
                 ploidy    = ploidy,
                 model     = "f1",
                 var_seq   = Inf,
                 var_bias  = Inf,
                 bias_init = exp(c(-2, -1, -0.5, 0, 0.5, 1, 2)),
                 verbose   = FALSE)

unew2 <- flexdog(refvec    = refvec,
                 sizevec   = sizevec,
                 ploidy    = ploidy,
                 model     = "f1",
                 var_bias  = Inf,
                 bias_init = exp(c(-2, -1, -0.5, 0, 0.5, 1, 2)),
                 verbose   = FALSE)


unew3 <- flexdog(refvec    = refvec,
                 sizevec   = sizevec,
                 ploidy    = ploidy,
                 model     = "f1",
                 var_seq   = 0.1,
                 bias_init = exp(c(-2, -1, -0.5, 0, 0.5, 1, 2)),
                 verbose   = FALSE)


## No labeling -----------------------------------------------------------------------
dfdat <- data_frame(A = refvec, a = sizevec - refvec)
dfdat$geno <- NA
dfdat$type <- "A"

pk <- updog:::xi_fun(p = (0:ploidy) / ploidy, eps = 0, h = 1)
slopevec <- pk/(1 - pk)
xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
df_lines <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
df_lines$type <- "A"

## No penalty on either --------------------------------------------------------------
dftemp <- data_frame(A = unew1$input$refvec, a = unew1$input$sizevec - unew1$input$refvec, geno = unew1$geno)
dftemp$type <- "B"
dfdat <- bind_rows(dfdat, dftemp)

pk <- updog:::xi_fun(p = (0:ploidy) / ploidy, h = unew1$bias, eps = unew1$seq)
slopevec <- pk/(1 - pk)
xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
df_lines_temp <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
df_lines_temp$type <- "B"
df_lines <- bind_rows(df_lines, df_lines_temp)

## Penalty only on sequencing error rate ---------------------------------------------
# dftemp <- data_frame(A = unew2$input$refvec, a = unew2$input$sizevec - unew2$input$refvec, geno = unew2$geno)
# dftemp$type <- "E"
# dfdat <- bind_rows(dfdat, dftemp)
#
# pk <- updog:::xi_fun(p = (0:ploidy) / ploidy, h = unew2$bias, eps = unew2$seq)
# slopevec <- pk/(1 - pk)
# xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
# yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
# df_lines_temp <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
# df_lines_temp$type <- "E"
# df_lines <- bind_rows(df_lines, df_lines_temp)

## Penalty on both bias and sequencing error rate -----------------------------------
dftemp <- data_frame(A = unew3$input$refvec, a = unew3$input$sizevec - unew3$input$refvec, geno = unew3$geno)
dftemp$type <- "C"
dfdat <- bind_rows(dfdat, dftemp)

pk <- updog:::xi_fun(p = (0:ploidy) / ploidy, h = unew3$bias, eps = unew3$seq)
slopevec <- pk/(1 - pk)
xend <- pmin(rep(maxcount, ploidy + 1), maxcount/slopevec)
yend <- pmin(rep(maxcount, ploidy + 1), maxcount * slopevec)
df_lines_temp <- data_frame(x = rep(0, ploidy + 1), y = rep(0, ploidy + 1), xend = xend, yend = yend)
df_lines_temp$type <- "C"
df_lines <- bind_rows(df_lines, df_lines_temp)

## Plot Results ----------------------------------------------------------------------
dfdat$geno <- factor(dfdat$geno, levels = 0:4)

dfdat %>%
ggplot(mapping = aes(x = a, y = A, col = geno)) +
  facet_grid(. ~ type) +
  geom_point(size = 0.2) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Counts a") +
  ylab("Counts A") +
  scale_color_hue(drop = FALSE) +
  # ggthemes::scale_color_colorblind(drop = FALSE) +
  geom_segment(data = df_lines, mapping = aes(x = x, y = y, xend = xend, yend = yend),
               lty = 2, color = "grey50", size = 0.5) +
  guides(colour = guide_legend(override.aes = list(size=1.5),
                               title = "Genotype")) ->
  pl

setEPS()
postscript(file = "./Output/fig/ident_prob.eps",
           family = "Times",
           colormodel = "cmyk",
           width = 6.5,
           height = 2.1,
           paper = "special",
           horizontal = FALSE)
print(pl)
dev.off()


cat(file = "./Output/text/hypo_text.txt",
    "Unpenalized Sequencing Error Rate:", unew1$seq_error, "\n",
    "Unpenalized Bias:", unew1$bias_val, "\n",
    "Seq error rate when only seq error rate is penalized:", unew2$seq_error, "\n",
    "Bias when only sequencing rate is penalized:", unew2$bias_val, "\n",
    "Penalized sequencing error rate:", unew3$seq_error, "\n",
    "Penalized bias:", unew3$bias_val)

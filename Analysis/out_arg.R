## Argument for outliers
library(updog)
uout <- readRDS("./Output/updog_fits/uout3.RDS")
pvec <- updog:::xi_fun(p   = (0:uout$input$ploidy) / uout$input$ploidy, 
                       h   = uout$bias, 
                       eps = uout$seq)
outcount <- uout$input$refvec[which.max(uout$prob_outlier)]
outsize <- uout$input$sizevec[which.max(uout$prob_outlier)]
pvalue <- updog::pbetabinom(q     = outcount, 
                            size  = outsize, 
                            mu    = pvec[5], 
                            rho   = uout$od, 
                            log_p = FALSE)
cat(file = "./Output/text/out_prob.txt",
    "Probability of 'outlying point' having as few or fewer A counts given genotype of AAAAaa:", pvalue)


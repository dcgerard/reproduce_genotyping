## Argument for outliers
library(updog)
uout <- readRDS("./Output/updog_fits/uout3.RDS")
pvec <- updog::get_pvec(ploidy = uout$input$ploidy, bias_val = uout$bias_val,
                        seq_error = uout$seq_error)
outcount <- uout$input$ocounts[which.min(uout$prob_ok)]
outsize <- uout$input$osize[which.min(uout$prob_ok)]
pvalue <- rmutil::pbetabinom(q = outcount, size = outsize, m = pvec[5],
                             s = (1 - uout$od_param) / uout$od_param)
cat(file = "./Output/text/out_prob.txt",
    "Probability of 'outlying point' having as few or fewer A counts given genotype of AAAAaa:", pvalue)


suppressMessages(library(vcfR))

dat <- vcfR::read.vcfR(file = "./Data/KDRIsweetpotatoXushu18S1LG2017.vcf.gz")
alt_mat <- apply(extract.gt(x = dat, element = "AD"), 2, as.numeric)
ref_mat <- apply(extract.gt(x = dat, element = "RD"), 2, as.numeric)

var_ave <- rowMeans(ref_mat + alt_mat, na.rm = TRUE)
order_vec <- order(var_ave, decreasing = TRUE)

## Extract the SNPs that would be good case studies
maxsnp <- 1000
example_snps <- c(506, 519, 517) ## OD, bias, and outlier, in that order
snp_nums <- c(order_vec[example_snps], order_vec[(1:maxsnp)[-example_snps]])

example_ref <- t(ref_mat[snp_nums, ])
example_alt <- t(alt_mat[snp_nums, ])
example_tot <- example_ref + example_alt
colnames(example_ref) <- paste0("SNP", 1:maxsnp)
colnames(example_tot) <- paste0("SNP", 1:maxsnp)

write.csv(example_ref, file = "./Output/shirasawa_snps/example_refcounts.csv", row.names = TRUE)
write.csv(example_tot, file = "./Output/shirasawa_snps/example_readcounts.csv", row.names = TRUE)

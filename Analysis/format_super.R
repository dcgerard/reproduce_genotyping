## format data for SuperMASSA

ploidy      <- 6
osize_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)[, 1:3]
ocounts_mat <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)[, 1:3]

for (index in 1:ncol(ocounts_mat)) {
  ocounts <- ocounts_mat[, index]
  osize   <- osize_mat[, index]
  rownames_vec <- rownames(osize_mat)
  rownames_vec <- rownames_vec[!is.na(ocounts)]
  ocounts <- ocounts[!is.na(ocounts)]
  osize   <- osize[!is.na(osize)]
  dat <- data.frame(A = ocounts[-1], a = osize[-1] - ocounts[-1])
  pdat <- data.frame(A = ocounts[1], a = osize[1] - ocounts[1])
  rownames(dat) <- rownames_vec[-1]
  rownames(pdat) <- rownames_vec[1]
  write.table(x = dat, file = paste0("./Output/supermassa_formatted_data/osnp", index, ".txt"),
              col.names = FALSE, quote = FALSE)
  write.table(x = pdat, file = paste0("./Output/supermassa_formatted_data/psnp", index, ".txt"),
              col.names = FALSE, quote = FALSE)
}

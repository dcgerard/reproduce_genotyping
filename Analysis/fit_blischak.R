## Extract sequencing error rates so that Blischak can use them.

ufiles <- list.files("./Output/updog_fits/")
seq_vec <- matrix(NA, nrow = length(ufiles), ncol = 1)
for (index in 1:length(ufiles)) {
  uout <- readRDS(file = paste0("./Output/updog_fits/", ufiles[index]))
  seq_vec[index, 1] <- uout$seq
}
write.table(x = seq_vec, file = "./Output/blischak_formatted_data/seq_error.txt", row.names = FALSE, col.names = FALSE)


counts_mat <- read.csv("./Output/shirasawa_snps/example_refcounts.csv", row.names = 1)
size_mat   <- read.csv("./Output/shirasawa_snps/example_readcounts.csv", row.names = 1)
ploidy <- 6

write.table(x = counts_mat, file = "./Output/blischak_formatted_data/counts_mat.txt",
            row.names = FALSE, col.names = FALSE, na = "-9") ## missing values are coded with a -9
write.table(x = size_mat, file = "./Output/blischak_formatted_data/size_mat.txt",
            row.names = FALSE, col.names = FALSE, na = "-9") ## missing values are coded with a -9

## Now run Blischak code
command_text <- paste0("ebg diseq -p ", ploidy, " -n ", nrow(counts_mat), " -l ", ncol(counts_mat),
                       " -r ./Output/blischak_formatted_data/counts_mat.txt ",
                       "-t ./Output/blischak_formatted_data/size_mat.txt ",
                       "-e ./Output/blischak_formatted_data/seq_error.txt")
btime <- system.time({
  system(command_text)
})
system("mv diseq* ./Output/blischak_formatted_data/")
saveRDS(btime, "./Output/blischak_formatted_data/btime.RDS")

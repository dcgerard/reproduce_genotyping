shirasawa_snps = ./Output/shirasawa_snps/example_refcounts.csv \
		 ./Output/shirasawa_snps/example_readcounts.csv

all : ./Output/fig/od_arg.pdf

$(shirasawa_snps) : ./Data/KDRIsweetpotatoXushu18S1LG2017.vcf.gz
	mkdir -p ./Output/shirasawa_snps
	Rscript ./Analysis/parse_vcf.R

./Output/fig/od_arg.pdf : $(shirasawa_snps) ./Analysis/od_argument.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/od_argument.R


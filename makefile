shirasawa_snps = ./Output/shirasawa_snps/example_refcounts.csv \
		 ./Output/shirasawa_snps/example_readcounts.csv

$(shirasawa_snps) : ./Data/KDRIsweetpotatoXushu18S1LG2017.vcf.gz
	mkdir -p Output/shirasawa_snps
	Rscript ./Analysis/parse_vcf.R


shirasawa_snps = ./Output/shirasawa_snps/example_refcounts.csv \
		 ./Output/shirasawa_snps/example_readcounts.csv

od_output = ./Output/text/od_summaries.txt \
            ./Output/fig/od_arg.pdf

bias_output = ./Output/fig/bias_arg.pdf \
	      ./Output/text/bias_summaries.txt

ufits = ./Output/updog_fits/uout1.RDS \
	./Output/updog_fits/uout2.RDS \
	./Output/updog_fits/uout3.RDS

all : $(od_output) \
      $(bias_output) \
      ./Output/text/out_prob.txt \
      ./Output/fig/prior_quantiles.pdf \
      ./Output/fig/snp_examples.pdf \
      ./Output/fig/prob_plots.pdf \
      ./Output/fig/updog_fits.pdf

$(shirasawa_snps) : ./Data/KDRIsweetpotatoXushu18S1LG2017.vcf.gz ./Analysis/parse_vcf.R
	mkdir -p ./Output/shirasawa_snps
	Rscript ./Analysis/parse_vcf.R

$(od_output) : $(shirasawa_snps) ./Analysis/od_argument.R
	mkdir -p ./Output/fig
	mkdir -p ./Output/text
	Rscript ./Analysis/od_argument.R

$(bias_output) : $(shirasawa_snps) ./Analysis/bias_arg.R
	mkdir -p ./Output/fig
	mkdir -p ./Output/text
	Rscript ./Analysis/bias_arg.R

./Output/text/out_prob.txt : $(ufits) ./Analysis/out_arg.R
	mkdir -p ./Output/text
	Rscript ./Analysis/out_arg.R

./Output/fig/prior_quantiles.pdf : $(shirasawa_snps) ./Analysis/plot_quantiles.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_quantiles.R

./Output/fig/snp_examples.pdf : $(shirasawa_snps) ./Analysis/plot_raw.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_raw.R

./Output/fig/prob_plots.pdf : $(shirasawa_snps) ./Analysis/possible_probs.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/possible_probs.R

$(ufits) : $(shirasawa_snps) ./Analysis/fit_all_updog.R
	mkdir -p ./Output/updog_fits
	Rscript ./Analysis/fit_all_updog.R

./Output/fig/updog_fits.pdf : $(ufits) ./Analysis/plot_updog_fits.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_updog_fits.R

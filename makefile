shirasawa_snps = ./Output/shirasawa_snps/example_refcounts.csv \
		 ./Output/shirasawa_snps/example_readcounts.csv

od_output = ./Output/text/od_summaries.txt \
            ./Output/fig/od_arg.pdf

bias_output = ./Output/fig/bias_arg.pdf \
	      ./Output/text/bias_summaries.txt

ufits = ./Output/updog_fits/uout1.RDS \
	./Output/updog_fits/uout2.RDS \
	./Output/updog_fits/uout3.RDS

blischak_fits = ./Output/blischak_formatted_data/counts_mat.txt \
		./Output/blischak_formatted_data/size_mat.txt \
		./Output/blischak_formatted_data/seq_error.txt \
		./Output/blischak_formatted_data/diseq-F.txt \
		./Output/blischak_formatted_data/diseq-freqs.txt \
		./Output/blischak_formatted_data/diseq-genos.txt

all : $(od_output) \
      $(bias_output) \
      ./Output/text/out_prob.txt \
      ./Output/fig/prior_quantiles.pdf \
      ./Output/fig/snp_examples.pdf \
      ./Output/fig/prob_plots.pdf \
      ./Output/fig/updog_fits.pdf \
      $(blischak_fits)

# Extract the reference and alternative counts from the shirasawa et al data.
$(shirasawa_snps) : ./Data/KDRIsweetpotatoXushu18S1LG2017.vcf.gz ./Analysis/parse_vcf.R
	mkdir -p ./Output/shirasawa_snps
	Rscript ./Analysis/parse_vcf.R

# Fit updog on the top 1000 covered SNPs.
$(ufits) : $(shirasawa_snps) ./Analysis/fit_all_updog.R
	mkdir -p ./Output/updog_fits
	Rscript ./Analysis/fit_all_updog.R

# Argument for overdispersion.
$(od_output) : $(shirasawa_snps) ./Analysis/od_argument.R
	mkdir -p ./Output/fig
	mkdir -p ./Output/text
	Rscript ./Analysis/od_argument.R

# Argument for bias.
$(bias_output) : $(shirasawa_snps) ./Analysis/bias_arg.R
	mkdir -p ./Output/fig
	mkdir -p ./Output/text
	Rscript ./Analysis/bias_arg.R

# Argument for outliers
./Output/text/out_prob.txt : $(ufits) ./Analysis/out_arg.R
	mkdir -p ./Output/text
	Rscript ./Analysis/out_arg.R

# Prior quantiles for bias and sequencing error rates
./Output/fig/prior_quantiles.pdf : $(shirasawa_snps) ./Analysis/plot_quantiles.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_quantiles.R

# Raw plots of the three SNP's I picked out
./Output/fig/snp_examples.pdf : $(shirasawa_snps) ./Analysis/plot_raw.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_raw.R

# Showing the different combinations of sequencing error rate and bias
./Output/fig/prob_plots.pdf : $(shirasawa_snps) ./Analysis/possible_probs.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/possible_probs.R

# Plot the updog fits of the three SNP's I picked out
./Output/fig/updog_fits.pdf : $(ufits) ./Analysis/plot_updog_fits.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_updog_fits.R

$(blischak_fits) : $(ufits) $(shirasawa_snps) ./Analysis/fit_blischak.R
	mkdir -p ./Output/blischak_formatted_data
	Rscript ./Analysis/fit_blischak.R

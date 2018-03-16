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

supermassa_fits = ./Output/supermassa_formatted_data/supermassa_out1.txt \
		  ./Output/supermassa_formatted_data/supermassa_out2.txt \
		  ./Output/supermassa_formatted_data/supermassa_out3.txt

sim_plots = ./Output/fig/allele_freq_est.pdf \
	    ./Output/fig/bias_v_af.pdf \
            ./Output/fig/od_v_af.pdf \
            ./Output/fig/od_v_bias.pdf \
            ./Output/fig/param_ests.pdf \
            ./Output/fig/prop_correct_lines.pdf \
            ./Output/fig/seq_v_af.pdf

all : sweet_potato simulations

.PHONY : sweet_potato
sweet_potato : $(od_output) \
      $(bias_output) \
      ./Output/text/out_prob.txt \
      ./Output/fig/prior_quantiles.pdf \
      ./Output/fig/snp_examples.pdf \
      ./Output/fig/prob_plots.pdf \
      ./Output/fig/updog_fits.pdf \
      $(blischak_fits) \
      ./Output/fig/blischak_fits.pdf \
      ./Output/fig/ident_prob.pdf \
      ./Output/fig/supermassa_fits.pdf \
      ./Output/fig/real_data_plots.pdf \
      ./Output/fig/ufit_features.pdf \
      ./Output/fig/seq_error_example.pdf \
      ./Output/fig/comp_time.pdf

.PHONY : simulations
simulations : ./Output/sims_out/sims_out.csv \
	      $(sim_plots) \
	      ./Output/hwe_f1_sims/hwe_f1_sims_out.csv \
	      ./Output/fig/hwe_s1_prop_correct.pdf


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

# Argument for sequencing error rate
./Output/fig/seq_error_example.pdf : $(shirasawa_snps) ./Analysis/seq_arg.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/seq_arg.R

# Showing the different combinations of sequencing error rate and bias
./Output/fig/prob_plots.pdf : $(shirasawa_snps) ./Analysis/possible_probs.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/possible_probs.R

# Plot the updog fits of the three SNP's I picked out
./Output/fig/updog_fits.pdf : $(ufits) ./Analysis/plot_updog_fits.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_updog_fits.R

# Plot computation time of all 1000 updog fits
./Output/fig/comp_time.pdf : $(ufits) ./Analysis/computation.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/computation.R

# Run ebg on the data using sequencing error rates estimated from updog
$(blischak_fits) : $(ufits) $(shirasawa_snps) ./Analysis/fit_blischak.R
	mkdir -p ./Output/blischak_formatted_data
	Rscript ./Analysis/fit_blischak.R

# plot blischak output
./Output/fig/blischak_fits.pdf : $(blischak_fits) ./Analysis/plot_blischak.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_blischak.R

# Hypothetical data with identifiability issue
./Output/fig/ident_prob.pdf : $(ufits) ./Analysis/hypothetical_problem.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/hypothetical_problem.R

# Plot SuperMASSA fits from http://statgen.esalq.usp.br/SuperMASSA/
./Output/fig/supermassa_fits.pdf : $(supermassa_fits) ./Analysis/plot_supermassa.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_supermassa.R

# Plot Supermassa, updog, and Blischak on the SNPs I chose
./Output/fig/real_data_plots.pdf : $(blischak_fits) $(ufits) $(supermassa_fits) ./Analysis/plot_combined.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_combined.R

# Run Simulations where simulate and fit under HWE.
./Output/sims_out/sims_out.csv : ./Output/shirasawa_snps/example_readcounts.csv ./Analysis/run_sims.R
	mkdir -p ./Output/sims_out
	mkdir -p ./Output/blischak_formatted_data_sims/
	Rscript ./Analysis/run_sims.R

# Plot simulations where simulate and fit under HWE
$(sim_plots) : ./Output/sims_out/sims_out.csv ./Analysis/plot_sims.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/plot_sims.R

# Run Simulations where simulate under S1 and fit under both S1 and HWE
./Output/hwe_f1_sims_out.csv : ./Output/shirasawa_snps/example_readcounts.csv ./Analysis/hwe_f1_sims.R
	mkdir -p ./Output/hwe_f1_sims
	Rscript ./Analysis/hwe_f1_sims.R

# Plot Simulations where simulate under S1 and fit under both S1 and HWE
./Output/fig/hwe_s1_prop_correct.pdf : ./Output/hwe_f1_sims/hwe_f1_sims_out.csv ./Analysis/hwe_f1_plots.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/hwe_f1_plots.R

# Histograms of updog estimates on Shirasawa data
./Output/fig/ufit_features.pdf : $(ufits) ./Analysis/shirasawa_data_features.R
	mkdir -p ./Output/fig
	Rscript ./Analysis/shirasawa_data_features.R

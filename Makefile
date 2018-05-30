# R scripting front-end. Note that makeCluster sometimes fails to
# connect to a socker when using Rscript, so we are using the "R CMD
# BATCH" interface instead.
rexec = R CMD BATCH --no-save --no-restore

# AVOID EDITING ANYTHING BELOW THIS LINE
# --------------------------------------

# These variables contain csv files for readcounts from the
# Shirasawa dataset.
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

# Figures output from plot_run_sims.R
sim_plots = ./Output/fig/allele_freq_est.pdf \
      ./Output/fig/diff.pdf \
      ./Output/fig/param_ests.pdf \
      ./Output/fig/bias_v_af.pdf \
      ./Output/fig/od_v_af.pdf \
      ./Output/fig/prop_correct.pdf \
      ./Output/fig/corr.pdf \
      ./Output/fig/od_v_bias.pdf \
      ./Output/fig/seq_v_af.pdf

# Figures output from plot_run_pp_sims.R
pp_sim_plots = ./Output/fig/est_weight_pp.pdf \
      ./Output/fig/prop_correct_pp.pdf

# Figures for oracle estimator
oracle_dir = ./Output/oracle/fig
oracle_figs = $(oracle_dir)/cor_n_for_05_alpha50_ploidy2.pdf \
	$(oracle_dir)/cor_n_for_05_alpha50_ploidy4.pdf \
	$(oracle_dir)/cor_n_for_05_alpha50_ploidy6.pdf \
	$(oracle_dir)/cor_n_for_05_alpha60_ploidy2.pdf \
	$(oracle_dir)/cor_n_for_05_alpha60_ploidy4.pdf \
	$(oracle_dir)/cor_n_for_05_alpha60_ploidy6.pdf \
	$(oracle_dir)/cor_n_for_05_alpha70_ploidy2.pdf \
	$(oracle_dir)/cor_n_for_05_alpha70_ploidy4.pdf \
	$(oracle_dir)/cor_n_for_05_alpha70_ploidy6.pdf \
	$(oracle_dir)/cor_n_for_05_alpha80_ploidy2.pdf \
	$(oracle_dir)/cor_n_for_05_alpha80_ploidy4.pdf \
	$(oracle_dir)/cor_n_for_05_alpha80_ploidy6.pdf \
	$(oracle_dir)/cor_n_for_05_alpha90_ploidy2.pdf \
	$(oracle_dir)/cor_n_for_05_alpha90_ploidy4.pdf \
	$(oracle_dir)/cor_n_for_05_alpha90_ploidy6.pdf \
	$(oracle_dir)/cor_n_for_05_alpha95_ploidy2.pdf \
	$(oracle_dir)/cor_n_for_05_alpha95_ploidy4.pdf \
	$(oracle_dir)/cor_n_for_05_alpha95_ploidy6.pdf \
	$(oracle_dir)/n_for_05_alpha50_ploidy2.pdf \
	$(oracle_dir)/n_for_05_alpha50_ploidy4.pdf \
	$(oracle_dir)/n_for_05_alpha50_ploidy6.pdf \
	$(oracle_dir)/n_for_05_alpha60_ploidy2.pdf \
	$(oracle_dir)/n_for_05_alpha60_ploidy4.pdf \
	$(oracle_dir)/n_for_05_alpha60_ploidy6.pdf \
	$(oracle_dir)/n_for_05_alpha70_ploidy2.pdf \
	$(oracle_dir)/n_for_05_alpha70_ploidy4.pdf \
	$(oracle_dir)/n_for_05_alpha70_ploidy6.pdf \
	$(oracle_dir)/n_for_05_alpha80_ploidy2.pdf \
	$(oracle_dir)/n_for_05_alpha80_ploidy4.pdf \
	$(oracle_dir)/n_for_05_alpha80_ploidy6.pdf \
	$(oracle_dir)/n_for_05_alpha90_ploidy2.pdf \
	$(oracle_dir)/n_for_05_alpha90_ploidy4.pdf \
	$(oracle_dir)/n_for_05_alpha90_ploidy6.pdf \
	$(oracle_dir)/n_for_05_alpha95_ploidy2.pdf \
	$(oracle_dir)/n_for_05_alpha95_ploidy4.pdf \
	$(oracle_dir)/n_for_05_alpha95_ploidy6.pdf

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
	      $(pp_sim_plots) \
	      $(oracle_figs) \
	      ./Output/hwe_f1_sims/hwe_f1_sims_out.csv \
	      ./Output/fig/hwe_s1_prop_correct.pdf \
	      ./Output/sims_out/sims_out_pp.csv


# Extract the reference and alternative counts from the shirasawa et al data.
$(shirasawa_snps) : ./Analysis/parse_vcf.R ./Data/KDRIsweetpotatoXushu18S1LG2017.vcf.gz 
	mkdir -p ./Output/shirasawa_snps
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Fit updog on the top 1000 covered SNPs.
$(ufits) : ./Analysis/fit_all_updog.R $(shirasawa_snps) 
	mkdir -p ./Output/updog_fits
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Argument for overdispersion.
$(od_output) : ./Analysis/od_argument.R $(shirasawa_snps) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/text
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Argument for bias.
$(bias_output) : ./Analysis/bias_arg.R $(shirasawa_snps) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/text
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Argument for outliers
./Output/text/out_prob.txt : ./Analysis/out_arg.R $(ufits) 
	mkdir -p ./Output/text
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Prior quantiles for bias and sequencing error rates
./Output/fig/prior_quantiles.pdf : ./Analysis/plot_quantiles.R $(shirasawa_snps) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Raw plots of the three SNP's I picked out
./Output/fig/snp_examples.pdf : ./Analysis/plot_raw.R $(shirasawa_snps) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Argument for sequencing error rate
./Output/fig/seq_error_example.pdf : ./Analysis/seq_arg.R $(shirasawa_snps) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Showing the different combinations of sequencing error rate and bias
./Output/fig/prob_plots.pdf : ./Analysis/possible_probs.R $(shirasawa_snps) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Plot the updog fits of the three SNP's I picked out
./Output/fig/updog_fits.pdf : ./Analysis/plot_updog_fits.R $(ufits) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Plot computation time of all 1000 updog fits
./Output/fig/comp_time.pdf : ./Analysis/computation.R $(ufits) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Run ebg on the data using sequencing error rates estimated from updog
$(blischak_fits) : ./Analysis/fit_blischak.R $(ufits) $(shirasawa_snps) 
	mkdir -p ./Output/blischak_formatted_data
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# plot blischak output
./Output/fig/blischak_fits.pdf : ./Analysis/plot_blischak.R $(blischak_fits) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Hypothetical data with identifiability issue
./Output/fig/ident_prob.pdf : ./Analysis/hypothetical_problem.R $(ufits) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Plot SuperMASSA fits from http://statgen.esalq.usp.br/SuperMASSA/
./Output/fig/supermassa_fits.pdf : ./Analysis/plot_supermassa.R $(supermassa_fits) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Plot Supermassa, updog, and Blischak on the SNPs I chose
./Output/fig/real_data_plots.pdf : ./Analysis/plot_combined.R $(blischak_fits) $(ufits) $(supermassa_fits) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Run Simulations where simulate and fit under HWE.
./Output/sims_out/sims_out.csv : ./Analysis/run_sims.R $(shirasawa_snps)
	mkdir -p ./Output/sims_out
	mkdir -p ./Output/blischak_formatted_data_sims/
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Run Simulations where simulate under preferential pairing
./Output/sims_out/sims_out_pp.csv : ./Analysis/run_pp_sims.R $(shirasawa_snps)
	mkdir -p ./Output/sims_out
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Plot simulations where simulate and fit under HWE
$(sim_plots) : ./Analysis/plot_run_sims.R ./Output/sims_out/sims_out.csv 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Plot simulations where simulate under preferential pairing
$(pp_sim_plots) : ./Analysis/plot_run_pp_sims.R ./Output/sims_out/sims_out_pp.csv 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Run Simulations where simulate under S1 and fit under both S1 and HWE
./Output/hwe_f1_sims_out.csv : ./Analysis/hwe_f1_sims.R $(shirasawa_snps)
	mkdir -p ./Output/hwe_f1_sims
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Plot Simulations where simulate under S1 and fit under both S1 and HWE
./Output/fig/hwe_s1_prop_correct.pdf : ./Analysis/hwe_f1_plots.R ./Output/hwe_f1_sims/hwe_f1_sims_out.csv 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Histograms of updog estimates on Shirasawa data
./Output/fig/ufit_features.pdf : ./Analysis/shirasawa_data_features.R $(ufits) 
	mkdir -p ./Output/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Joint distributions for oracle estimator and true genotype
./Output/oracle/odat.RDS : ./Analysis/oracle_explore.R
	mkdir -p ./Output/oracle
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

# Create figures from output of oracle_explore.R
$(oracle_figs) : ./Analysis/plot_oracle.R ./Output/oracle/odat.RDS
	mkdir -p ./Output/oracle
	mkdir -p ./Output/oracle/fig
	mkdir -p ./Output/rout
	$(rexec) $< Output/rout/$(basename $(notdir $<)).Rout

Reproducing Results of Gerard, Ferrao, and Stephens (2017)
================

Introduction
============

This repository contains code to reproduce the empirical evaluations of Gerard, Ferrão, and Stephens (2017). The new methods can be found in the [updog](https://github.com/dcgerard/updog) package.

If you are having trouble reproducing these results, it might be that you need to update some of your R packages. These are the versions that I used:

``` r
sessionInfo()
```

    ## R version 3.3.2 (2016-10-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.2 LTS
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] rmutil_1.1.0    updog_0.1.0     dplyr_0.5.0     purrr_0.2.2    
    ## [5] readr_1.0.0     tidyr_0.6.1     tibble_1.2      ggplot2_2.2.1  
    ## [9] tidyverse_1.1.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.11     plyr_1.8.4       forcats_0.2.0    tools_3.3.2     
    ##  [5] digest_0.6.12    jsonlite_1.3     lubridate_1.6.0  evaluate_0.10   
    ##  [9] nlme_3.1-131     gtable_0.2.0     lattice_0.20-34  psych_1.6.12    
    ## [13] DBI_0.6          yaml_2.1.14      parallel_3.3.2   haven_1.0.0     
    ## [17] xml2_1.1.1       stringr_1.2.0    httr_1.2.1       knitr_1.15.1    
    ## [21] hms_0.3          rprojroot_1.2    grid_3.3.2       R6_2.2.0        
    ## [25] readxl_0.1.1     foreign_0.8-67   rmarkdown_1.3    modelr_0.1.0    
    ## [29] reshape2_1.4.2   magrittr_1.5     backports_1.0.5  scales_0.4.1    
    ## [33] htmltools_0.3.5  rvest_0.3.2      assertthat_0.2.0 mnormt_1.5-5    
    ## [37] colorspace_1.3-2 stringi_1.1.2    lazyeval_0.2.0   munsell_0.4.3   
    ## [41] broom_0.4.2

As you can see above, I've also only tried this out on Ubuntu.

If you find a bug, please create an [issue](https://github.com/dcgerard/reproduce_genotyping/issues).

Instructions
============

To reproduce the results in Gerard, Ferrão, and Stephens (2017), you need to (1) download the appropriate R packages, (2) obtain the appropriate data, (3) run `make`, and (4) get coffee.

Install R Packages
------------------

To install the needed R packages, run the following in R

``` r
install.packages("tidyverse")
devtools::install_github("dcgerard/updog")
```

Get Data
--------

Place the following files in the Data folder:

1.  [KDRIsweetpotatoXushu18S1LG2017.vcf.gz](http://sweetpotato-garden.kazusa.or.jp/)

The actual URL is: <ftp://ftp.kazusa.or.jp/pub/sweetpotato/GeneticMap/>

Run Make
--------

To reproduce all of the results in (<span class="citeproc-not-found" data-reference-id="gerard2016unifying">**???**</span>), simply run `make` from the terminal.

Get Coffee
----------

You should get some coffee. Here is a list of some of my favorite places:

-   Chicago
    -   [Sawada Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
    -   [Plein Air Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
-   Seattle
    -   [Bauhaus Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
    -   [Cafe Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
-   Columbus
    -   [Yeah, Me Too](https://www.yelp.com/biz/yeah-me-too-columbus)
    -   [Stauf's Coffee Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)
    -   [Caffe Apropos](https://www.yelp.com/biz/caff%C3%A9-apropos-columbus-2)

References
==========

Gerard, David, Luis Felipe Ventorim Ferrão, and Matthew Stephens. 2017. “Harnessing Mendelian Segregation and Empirical Bayes for Genotyping Autopolyploids.” *Overleaf Preprint*.

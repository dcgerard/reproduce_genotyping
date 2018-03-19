Reproducing Results of Gerard et al. (2018)
================

[![DOI](https://zenodo.org/badge/94825101.svg)](https://zenodo.org/badge/latestdoi/94825101)

Introduction
============

This repository contains code to reproduce the empirical evaluations of Gerard et al. (2018). The new methods can be found in the [updog](https://github.com/dcgerard/updog) package.

If you are having trouble reproducing these results, it might be that you need to update some of your R packages. These are the versions that I used:

``` r
sessionInfo()
```

    ## R version 3.4.3 (2017-11-30)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 15063)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3   ggthemes_3.4.0  rmutil_1.1.0    updog_0.1.0    
    ##  [5] forcats_0.2.0   stringr_1.3.0   dplyr_0.7.4     purrr_0.2.4    
    ##  [9] readr_1.1.1     tidyr_0.7.2     tibble_1.4.2    ggplot2_2.2.1  
    ## [13] tidyverse_1.2.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] reshape2_1.4.3   haven_1.1.0      lattice_0.20-35  colorspace_1.3-2
    ##  [5] htmltools_0.3.6  yaml_2.1.18      rlang_0.2.0      pillar_1.2.1    
    ##  [9] foreign_0.8-69   glue_1.2.0       modelr_0.1.1     readxl_1.0.0    
    ## [13] bindrcpp_0.2     bindr_0.1        plyr_1.8.4       munsell_0.4.3   
    ## [17] gtable_0.2.0     cellranger_1.1.0 rvest_0.3.2      psych_1.7.8     
    ## [21] evaluate_0.10.1  knitr_1.20       broom_0.4.3      Rcpp_0.12.16    
    ## [25] scales_0.5.0     backports_1.1.2  jsonlite_1.5     mnormt_1.5-5    
    ## [29] hms_0.4.0        digest_0.6.15    stringi_1.1.6    grid_3.4.3      
    ## [33] rprojroot_1.3-2  cli_1.0.0        tools_3.4.3      magrittr_1.5    
    ## [37] lazyeval_0.2.1   crayon_1.3.4     pkgconfig_2.0.1  xml2_1.2.0      
    ## [41] lubridate_1.7.1  assertthat_0.2.0 rmarkdown_1.9    httr_1.3.1      
    ## [45] rstudioapi_0.7   R6_2.2.2         nlme_3.1-131     compiler_3.4.3

I've also only tried this out on Ubuntu.

If you find a bug, please create an [issue](https://github.com/dcgerard/reproduce_genotyping/issues).

Instructions
============

To reproduce the results in Gerard et al. (2018), you need to (1) download and install the appropriate packages, (2) obtain the appropriate data, (3) run `make`, and (4) get coffee.

Install Packages
----------------

To install the needed R packages, run the following in R

``` r
install.packages(c("tidyverse", "rmutil", "snow", "parallel", 
                   "ggthemes", "gridExtra", "devtools"))
devtools::install_github("dcgerard/updog")
```

Please follow the directions [here](https://github.com/pblischak/polyploid-genotyping/tree/master/ebg) to install `ebg`.

Get Data
--------

Place [KDRIsweetpotatoXushu18S1LG2017.vcf.gz](http://sweetpotato-garden.kazusa.or.jp/) in the Data folder. The direct url is: <ftp://ftp.kazusa.or.jp/pub/sweetpotato/GeneticMap/>

Run Make
--------

To reproduce all of the results in Gerard et al. (2018), simply run `make` from the terminal. To reproduce the real-data analysis, run

``` bash
make sweet_potato
```

To reproduce the simulations, run

``` bash
make simulations
```

Get Coffee
----------

The simulations should take a few hours. You should get some coffee. Here is a list of some of my favorite places:

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

Note on SuperMASSA
==================

The [source](https://bitbucket.org/orserang/supermassa) is available for SuperMASSA, but it isn't too well documented, so I used the [web application](http://statgen.esalq.usp.br/SuperMASSA/) version during the empirical data analysis. As such, I've committed the inputs and fits in the [Output](https://github.com/dcgerard/reproduce_genotyping/tree/master/Output/supermassa_formatted_data) folder. If you want to generate these fits, you'll have to do it manually.

References
==========

Gerard, David, Luis Felipe Ventorim Ferrão, Antonio Augusto Franco Garcia, and Matthew Stephens. 2018. “Harnessing Empirical Bayes and Mendelian Segregation for Genotyping Autopolyploids from Messy Sequencing Data.” *bioRxiv*. Cold Spring Harbor Laboratory. doi:[10.1101/281550](https://doi.org/10.1101/281550).

Reproducing Results of Gerard et al. (2018)
================

[![DOI](https://zenodo.org/badge/94825101.svg)](https://zenodo.org/badge/latestdoi/94825101)

Introduction
============

This repository contains code to reproduce the empirical evaluations of Gerard et al. (2018). The new methods can be found in the [updog](https://cran.r-project.org/package=updog) package on CRAN.

If you are having trouble reproducing these results, it might be that you need to update some of your R packages. These are the versions that I used:

``` r
sessionInfo()
```

    ## R version 3.5.0 (2018-04-23)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 16.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/local/lib/R/lib/libRblas.so
    ## LAPACK: /usr/local/lib/R/lib/libRlapack.so
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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] fitPoly_3.0.0    updogAlpha_1.0.1 gridExtra_2.3    ggthemes_4.0.0  
    ##  [5] rmutil_1.1.1     updog_1.0.1      forcats_0.3.0    stringr_1.3.1   
    ##  [9] dplyr_0.7.4      purrr_0.2.5      readr_1.1.1      tidyr_0.8.0     
    ## [13] tibble_1.4.2     ggplot2_3.0.0    tidyverse_1.2.1 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] reshape2_1.4.3            haven_1.1.1              
    ##  [3] lattice_0.20-35           colorspace_1.3-2         
    ##  [5] htmltools_0.3.6           yaml_2.1.19              
    ##  [7] rlang_0.2.1               pillar_1.2.2             
    ##  [9] foreign_0.8-70            glue_1.2.0               
    ## [11] withr_2.1.2               modelr_0.1.2             
    ## [13] readxl_1.1.0              bindrcpp_0.2.2           
    ## [15] foreach_1.4.4             bindr_0.1.1              
    ## [17] plyr_1.8.4                munsell_0.4.3            
    ## [19] gtable_0.2.0              cellranger_1.1.0         
    ## [21] rvest_0.3.2               codetools_0.2-15         
    ## [23] psych_1.8.4               evaluate_0.10.1          
    ## [25] knitr_1.20                RcppArmadillo_0.8.600.0.0
    ## [27] doParallel_1.0.11         broom_0.4.4              
    ## [29] Rcpp_0.12.17              scales_0.5.0             
    ## [31] backports_1.1.2           jsonlite_1.5             
    ## [33] mnormt_1.5-5              hms_0.4.2                
    ## [35] digest_0.6.15             stringi_1.2.2            
    ## [37] grid_3.5.0                rprojroot_1.3-2          
    ## [39] cli_1.0.0                 tools_3.5.0              
    ## [41] magrittr_1.5              lazyeval_0.2.1           
    ## [43] crayon_1.3.4              pkgconfig_2.0.1          
    ## [45] xml2_1.2.0                lubridate_1.7.4          
    ## [47] iterators_1.0.9           assertthat_0.2.0         
    ## [49] rmarkdown_1.9             httr_1.3.1               
    ## [51] rstudioapi_0.7            R6_2.2.2                 
    ## [53] nlme_3.1-137              compiler_3.5.0

I've also only tried this out on Ubuntu.

If you find a bug, please create an [issue](https://github.com/dcgerard/reproduce_genotyping/issues).

Instructions
============

To reproduce the results in Gerard et al. (2018), you need to (1) download and install the appropriate packages, (2) obtain the appropriate data, (3) run `make`, and (4) get coffee.

Install Packages
----------------

To install the needed R packages, run the following in R

``` r
install.packages(c("updog", "tidyverse", "rmutil", "snow", "parallel", 
                   "ggthemes", "gridExtra", "devtools", "fitPoly"))
devtools::install_github("dcgerard/updogAlpha")
```

If the most recent CRAN version of `updog` doesn't seem to work, you can retry with the version that I last used to reproduce these results:

``` r
devtools::install_github("dcgerard/updog", 
                         ref = "76c72eb717e18061576fc20b3c07f9da71b67263")
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

-   Washington, D.C.
    -   [Politics and Prose](https://www.yelp.com/biz/politics-and-prose-washington)
    -   [Colony Club](https://www.yelp.com/biz/colony-club-washington)
-   Chicago
    -   [Sawada Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
    -   [Plein Air Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
-   Seattle
    -   [Bauhaus Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
    -   [Cafe Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
-   Columbus
    -   [Yeah, Me Too](https://www.yelp.com/biz/yeah-me-too-columbus)
    -   [Fox in the Snow](https://www.yelp.com/biz/fox-in-the-snow-cafe-columbus-2)
    -   [Stauf's Coffee Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)
    -   [Caffe Apropos](https://www.yelp.com/biz/caff%C3%A9-apropos-columbus-2)

Note on SuperMASSA
==================

The [source](https://bitbucket.org/orserang/supermassa) is available for SuperMASSA, but it isn't too well documented, so I used the [web application](http://statgen.esalq.usp.br/SuperMASSA/) version during the empirical data analysis. As such, I've committed the inputs and fits in the [Output](https://github.com/dcgerard/reproduce_genotyping/tree/master/Output/supermassa_formatted_data) folder. If you want to generate these fits, you'll have to do it manually.

References
==========

Gerard, David, Luis Felipe Ventorim Ferrão, Antonio Augusto Franco Garcia, and Matthew Stephens. 2018. “Harnessing Empirical Bayes and Mendelian Segregation for Genotyping Autopolyploids from Messy Sequencing Data.” *bioRxiv*. Cold Spring Harbor Laboratory. doi:[10.1101/281550](https://doi.org/10.1101/281550).

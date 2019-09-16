Reproducing Results of Gerard et al. (2018)
================

[![DOI](https://zenodo.org/badge/94825101.svg)](https://zenodo.org/badge/latestdoi/94825101)

# Introduction

This repository contains code to reproduce the empirical evaluations of
Gerard et al. (2018). The new methods can be found in the
[updog](https://cran.r-project.org/package=updog) package on CRAN.

If you are having trouble reproducing these results, it might be that
you need to update some of your R packages. These are the versions that
I used:

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
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
    ##  [1] fitPoly_3.0.0    updogAlpha_1.0.1 gridExtra_2.3    ggthemes_4.2.0  
    ##  [5] updog_1.1.1      forcats_0.4.0    stringr_1.4.0    dplyr_0.8.3     
    ##  [9] purrr_0.3.2      readr_1.3.1      tidyr_1.0.0      tibble_2.1.3    
    ## [13] ggplot2_3.2.1    tidyverse_1.2.1 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.5          xfun_0.9                 
    ##  [3] haven_2.1.1               lattice_0.20-38          
    ##  [5] colorspace_1.4-1          vctrs_0.2.0              
    ##  [7] generics_0.0.2            htmltools_0.3.6          
    ##  [9] yaml_2.2.0                rlang_0.4.0              
    ## [11] pillar_1.4.2              glue_1.3.1               
    ## [13] withr_2.1.2               modelr_0.1.5             
    ## [15] readxl_1.3.1              foreach_1.4.7            
    ## [17] lifecycle_0.1.0           munsell_0.5.0            
    ## [19] gtable_0.3.0              cellranger_1.1.0         
    ## [21] rvest_0.3.4               codetools_0.2-16         
    ## [23] evaluate_0.14             knitr_1.24               
    ## [25] RcppArmadillo_0.9.700.2.0 doParallel_1.0.15        
    ## [27] broom_0.5.2               Rcpp_1.0.2               
    ## [29] scales_1.0.0              backports_1.1.4          
    ## [31] jsonlite_1.6              hms_0.5.1                
    ## [33] digest_0.6.20             stringi_1.4.3            
    ## [35] grid_3.6.1                cli_1.1.0                
    ## [37] tools_3.6.1               magrittr_1.5             
    ## [39] lazyeval_0.2.2            crayon_1.3.4             
    ## [41] pkgconfig_2.0.2           zeallot_0.1.0            
    ## [43] xml2_1.2.2                lubridate_1.7.4.9000     
    ## [45] iterators_1.0.12          assertthat_0.2.1         
    ## [47] rmarkdown_1.15            httr_1.4.1               
    ## [49] rstudioapi_0.10           R6_2.4.0                 
    ## [51] nlme_3.1-141              compiler_3.6.1

I’ve also only tried this out on Ubuntu.

If you find a bug, please create an
[issue](https://github.com/dcgerard/reproduce_genotyping/issues).

# Instructions

To reproduce the results in Gerard et al. (2018), you need to (1)
download and install the appropriate packages, (2) obtain the
appropriate data, (3) run `make`, and (4) get coffee.

## Install Packages

To install the needed R packages, run the following in R

``` r
install.packages(c("updog", "tidyverse", "rmutil", "snow", "parallel", 
                   "ggthemes", "gridExtra", "devtools", "fitPoly"))
devtools::install_github("dcgerard/updogAlpha")
```

If the most recent CRAN version of `updog` doesn’t seem to work, you can
retry with the version that I last used to reproduce these results:

``` r
devtools::install_github("dcgerard/updog", 
                         ref = "76c72eb717e18061576fc20b3c07f9da71b67263")
```

Please follow the directions
[here](https://github.com/pblischak/polyploid-genotyping/tree/master/ebg)
to install `ebg`.

## Get Data

Place
[KDRIsweetpotatoXushu18S1LG2017.vcf.gz](http://sweetpotato-garden.kazusa.or.jp/)
in the Data folder. The direct url is:
<ftp://ftp.kazusa.or.jp/pub/sweetpotato/GeneticMap/>

## Run Make

To reproduce all of the results in Gerard et al. (2018), simply run
`make` from the terminal. To reproduce the real-data analysis, run

``` bash
make sweet_potato
```

To reproduce the simulations, run

``` bash
make simulations
```

## Get Coffee

The simulations should take a few hours. You should get some coffee.
Here is a list of some of my favorite places:

  - Washington, D.C.
      - [Politics and
        Prose](https://www.yelp.com/biz/politics-and-prose-washington)
      - [Colony Club](https://www.yelp.com/biz/colony-club-washington)
  - Chicago
      - [Sawada Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
      - [Plein Air
        Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
  - Seattle
      - [Bauhaus
        Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
      - [Cafe Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
  - Columbus
      - [Yeah, Me Too](https://www.yelp.com/biz/yeah-me-too-columbus)
      - [Fox in the
        Snow](https://www.yelp.com/biz/fox-in-the-snow-cafe-columbus-2)
      - [Stauf’s Coffee
        Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)
      - [Caffe
        Apropos](https://www.yelp.com/biz/caff%C3%A9-apropos-columbus-2)

# Note on SuperMASSA

The [source](https://bitbucket.org/orserang/supermassa) is available for
SuperMASSA, but it isn’t too well documented, so I used the [web
application](http://statgen.esalq.usp.br/SuperMASSA/) version during the
empirical data analysis. As such, I’ve committed the inputs and fits in
the
[Output](https://github.com/dcgerard/reproduce_genotyping/tree/master/Output/supermassa_formatted_data)
folder. If you want to generate these fits, you’ll have to do it
manually.

# References

<div id="refs" class="references">

<div id="ref-gerard2018genotyping">

Gerard, David, Luís Felipe Ventorim Ferrão, Antonio Augusto Franco
Garcia, and Matthew Stephens. 2018. “Genotyping Polyploids from Messy
Sequencing Data.” *Genetics* 210 (3). Genetics: 789–807.
<https://doi.org/10.1534/genetics.118.301468>.

</div>

</div>

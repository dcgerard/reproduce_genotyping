Reproducing Results of Gerard et al. (2018)
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

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
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
    ##  [1] gridExtra_2.3   ggthemes_4.2.4  rmutil_1.1.9    updog_2.1.3    
    ##  [5] forcats_0.5.1   stringr_1.4.0   dplyr_1.0.9     purrr_0.3.4    
    ##  [9] readr_2.1.2     tidyr_1.2.0     tibble_3.1.7    ggplot2_3.3.6  
    ## [13] tidyverse_1.3.1
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.8.3             lubridate_1.8.0          listenv_0.8.0           
    ##  [4] assertthat_0.2.1         digest_0.6.29            foreach_1.5.2           
    ##  [7] utf8_1.2.2               parallelly_1.31.1        R6_2.5.1                
    ## [10] cellranger_1.1.0         backports_1.4.1          reprex_2.0.1            
    ## [13] evaluate_0.15            httr_1.4.3               pillar_1.7.0            
    ## [16] rlang_1.0.2              readxl_1.4.0             rstudioapi_0.13         
    ## [19] rmarkdown_2.14           munsell_0.5.0            broom_0.8.0             
    ## [22] compiler_4.2.0           modelr_0.1.8             xfun_0.31               
    ## [25] pkgconfig_2.0.3          globals_0.15.0           htmltools_0.5.2         
    ## [28] tidyselect_1.1.2         codetools_0.2-18         doFuture_0.12.2         
    ## [31] fansi_1.0.3              future_1.26.1            crayon_1.5.1            
    ## [34] tzdb_0.3.0               dbplyr_2.1.1             withr_2.5.0             
    ## [37] grid_4.2.0               jsonlite_1.8.0           gtable_0.3.0            
    ## [40] lifecycle_1.0.1          DBI_1.1.2                magrittr_2.0.3          
    ## [43] scales_1.2.0             cli_3.3.0                stringi_1.7.6           
    ## [46] doRNG_1.8.2              RcppArmadillo_0.11.1.1.0 fs_1.5.2                
    ## [49] xml2_1.3.3               ellipsis_0.3.2           generics_0.1.2          
    ## [52] vctrs_0.4.1              iterators_1.0.14         tools_4.2.0             
    ## [55] glue_1.6.2               rngtools_1.5.2           hms_1.1.1               
    ## [58] fastmap_1.1.0            yaml_2.3.5               colorspace_2.0-3        
    ## [61] rvest_1.0.2              knitr_1.39               haven_2.5.0

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
[KDRIsweetpotatoXushu18S1LG2017.vcf.gz](https://github.com/dcgerard/KDRIsweetpotatoXushu18S1LG2017)
in the Data folder. The direct url is:
<https://github.com/dcgerard/KDRIsweetpotatoXushu18S1LG2017/raw/main/KDRIsweetpotatoXushu18S1LG2017.vcf.gz>

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

-   Washington, D.C.
    -   [Politics and
        Prose](https://www.yelp.com/biz/politics-and-prose-washington)
    -   [Colony Club](https://www.yelp.com/biz/colony-club-washington)
-   Chicago
    -   [Sawada Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
    -   [Plein Air
        Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
-   Seattle
    -   [Bauhaus
        Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
    -   [Cafe Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
-   Columbus
    -   [Yeah, Me Too](https://www.yelp.com/biz/yeah-me-too-columbus)
    -   [Fox in the
        Snow](https://www.yelp.com/biz/fox-in-the-snow-cafe-columbus-2)
    -   [Stauf’s Coffee
        Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)
    -   [Caffe
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

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-gerard2018genotyping" class="csl-entry">

Gerard, David, Luís Felipe Ventorim Ferrão, Antonio Augusto Franco
Garcia, and Matthew Stephens. 2018. “Genotyping Polyploids from Messy
Sequencing Data.” *Genetics* 210 (3): 789–807.
<https://doi.org/10.1534/genetics.118.301468>.

</div>

</div>

#!/bin/bash -e

R --slave -e 'install.packages(c("cowplot", "tidyverse"), repos="https://cran.rstudio.com/")' && \
R --slave -e 'if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager", repos="https://cran.rstudio.com/")' && \
R --slave -e 'BiocManager::install("dada2")'

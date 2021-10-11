BootStrap: docker 
From: continuumio/miniconda3
#From: centos:7
#OSVersion: 7
#MirrorURL: http://mirror.centos.org/centos-%{OSVERSION}/%{OSVERSION}/os/$basearch/
#Include: yum wget


%setup


%environment


%post
	apt-get install build-essential
	
	conda install -y -c conda-forge mamba
	mamba install 'r-base>4'
	# conda install -y -c conda-forge 'r-base>4'

	R --slave -e 'install.packages(c("cowplot", "tidyverse"), repos="https://cran.rstudio.com/")' && \
	R --slave -e 'if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager", repos="https://cran.rstudio.com/")' && \
	R --slave -e 'BiocManager::install("dada2")'	

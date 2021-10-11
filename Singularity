BootStrap: docker 
From: continuumio/miniconda3
#From: centos:7
#OSVersion: 7
#MirrorURL: http://mirror.centos.org/centos-%{OSVERSION}/%{OSVERSION}/os/$basearch/
#Include: yum wget


%setup


%environment


%post
	rm -vf /var/lib/apt/lists/*
	apt-get update -y
	apt-get install build-essential libreadline-dev libncurses5 r-base r-base-core r-recommended r-base-dev -y
	cd /lib/x86_64-linux-gnu/
	ln -s libreadline.so.7.0 libreadline.so.6



	
	#conda install -y -c conda-forge mamba
	#mamba install 'r-base'
	# conda install -y -c conda-forge 'r-base>4'

	R --version

	R --slave -e 'install.packages(c("cowplot", "tidyverse"), repos="https://cran.rstudio.com/")' && \
	R --slave -e 'if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager", repos="https://cran.rstudio.com/")' && \
	R --slave -e 'BiocManager::install("dada2")'	

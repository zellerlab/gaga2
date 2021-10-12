Bootstrap: docker
OSVersion: ubuntu:20.04
# MirrorURL: http://archive.ubuntu.com/ubuntu/


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

	#R --slave -e 'install.packages(c("devtools"))
	#R --slave -e 'library("devtools"); devtools::install_github("benjjneb/dada2", ref="v1.20"
	#	

	R --slave -e 'install.packages(c("devtools", "cowplot", "tidyverse"), repos="https://cran.rstudio.com/")' && \
	R --slave -e 'if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager", repos="https://cran.rstudio.com/")' && \
	R --slave -e 'BiocManager::install("dada2", version = "3.11")'	

BootStrap: docker 
From: centos:7
OSVersion: 7
MirrorURL: http://mirror.centos.org/centos-%{OSVERSION}/%{OSVERSION}/os/$basearch/
Include: yum wget


%setup


%environment


%post


    mkdir -p /opt/software

    ### Install your packages ###

    # update yum
    yum makecache fast && \
    yum update -y

    yum -y install git bzip2 wget which sudo vi source zlib-devel xz-devel bzip2-devel
    yum -y group install "Development Tools"


	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/miniconda
    # echo "PATH=/opt/miniconda/bin:\$PATH" >> /root/.bashrc
    # echo "export PATH" >> /root/.bashrc
    rm Miniconda3-latest-Linux-x86_64.sh
    # source /root/.bashrc

	conda install -y -c conda-forge 'r-base>4'

	R --slave -e 'install.packages(c("cowplot", "tidyverse"), repos="https://cran.rstudio.com/")' && \
	R --slave -e 'if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager", repos="https://cran.rstudio.com/")' && \
	R --slave -e 'BiocManager::install("dada2")'	

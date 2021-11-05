Bootstrap: docker
From: ubuntu:18.04
IncludeCmd: yes

%environment
R_VERSION=4.1
export R_VERSION
R_CONFIG_DIR=/etc/R/
export R_CONFIG_DIR
export LC_ALL=C
export PATH=$PATH

%post
  apt-get update
  apt-get install -y apt-transport-https apt-utils software-properties-common
  apt-get install -y add-apt-key
  export DEBIAN_FRONTEND=noninteractive
  ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
  apt-get install -y tzdata
  dpkg-reconfigure --frontend noninteractive tzdata

  #add CRAN/Ubuntu repo, add key, then refresh
  apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
  add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
  apt-get update

  apt-get install -y wget python3-pip git nano libblas3 libblas-dev liblapack-dev liblapack3 curl gcc fort77 aptitude g++ xorg-dev libreadline-dev gfortran libssl-dev libxml2-dev libpcre3-dev liblzma-dev libbz2-dev libcurl4-openssl-dev libhdf5-dev hdf5-helpers libmariadb-client-lgpl-dev
#  apt-get install -y wget nano
#  apt-get install -y libblas3 libblas-dev liblapack-dev liblapack3 curl
#  apt-get install -y gcc fort77 aptitude
#  aptitude install -y g++
#  aptitude install -y xorg-dev
#  aptitude install -y libreadline-dev
#  aptitude install -y gfortran
#  gfortran --version
#  apt-get install -y libssl-dev libxml2-dev libpcre3-dev liblzma-dev libbz2-dev libcurl4-openssl-dev 
#  apt-get install -y libhdf5-dev hdf5-helpers libmariadb-client-lgpl-dev
#
  

  git clone https://github.com/cschu/figaro.git  # Zymo-Research repo does not yet allow 0 length primers
  cd figaro
  git checkout fix/primer_length_check
  python3 setup.py bdist_wheel
  pip3 install --force-reinstall dist/*.whl

  apt-get install -y r-base r-base-dev
  R --version
  
  # installing packages from cran
  R --slave -e 'install.packages(c("devtools", "tidyverse", "cowplot"),repos="https://cran.rstudio.com/")'

  # installing from bioc
  R --slave -e 'if (!requireNamespace("BiocManager",quietly=TRUE)) install.packages("BiocManager")'
  R --slave -e 'BiocManager::install(version = "3.13", ask = FALSE)'
  R --slave -e 'BiocManager::install(c("dada2"))'

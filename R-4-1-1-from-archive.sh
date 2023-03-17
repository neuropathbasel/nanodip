#!/bin/bash
# installation of R
# J. Hench, IfP 2020,2021,2022
# ubuntu 18.04 apt-based installations, work out-of-the-box (outcomment once installated)

## install packages required for python and GPUMAP
#sudo apt-get install -y libopenblas-dev libomp-dev
#sudo apt-get -y install python3.7 python3.7-venv python3.7-dev
#sudo apt-get -y install ncbi-blast+
#sudo apt-get -y install r-base
## R compilation dependencies (https://github.com/Jiefei-Wang/Painless-R-compilation-and-installation-on-Ubuntu)
# gcc-multilib not avaiable for AGX # sudo apt-get -y install gcc-multilib
#sudo apt-get -y install build-essential fort77 xorg-dev liblzma-dev libblas-dev gfortran gobjc++ aptitude libreadline-dev libbz2-dev libpcre2-dev libcurl4 libcurl4-openssl-dev default-jre default-jdk openjdk-8-jdk openjdk-8-jre  texinfo texlive texlive-fonts-extra -y libssl-dev -y libxml2-dev 
## GenoGAM dependency
# sudo apt-get -y install libjpeg-dev

venvpath=/applications/
myR="R-4.1.1"
myRbranch="R-4"
#startDir=`pwd`
cd $venvpath
#myRstamp=$myR"_"`date +%s`
#mkdir -p $myRstamp 
#cd $myRstamp
# local R installation
wget `echo "https://cran.microsoft.com/snapshot/2022-01-01/src/base/"$myRbranch"/"$myR".tar.gz"`
tar -xzvf `echo "./"$myR".tar.gz"`
cd $myR
./configure
make
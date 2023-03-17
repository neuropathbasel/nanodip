#!/bin/bash
# Setup of NanoDiP venv with Python 3.7.x
# J. Hench, IfP Basel, 2021
# jupyter notebook and minknow_api won't work togehter in python 3.8.8
# Ubuntu 18.04 jetpack requirements:
sudo apt install python3.7 python3.7-dev python3.7-venv python3-pip


venvname=nanodipVenv01
venvpath=/applications
mypython=python3.7
startDir=`pwd`
mkdir -p $venvpath/$venvname

cd $venvpath
$mypython -m venv $venvname
source $venvpath/$venvname/bin/activate


# essential components
pip install --upgrade pip
pip install wheel
pip install jupyter

# convenience components for jupyter notebook
pip install jupyter_contrib_nbextensions
jupyter contrib nbextension install --user
pip install jupyter_nbextensions_configurator
jupyter nbextensions_configurator enable --user

# remaining components in alphabetical order
pip install cupy
pip install cherrypy
pip install grpcio # required for MinKNOW API, replaces old "pip install grpc" which does not work any longer (2021-10-11)
pip install -U kaleido
pip install matplotlib
pip install numpy
pip install numba
pip install openpyxl
pip install pandas
pip install plotly
pip install psutil
pip install pysam
pip install reportlab==3.6.1 # 20220308 discovered bug in current reportlab version, hence the version pinning
pip install seaborn
pip install scikit-learn
pip install neuralnetwork
pip install tqdm
pip install umap-learn
pip install xhtml2pdf==0.2.5 # 20220308 discovered bug in current reportlab version, hence the version pinning

# install a patched version of the minkow API:
# AND NOT! pip install minknow_api #install from github https://github.com/neuropathbasel/minknow_api
cd $venvpath/$venvname
git clone https://github.com/neuropathbasel/minknow_api
cd minknow_api
python setup.py install

# launch the jupyter notebook
# cd /applications/nanodip
# python -m notebook

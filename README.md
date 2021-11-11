# NanoDiP
Nanopore Digital Pathology - NanoDiP, the Nanopore companion of EpiDiP (http://www.epidip.org)

NanoDiP operates Nanopore (ONT) MinION sequencers through the MinKNOW API and performs UMAP-based dimension reduction of the obtained methylation patterns against large numbers of reference datasets derived from the Illumina Infinium 450K/EPIC platform. It also generates copy number plots from aligned reads. Analysis is performed in a live manner during the sequencing. It runs on Ubuntu 18.04 Linux computers with 32GB RAM, 8 cores CPU, and a CUDA GPU with >500 cores.

*This is work in progress*. The software within this repository is intended for development purpose only. Please contact the developer if you need more information. Currently supported platforms are the NVIDIA Jetson AGX Xavier and the x64_64 + NVIDIA CUDA GPU PCs running Ubuntu 18.04. Other platforms are not supported and are not planned for support in the near future.

![NanoDiP_v20_20211014](https://user-images.githubusercontent.com/59837805/137954721-ff179777-ac40-4e14-898f-c23ffea3d8a5.gif)

# Disclaimer
Implementation of this software in a diagnostic settings occurs in the sole responsibility of the treating physician. Usage of this software occurs at the risk of the user. The authors may not be held liable for any damage (including hardware) this software might cause.

# Installation
NanoDiP requires a distinct setup of software / hardware. Instructions on how to install NanoDiP are currently being prepared. Installation scripts that (should) compile all required dependencies are available at https://github.com/neuropathbasel/nanodip_dependencies. 

The intended location for the installation is:

`/applications/nanodip` for this repository
`/applications/nanodip_dependencies` for the installation routines, which place compiled packages in `/applications/`.

The default and preferably *only* user name for the setup is `minit` on an ARMv8 platform and `minknow` on an x86_64 computer. Data acquisition and post-hoc data analysis should be possible without GPUs, too. While MinION sequencers seem to work with USB2, they are intended to communicate through USB3. USB-C>>USB3 adaptors work well. SSD drives (e.g., 1TB) are strongly recommended for the OS and immediate data storage locations. Such modules (e.g., M2 keys) can be connected to the AGX through PCIe adapters; do not use eSATA or USB3 interfaces (on the AGX). We use conventional (large) USB3 hard drives for post-hoc local data archives.

# Usage
At present, there is no user manual. The code can either be executed as jupyter notebook, or as python script. The python script in this repository was gerated from the jupyter notebook version by downloading it "as .py". In order to understand how things work under the hood, please examine the jupyter notebook (also while executing the code). The script launches a CherryPy web server listetning on port 8080 of "localhost". There is no login prevention, and we advise to only access the web-based UI locally with chromium-browser. Firefox works except for WebGL rendering on some hardware, which negatively impacts on interactive plotly plots. In order to run, a python virtual environment is suggested. This venv can be generated with the respective script in this repository.

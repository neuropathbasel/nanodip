# NanoDiP
Nanopore Digital Pathology - NanoDiP, the Nanopore companion of EpiDiP (http://www.epidip.org)

NanoDiP operates Nanopore (ONT) MinION sequencers through the MinKNOW API and performs UMAP-based dimension reduction of the obtained methylation patterns against large numbers of reference datasets derived from the Illumina Infinium 450K/EPIC platform. It also generates copy number plots from aligned reads. Analysis is performed in a live manner during the sequencing. It runs on Ubuntu 18.04 Linux computers with 32GB RAM, 8 cores CPU, and a GPU with >500 cores.

This is work in progress. The software within this repository is intended for development purpose only. Please contact the developer if you need more information. Currently supported platforms are the NVIDIA Jetson AGX Xavier and the x64_64 / NVIDIA CUDA GPU Platforms running Ubuntu 18.04. Other platforms are not supported and are not planned for support in the near future.

![NanoDiP_v20_20211014](https://user-images.githubusercontent.com/59837805/137954721-ff179777-ac40-4e14-898f-c23ffea3d8a5.gif)

# Disclaimer
Implementation of this software in a diagnostic settings occurs in the sole responsibility of the treating physician.

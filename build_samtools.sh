#!/bin/bash

sudo apt install libhts-dev libbz2-dev libncurses-dev


# download sources from https://sourceforge.net/projects/samtools/files/samtools/1.12/ for htslib and samtools
# unzip to /applications/htslib, /applications/samtools


cd /applications/htslib
./configure    # Optional but recommended, for choosing extra functionality
make

cd /applications/samtools
./configure
make




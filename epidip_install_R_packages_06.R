# testing with R 4.1.1 - works on Jetson AGX 2022-05-2
# J. Hench, IfP Basel 2021-2022

# sudo apt -y install libgit2-dev (required for devtools on AGX)

# CRAN mirror to use:
local({r <- getOption("repos")
       r["CRAN"] <- "https://cran.microsoft.com/snapshot/2022-01-01/" 
       options(repos=r)
})

# install CRAN packages
install.packages(c(
 "doParallel",
 "plyr",
 "devtools",
 "xlsx",
 "locfit",
 "BiocManager"))

library(BiocManager)
repositories(version="3.14",site_repository="https://bioconductor.statistik.tu-dortmund.de/packages/3.14/")
BiocManager::install(version = "3.14", ask=F, update=F, pkgs=c(
 "minfi",
 "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
 "IlluminaHumanMethylationEPICanno.ilm10b3.hg19",
 "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
 "IlluminaHumanMethylation450kanno.ilmn12.hg19",
 "IlluminaHumanMethylation27kanno.ilmn12.hg19",
 "IlluminaHumanMethylationEPICmanifest",
 "IlluminaHumanMethylation450kmanifest",
 "IlluminaHumanMethylation27kmanifest",
 "conumee",
 "GWASTools",
 "minfiData"
))

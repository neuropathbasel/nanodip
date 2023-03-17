# J. Hench, 2019-2022
# IfP Basel
# get filtered betas for single core
# 450K and/or EPIC Illumina Methylation arrays; only 450K-equivalent probes will be considered

getSentrix <- function(fullPath){
  a<-unlist(strsplit(fullPath,"/"))
  b<-a[length(a)]
  return(b)
}
message("program start: ",Sys.time())
# ------ LIBRARIES --------
library(doParallel)
library(plyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# ------ INPUT PARAMETERS --------
baseDir <-                "/applications/epidip_demo_data/data/demo_idat/" # baseDir <- "/imagesets/IDAT_450K_EPIC_symlinks/" #baseDir <- "/mnt/optane01/imagesets/TCGA_GSE90496_IFP_symlinks/" #baseDir <- "/mnt/optane01/imagesets/IFP_TCGA_mix_mini_broken/" #" # e.g., baseDir <- "/mnt/optane01/imagesets/IFP_idat/20190211/"
outputDir <-              "/applications/epidip_demo_data/data/demo_output/"
betaLogFile <-            "/applications/tmp/beta.log"
commonProbesCsv <-        "/applications/reference_data/infinium_conumee/featuresCommon.csv"
crossReactiveProbesCsv <- "/applications/reference_data/infinium_conumee/48639-non-specific-probes-Illumina450k.csv"
binIndex<-paste0(outputDir,"index.csv")
suffix <- "_Grn.idat" # to be clipped to create minfi-compatible "base pathes"
cpuCap <- 0.5 #execution cap (present_CPUs*cpuCap)
# ---- NO INPUT PARAMETERS BELOW --------
if (file.exists(betaLogFile)==T){
	file.remove(betaLogFile)
}
allFiles <- list.files(path=baseDir, recursive=TRUE, full.names=TRUE, include.dirs=TRUE)
idatFileDf <- data.frame(allFiles[grepl("_Grn.idat",allFiles)]) # collect list of IDAT pathes, linux shell version, should be re-coded in pure R to work without the command line interpreter
colnames(idatFileDf) <- c("basePath") # rename the idatFileDf colum to address it later on
idatFileDf$basePath <- gsub(suffix,"",idatFileDf$basePath) # remove the suffix (see INPUT PARAMETERS) since minfi requires base pathes
idatFileDf<-mdply(idatFileDf, function(basePath,outputDir){
  a<-unlist(strsplit(basePath,"/"))
  b<-a[length(a)]
  return(b)
})
colnames(idatFileDf) <- c("basePath","sentrixId") # rename the idatFileDf colum to address it later on
message("process IDAT to beta values (single): ",Sys.time())
cl<-makeCluster(round(detectCores()*cpuCap),outfile=betaLogFile) # create a cluster for parallel processing
registerDoParallel(cl) # parallel processing
iterations <- nrow(idatFileDf) # number of IDATs to process
message("number of IDATs to process: ",iterations)
startNo <- 1
endNo <-   iterations
keep1 <- (read.csv(file=commonProbesCsv, stringsAsFactors=FALSE)$x) # define filters: 1) keep features present in both EPIC and 450K data sets
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)# if your data includes males and females, remove probes on the sex chromosomes, read in 450k annotation data
keep2 <- keep1[!(keep1 %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])] # 2) exclude sex chromosomes
xReactiveProbes <- read.csv(file=crossReactiveProbesCsv, stringsAsFactors=FALSE)  # 3) exclude cross-reactive probes (based on a publication by Chen et al. (2013))
keep3 <- keep2[!(keep2 %in% xReactiveProbes$TargetID)] # 4) merge into "keep3 Large character" which contains a list of all probe names
foreach(i=startNo:endNo, .errorhandling = 'remove') %dopar%{ # need to reload libraries within parallel threads
	library(minfi)
	library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
	library(IlluminaHumanMethylation450kmanifest)
  betasBinFileName <- paste0(outputDir,idatFileDf$sentrixId[i],"_betas_filtered.bin")
  reanalyze<-F
  if (file.exists(betasBinFileName)){
    message("betas file exists, checking size... ")
    s<-file.info(betasBinFileName)$size
    if(s>3300000){ # if larger than 4 MB (will ignore 27K etc.)
      message("betas calculated already",idatFileDf$basePath[i])
      system (paste0("echo betas calculated already ",i,"/",iterations," ",idatFileDf$basePath[i]," ",Sys.time()))
    }else{
      message("betas file faulty ",betasBinFileName[i])
      system (paste0("echo betas file faulty ",betasBinFileName[i],Sys.time()))
      reanalyze<-T
    }
  }else{
    reanalyze<-T
  }
  if(reanalyze==T){
    rgSet <- read.metharray(basenames=idatFileDf$basePath[i],extended=FALSE,verbose=TRUE,force=TRUE)
    mSetSqFlt0 <- convertArray(mapToGenome(preprocessSWAN(rgSet)),"IlluminaHumanMethylation450k") # rgSet conversion seems to result in different numbers of elements
    keepFilter <- (featureNames(mSetSqFlt0) %in% keep3) # create logical filter to apply to methylSet based on probes actually contained in methylSet
    mSetSqFlt1 <- mSetSqFlt0[keepFilter,]
    mSetSqFlt2 <- dropLociWithSnps(mSetSqFlt1) # remove probes with SNPs at CpG site
    system (paste0("echo ",i,"/",iterations," ",Sys.time())) # works on linux to display status of individual threads
    singleBetas <- getBeta(mSetSqFlt2) # calculate beta values (return value)
    betas<-singleBetas[,1]
    betas[is.na(betas)]<-0.49 # replace NAs with 0.49
    writeBin(betas,betasBinFileName)
    if(file.exists(binIndex)==F){ # only once:
      writeLines(rownames(singleBetas),binIndex)
    }
  }
}
stopCluster(cl)
message("process IDAT to beta values (completed): ",Sys.time())
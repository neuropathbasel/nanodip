# generate copy number plots from Infinium 27K, 450K, and EPIC IDAT Methylome sets
# code modified from samples kindly provided by Martin Sill / Anne Schöler / David Capper 2018-5-25
# Modifications for Institut für Medizinische Genetik und Pathologie, USB, Basel, by Jürgen Hench, 2019-2022
# including extension to 27K, 450K, and EPIC datasets

# LIBRARIES for single core-----------------------------
library(doParallel)
# END LIBRARIES -----------------------------
# FUNCTIONS -----------------------------
getSentrix <- function(fullPath){
  a<-unlist(strsplit(fullPath,"/"))
  b<-a[length(a)]
  return(b)
}
# END FUNCTIONS ---------------------------
# INPUT PARAMETERS ---------------------------------------------------------------------------------------------------
searchDir <- "/applications/epidip_demo_data/data/demo_idat/"
outputDir <- "/applications/epidip_demo_data/data/demo_output/"
minIdatSize <- 7500000 # a typical 450K size is 8095252 bytes
yPlotRange <- 1.2
oneCopy <- 0.5
# INPUT PARAMETERS: if this is used as a module, these could be set as defaults
conumeeWorkingDirectory <-    "/applications/tmp"
conumeeRefEpicRdaPath <-      "/applications/reference_data/infinium_conumee/CNanalysis5_conumee_REF.2017-02-10.RData" # load refEPIC.data, provided by Martin Sill
conumeeAnnoEpicRdaPath <-     "/applications/reference_data/infinium_conumee/IlluminaArrayDBconumee_annotation_EPIC_B4.2017-06-07.RData" # load annoEPIC and annoEPICxy, provided by Martin Sill
conumeeEpicManifestRdaPath <- "/applications/reference_data/infinium_conumee/epicManifest.rda" #pre-loaded EPIC manifest CSV as dataframe below line 7
customAnnotationsCsv <-       "/applications/reference_data/infinium_conumee/annotationGenes_27K.csv"
ref450Kfile <-                "/applications/reference_data/infinium_conumee/refIDAT/ref450K_Red.idat"
cnvLogFile <-                 "/applications/tmp/cnv.log"
cpuCap <- 0.5 #execution cap (present_CPUs*cpuCap)
cpuPlot <- round(detectCores()*cpuCap)
cpuPlot <- 2 # override; RAM limitation!
# NO FURTHER MODIFYABLE / INPUT PARAMETERS BELOW ---------------------------------------------------------------------
grnSuffix <- "_Grn.idat" # collect all IDAT file locations within search directory based on the green files
idatGrn <- list.files(path=searchDir, pattern = paste0("*",grnSuffix), recursive = TRUE)
idatPathes <- ""
idatPathes <- data.frame(idatGrn)
idatPathes$idatPath <- with(idatPathes, paste0(searchDir,gsub(grnSuffix,"",idatPathes$idatGrn))) # cut away with suffix from IDAT path, leaving only Sentrix ID and plate position
message(">>> collected all files to plot ",Sys.time())
numFiles <- nrow(idatPathes)
message(">>> number of files: ",numFiles)
firstrun<-TRUE # will be set FALSE after initialzation for each core (or once in a normal loop)
if (cpuPlot > numFiles){
  cpuPlot <- numFiles
}
message(">>> launching ",cpuPlot," parallel tasks for plotting ",Sys.time())
cl<-makeCluster(cpuPlot,outfile=cnvLogFile) # select number of cores
registerDoParallel(cl) # parallel processing
#for (idatCounter in 1:length(idatPathes$idatPath)){
#for (idatCounter in 2:4){
foreach(idatCounter=1:length(idatPathes$idatPath),.errorhandling = 'remove') %dopar%{
  library(conumee)
  library(devtools)
  library(DNAcopy)
  #library(GenoGAM)
  library(GenomicRanges)
  library(GWASTools)
  library(IlluminaHumanMethylationEPICmanifest)
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylation27kmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b3.hg19)
  library(minfi)
  library(minfiData)
  if (firstrun==TRUE){   # during first run (per core, load reference data to harmonize annotations to reference IDAT)
    message(">>> loading reference IDAT to create custom annotation ",Sys.time())
    setwd(conumeeWorkingDirectory)
    data(centromeres.hg19) # load centromere data (for primitive plot)
    load(conumeeRefEpicRdaPath) # load Epic control set and annotation, provided by Martin Sill
    load(conumeeAnnoEpicRdaPath)
    message(">>> EPIC control set and annotation loaded ",Sys.time())
    targetAnnotations <- read.csv(customAnnotationsCsv,header=TRUE) # create annotations from csv
    detail<- GRanges(targetAnnotations$targetChromosome, ranges=IRanges(targetAnnotations$targetCoordinateStart,targetAnnotations$targetCoordinateEnd), strand=targetAnnotations$targetStrand, name=targetAnnotations$targetGene) #taken from ncbi
    genome(detail) <- "hg19"
    myAnno <- CNV.create_anno(array_type = "450k", chrXY = TRUE, detail_regions=detail) # create new annotation object; will be modified by smallest possible array, i.e. 27K default bin_minprobes = 15, min_bin_size = 50000, max_bin_size = 5e+06
    message(">>> custom annotation created from CSV ",Sys.time())
    load(conumeeEpicManifestRdaPath) # load preprocessed EPIC manifest file
    message(">>> EPIC manifest RDA loaded ",Sys.time())
    myRawIntensitySet <- read.metharray(ref450Kfile, force=TRUE) # load a 27K IDAT to adjust annotations to this level: read in raw intensity signals
    myMethylSet <- convertArray(preprocessIllumina(myRawIntensitySet),"IlluminaHumanMethylationEPIC") # Create Methyl Set
    message(">>> 450K reference: generated methylation set from raw intensity values ",Sys.time())
    myCustomMethylSet <- mapToGenome(myMethylSet) # find overlaps between custom annotation and Mset
    myAnno@probes <- subsetByOverlaps(myAnno@probes, granges(myCustomMethylSet)) # myAnno will be used for all plots
    message(">>> 450K reference: determined overlaps between custom annotation and Mset ",Sys.time())
    firstrun<-FALSE 
  } # END of firstrun code
  epicPath <- idatPathes$idatPath[idatCounter] # start processing sample
  sampleID <- getSentrix(epicPath)
  outputPath <- outputDir # the shared path prefix for all output files, if they shall be stored in the same place as the input file. Example for epicPath: "/mnt/mydrive/myepicfolder/mysamplefolder/202135260076_R06C01"
  cnvPdfPath <- paste0(outputPath,sampleID,"_CNV_IFPBasel_annotations.pdf")
  cnvRdsPath <- paste0(outputPath,sampleID,"_CNV_IFPBasel_annotations.rds")
  cnvBinPath <- paste0(outputPath,sampleID,"_CNV.bin")
  cnvCsvPath <- paste0(outputPath,sampleID,"_CNV.csv")
  cnvCsvGzPath <- paste0(outputPath,sampleID,"_CNV.csv.gz")
  cnvFitPath <- paste0(outputPath,sampleID,"_CNV_fit.csv")
  cnvFitGzPath <- paste0(outputPath,sampleID,"_CNV_fit.csv.gz")
  cnvFitBinPath <- paste0(outputPath,sampleID,"_CNV_fit.bin")
  cnvSegPath <- paste0(outputPath,sampleID,"_CNV_seg.csv")
  cnvSegGzPath <- paste0(outputPath,sampleID,"_CNV_seg.csv.gz")
  cnvIndexPath <- paste0(outputPath,"CNV_index.csv")
  message (paste0("### processing ",sampleID," in ",epicPath," ",Sys.time()))
  if(file.exists(cnvFitGzPath)==F){
		if (file.info(paste0(epicPath,grnSuffix))$size>minIdatSize){
		  myRawIntensitySet <- read.metharray(epicPath, force=TRUE)   # read in raw intensity signals
		  message(">>> raw intensity signals loaded ",Sys.time())
		  myMethylSet <- convertArray(preprocessIllumina(myRawIntensitySet),"IlluminaHumanMethylationEPIC")   # Create Methyl Set
		  message(">>> generated methylation set from raw intensity values ",Sys.time())
		  myCustomMethylSet <- mapToGenome(myMethylSet)   # find overlaps between custom annotation and Mset
		  message(">>> determined overlaps between custom annotation and Mset ",Sys.time())
		  myCNV.data <- CNV.load(myMethylSet) # load CNV data
		  mysampleCust.fit <- CNV.fit(myCNV.data[sampleID],refEPIC.data, myAnno)   # find CNVs (custom annotation)
		  mysampleCust.bin <- CNV.bin(mysampleCust.fit)
		  mysampleCust.CNVs <- CNV.segment(CNV.detail(mysampleCust.bin))
		  message(">>> determined CNVs ",Sys.time())
		  message(">>> plotting classical conumee annotations",Sys.time())
		  pdf(cnvPdfPath, height = 8, width = 17)
		  CNV.genomeplot(mysampleCust.CNVs)
		  dev.off()
		  myProbes <- CNV.write(mysampleCust.bin, what='bin') # obtain CNV values for each probe
		  colnames(myProbes)<-c("Chr","Start","End","LocID","binRatio") #rename columns for CSV output
		  myProbes<-data.frame(myProbes$Chr,myProbes$Start,myProbes$End,myProbes$binRatio) # reduce to essential columns
		  saveRDS(mysampleCust.CNVs, file = cnvRdsPath)
    	#writeBin(myProbes$binRatio,cnvBinPath)
    	ratios<-mysampleCust.CNVs@fit$ratio # prepare probe-level data to be saved as CSV file
    	#colnames(ratios)<-c("IlmnId","cnvRatio")
    	write.csv(ratios,file=cnvFitPath) # write CSV
    	system(paste("gzip -f",cnvFitPath)) #gzip-compress CSV
    	writeBin(ratios,cnvFitBinPath) # write binary (float) file
    	seg<-mysampleCust.CNVs@seg$p # prepare segment data to be saved as CSV file
    	segshort<-data.frame(seg$chrom,seg$loc.start,seg$loc.end,seg$seg.mean)
    	colnames(segshort)<-c("chrom","loc_start","loc_end","seg_mean")
    	write.csv(segshort,file=cnvSegPath,row.names=FALSE) # write CSV
    	system(paste("gzip -f",cnvSegPath)) #gzip-compress CSV
    	if(file.exists(cnvIndexPath)==F){ # only once:
      	writeLines(names(ratios),cnvIndexPath)
    	}
		  message(">>> done with sample ",Sys.time())
		}else{
		  message(">>> sample is not a 450K or EPIC idat (size too small) ",Sys.time())
		}
	}else{
		message(paste0(">>> file exists: ",cnvPdfPath,", skipping"),Sys.time())
	}
} # end of main for loop that goes through IDAT sets
stopCluster(cl)
message("done plotting ",Sys.time())
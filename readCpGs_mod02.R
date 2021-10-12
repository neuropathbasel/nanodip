# call overlap methylation from f5c (or nanopolish) output
# nanodip variant (adapted from Euskirchen et al. 2017)
# tested with R 4.0.3
# J. Hench, IfP Basel, 2021

# test if there are 4 additional command line argments, if not, then terminate
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
	message("usage: Rscript [THISSCRIPT.R] [methylation frequency f5c/nanopolish output tsv] [cgmap tsv] [output file methoverlap.tsv]")
} else {
	# read input parameters
	methylationcalls <- args[1]
	ilmncgmapfile <- args[2]
	overlapMethylationFile <- args[3]
	overlapMethylationCountFile <- args[4]

	message("methylationcalls            :", methylationcalls)
	message("ilmncgmapfile               :", ilmncgmapfile)
	message("overlapMethylationFile      :", overlapMethylationFile)
	message("overlapMethylationCountFile :",overlapMethylationCountFile)

	# read data and reference files
	ilmncgmap <- read.table(ilmncgmapfile, sep = "\t", header = F)
	colnames(ilmncgmap) <- c("ilmnid","chromosome","strand","start")
	methfreq <- read.table(methylationcalls, sep="\t", header=T)
	overlapcpgs <- merge(methfreq, ilmncgmap, by=c("chromosome","start"))

	# determine called singleton CpG sites and write out
	CpGcalls <- subset(overlapcpgs, num_cpgs_in_group == 1 & !(chromosome %in% c("chrX","chrY")))
	CpGcalls <- CpGcalls[!duplicated(CpGcalls$ilmnid),]
	CpGcalls$isMethylated <- (CpGcalls$methylated_frequency > 0.5) * 1
	rownames(CpGcalls) <- CpGcalls$ilmnid
	case <- data.frame(t(CpGcalls))
	case <- t(case["isMethylated",])

	write.table(case,file=overlapMethylationFile,col.names=FALSE,row.names=TRUE,quote=FALSE,sep="\t")
	writeLines(as.character(nrow(case)),con=overlapMethylationCountFile)
	message(nrow(case)," overlap CpGs found; done!")
}

#!/bin/bash
# loop through all folders of a MinION sequencing run that has been live-basecalled and sorted into barcode folders
# determine which barcode has most fast5 files (and assume that this is the target barcode)
# J. Hench, IfP Basel, 2021

# input and configuration parameters
sampleid=$1
dataroot=$2
resultroot=/data/nanodip_output

# prepare pathes
inputpath=`echo $dataroot/$sampleid`

# dependencies
outputpath="/data/nanodip_output"
refgenomefa="/applications/reference_data/minimap_data/hg19.fa"
refgenomemmi="/applications/reference_data/minimap_data/hg19_20201203.mmi"
ilmncgmapfile="/applications/reference_data/minimap_data/hg19_HumanMethylation450_15017482_v1-2_cgmap.tsv"

export PATH="/applications/f5c:${PATH}"
export PATH="/applications/nanopolish/minimap2/:${PATH}"
export PATH="/applications/samtools/:${PATH}"
export PATH="/applications/R-4.0.3/bin/:${PATH}"

# functions
gettargetbarcode () { # determine target barcode (needs to be adapted according to sequencing kit, here: RAD-SQK004)
	barcodes="barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12 unclassified"
	targetbarcode="unclassified" # default
	targetcount="0"
	for a in $barcodes; do
		barcodedirlist=`find $inputpath -name $a`
		for b in $barcodedirlist; do
			barcount=`find $b/ -name *.fast5 | wc -l`
			if [ "$barcount" -gt "$targetcount" ]; then
				targetbarcode=$a
				targetcount=$barcount
			fi 
		done
	done
	if [ "2" -gt "$targetcount" ]; then # barcode can not be determined yet; max. 1 fast5 file does not tell anything
		echo "Target barcode can not yet be determined (no more than 1 fast5 file max. per directory)."
	else
		echo "Target barcode determined as "$targetbarcode"."
	fi
}

loopthroughfast5s () {
	targetfast5=`find $inputpath -name *.fast5 | grep $targetbarcode`
	for f in $targetfast5; do 
		fast5filename=`echo $f | sed -e 's/\//\n/g' | tail -n1`    # get the last part of the path, i.e. the fast5 filename
		fast5fileid=`echo $fast5filename | sed -e 's/\.fast5//g'`  # get the file ID without the suffix
		fastqfile=`echo $f | sed -e 's/fast5/fastq/g'`		   # generate the expected MinKNOW path of the corresponding fastq file (from internal guppy)
		echo "f: "$f 
		echo "fast5file: "$fast5filename
		echo "fast5id: "$fast5fileid
		echo "fastqfile: "$fastqfile
		getoverlapmethylation $f $fast5filename $fast5fileid $fastqfile
	done
}

getoverlapmethylation () {
	f=$1; fast5filename=$2; fast5fileid=$3; fastqfile=$4
	if [ -f "$fastqfile" ]; then
		echo "$fastqfile exists."
		if [ -f "$resultroot/$sampleid/$fast5fileid/$fast5fileid-methoverlapcount.txt" ]; then
			echo "CpG output for $fast5fileid exists, skipping creation"
		else
			echo "creating symlinks for $fast5filename and associated fastq"
			mkdir -p $resultroot/$sampleid/$fast5fileid
			ln -fs $f $resultroot/$sampleid/$fast5fileid/
			ln -fs $fastqfile $resultroot/$sampleid/$fast5fileid/
			runf5cetc $resultroot/$sampleid/$fast5fileid $fastqfile $fast5fileid
		fi
	fi
}

runf5cetc () {
	fast5datadir=$1; fastqfile=$2; fast5fileid=$3
	echo "analyzing fast5 and fastq to determine overlap CpGs for "$datapath
	starttime=`date +%s`


	#index, call methylation and get methylation frequencies
	f5c index -t 1 --iop 100 -d $fast5datadir $fastqfile
	# get sorted BAM (4 threads)
	minimap2 -a -x map-ont $refgenomemmi $fastqfile -t 4 | samtools sort -T tmp -o $fast5datadir/$fast5fileid-reads_sorted.bam
	samtools index $fast5datadir/$fast5fileid-reads_sorted.bam

	# set B to 2 megabases (GPU) and 0.4 kreads
	f5c call-methylation -B2000000 -K400 -b $fast5datadir/$fast5fileid-reads_sorted.bam -g $refgenomefa -r $fastqfile  > $fast5datadir/$fast5fileid-result.tsv
	f5c meth-freq -c 2.5 -s -i $fast5datadir/$fast5fileid-result.tsv > $fast5datadir/$fast5fileid-freq.tsv

	# detect overlap CpG sites and their methylation states
	methylationcalls=$fast5datadir/$fast5fileid-freq.tsv
	overlapMethylationFile=$fast5datadir/$fast5fileid-methoverlap.tsv
	overlapMethylationCountFile=$fast5datadir/$fast5fileid-methoverlapcount.txt

	Rscript readCpGs_mod02.R $methylationcalls $ilmncgmapfile $overlapMethylationFile $overlapMethylationCountFile

	endtime=`date +%s`
	echo "------------------------------"
	echo "execution time "`expr $endtime - $starttime`"s"
	echo 
}

# main
mkdir -p $resultroot/$sampleid
gettargetbarcode
loopthroughfast5s


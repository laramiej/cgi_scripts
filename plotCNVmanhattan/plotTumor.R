#!/usr/bin/Rscript

#Make sure the CNVmanhattan.R file is in the same directory
source('CNVmanhattan.R')

#find and store level (ploidy) information from the TumorDetails file
#this is a manual process of inspecting the MEAN_LEVEL fields in the header
#first load rows containing these fields (modify the nrows= depending on the number of levels)
ploidy<-read.table("cnvTumorDetailsBeta-XXX-ASM.tsv", header=T, sep="\t", na.string="N", nrows=10, comment.char="")
#change the number of rows to remove at the begining to only capture the the MEAN_LEVEL fields
ploidy<-ploidy[-(1:4),]
#change column names
names(ploidy)<-c("ploidy","relativeCvg")
#number according to levels
ploidy$ploidy<-0:5 #change this to reflect how many levels are in the header
ploidy$relativeCvg<-as.numeric(as.character(ploidy$relativeCvg))

#Now read data again
cnvDetails<-read.table("cnvTumorDetailsBeta-XXX-ASM.tsv", header=T, sep="\t", na.string="N")

# export PNG
create.mht.plot(cnvDetails, cell.type="T", plot.name="tumorCNVmanhattan.png", max.ploidy=10, chr=NA, title="Tumor Example", beg=NA, end=NA, ploidy=ploidy, plot.type="png")

# export PDF
create.mht.plot(cnvDetails, cell.type="T", plot.name="tumorCNVmanhattan.pdf", max.ploidy=10, chr=NA, title="Tumor Example", beg=NA, end=NA, ploidy=ploidy, plot.type="pdf")


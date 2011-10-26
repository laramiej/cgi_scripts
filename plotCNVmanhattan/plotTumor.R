#!/usr/bin/Rscript

source('plotCNVmanhattan.R')
ploidy<-read.table("/Users/gtyrelle/FAS/projects/customer/emc_dinjens/cnvTumorDetailsBeta-GS000004997-ASM.tsv", header=T, sep="\t", na.string="N", nrows=10, comment.char="")
#change the number of rows read (nrows) to capture all of the MEAN_LEVEL fields
ploidy<-ploidy[-(1:4),]
names(ploidy)<-c("ploidy","relativeCvg")
ploidy$ploidy<-0:5 #change this to reflect how many levels are in the header
ploidy$relativeCvg<-as.numeric(as.character(ploidy$relativeCvg))

cnvDetails<-read.table("/Users/gtyrelle/FAS/projects/customer/emc_dinjens/cnvTumorDetailsBeta-GS000004997-ASM.tsv", header=T, sep="\t", na.string="N")

# export PNG
create.mht.plot(cnvDetails, cell.type="T", plot.name="tumorCNVmanhattan.png", max.ploidy=10, chr=NA, title="Tumor Example", beg=NA, end=NA, ploidy=ploidy, plot.type="png")

# export PDF
create.mht.plot(cnvDetails, cell.type="T", plot.name="tumorCNVmanhattan.pdf", max.ploidy=10, chr=NA, title="Tumor Example", beg=NA, end=NA, ploidy=ploidy, plot.type="pdf")


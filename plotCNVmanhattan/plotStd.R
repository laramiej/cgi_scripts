source('CNVmanhattan.R')

cnvDetails<-read.table("cnvTumorDetailsBeta-XXX-ASM.tsv", header=T, sep="\t", na.string="N")

create.mht.plot(cnvDetails, cell.type="N", plot.name="cnvnormal", max.ploidy=10, chr=NA, title="CNV Manhattan for Normal", beg=NA, end=NA, ploidy=, plot.type="pdf")


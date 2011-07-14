# Jason M Laramie - jlaramie@completegenomics.com
# plotCNVmanhattan.R will create a manhattan style plot using the relative coverage (see attached plot from NA19240)
# Create CNV plots across a whole genome or by a single chromosome
# Need to read in CNV-details file as:

This script has a couple steps that need to be run to utilize it.

Step 1:
The cnvDetails file needs to be read into R using a command similar to for non-tumor samples:
	
cnvDetails<-read.table("cnvDetailsBeta-*-ASM.tsv", header=T, sep="\t", na.string="N")

#if the chromosome is not specified a whole genome plot will be created

For tumors you will need to create a dataframe from the header of the cnvTumorDetailsBeta file
#this header looks like:

##ASSEMBLY_ID    GSXXXXXXXX-ASM
#FORMAT_VERSION 1.5
#GENERATED_AT   X
#GENERATED_BY   CallCNVs
#GENOME_REFERENCE       NCBI build 37
#MEAN_LEVEL_0   0.039
#MEAN_LEVEL_1   0.378
#MEAN_LEVEL_2   0.723
#MEAN_LEVEL_3   0.922
#MEAN_LEVEL_4   1.055
#MEAN_LEVEL_5   1.239
#MEAN_LEVEL_6   1.400
#MEAN_LEVEL_7   1.834
#MEAN_LEVEL_8   3.008
#

Use a command like:

ploidy<-read.table("cnvTumorDetailsBeta-GSXXXXXXXXXX-ASM.tsv", header=T, sep="\t", na.string="N", nrows=13, comment.char="")
# change the number of rows read (nrows) to read in all of the MEAN_LEVEL fields (see above)
#remove the top 4 rows
ploidy<-ploidy[-(1:4),]
#add column names
names(ploidy)<-c("ploidy","relativeCvg")
#add levels to the dataframe
ploidy$ploidy<-0:8 #change this to reflect how many levels are in the header
#change relativeCvg from a factor to numeric
ploidy$relativeCvg<-as.numeric(as.character(ploidy$relativeCvg))


#Step 2:
Now you can run the R script using:
create.mht.plot(dat, cell.type="N", plot.name, max.ploidy=10, chr=NA, title="", beg=NA, end=NA, ploidy=NA, plot.type="pdf")

Where the arguments are as follows:
#dat = is the cnvDetailsBeta-*-ASM.tsv for normal samples or cnvTumorDetailsBeta-GSXXXXXXXX-ASM.tsv read into R using read.table and using the parameter na.strings="N" read in using step 1

#cell.type = "N" for normal or "T" for tumor
#plot.name = what name do you want to give to the outputted plot in quotes
#max.ploidy = this number is given in the header of the cnvDetailsBeta file.  Not needed for tumors
#chr = if you want to make a plot of just a single chr specify it here as "chr2"
#title = what title do you want to appear at the top of the plot, specified in quotes
#beg = chromosome position at the is the start of a horizontal red line 
#end = chromosome position at the is the end of a horizontal red line.
#		beg and end are used if one is trying to visualize a specific area (e.g. SV) overlaid on the CNV data
#ploidy = this is only needed for tumors and is a dataframe created using the header of the cnvTumorDetailsBeta as specified above
#plot.type = default is to create a pdf.  Will also take "png" but I think this doesn't work due to margins.  
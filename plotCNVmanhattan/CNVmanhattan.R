library(snpMatrix)
library(RColorBrewer)
library(hgu95av2.db)

# Jason M Laramie - jlaramie@completegenomics.com
# plotCNVmanhattan.R will create a manhattan style plot using the relative coverage (see attached plot from NA19240)
# Create CNV plots across a whole genome or by a single chromosome
# Need to read in CNV-details file as:
# cnvDetails<-read.table("cnvDetailsBeta-*-ASM.tsv", header=T, sep="\t", na.string="N")
#
#if the chromosome is not specified a whole genome plot will be created

#for tumors you will need to create a dataframe from the header of the cnvTumorDetailsBeta file
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
#
# ploidy<-read.table("/Users/jlaramie/Desktop/cnvTumorDetailsBeta-GSXXXXXXXXXX-ASM.tsv", header=T, sep="\t", na.string="N", nrows=13, comment.char="")
# change the number of rows read (nrows) to capture all of the MEAN_LEVEL fields
# ploidy<-ploidy[-(1:4),]
# names(ploidy)<-c("ploidy","relativeCvg")
# ploidy$ploidy<-0:8 #change this to reflect how many levels are in the header
# ploidy$relativeCvg<-as.numeric(as.character(ploidy$relativeCvg))
#

create.mht.plot<-function(dat, cell.type="N", plot.name, max.ploidy=10, chr=NA, title="", beg=NA, end=NA, ploidy=NA, plot.type="pdf"){

#dat = is the cnvDetailsBeta-*-ASM.tsv for normal samples or cnvTumorDetailsBeta-GSXXXXXXXX-ASM.tsv read into R using read.table and using the parameter na.strings="N"
#		this can be created using the above code for creating the dataframe
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

	names(dat)[1]<-"chromosome"
	dat$cex<-0.5

	if(cell.type=="N"){
		#column name for first datapoints
		y.axis1="avgNormalizedCvg"
		#column name for second datapoints
		y.axis2="calledPloidy"
	
		#for the second axis we need to know the average value for the first datapoints
		cnv.means<-unlist(by(dat[,y.axis1], dat[,y.axis2], mean))
		print(cnv.means)
		#set plotting point sizes
	}else{
		#column name for first datapoints
		y.axis1="relativeCvg"
		#column name for second datapoints
		y.axis2="calledLevel"
	}
	
	
	chrs.num<-c(1:22,"X","Y","M")
	dat$chr.num<-sub("chr","",dat$chromosome)
	dat$chrfactnum<-factor(dat$chr.num,levels=chrs.num)
	chrs<-c(paste("chr",seq(1,22),sep=""),"chrX","chrY","chrM")
	dat$chrfact<-factor(dat$chromosome,levels=chrs)

	dat$cols<-"grey25"

	#alternate between to colors for better clarity between chromosomes
	dat$cols[which((as.numeric(dat$chrfactnum) %% 2) == FALSE)] <-"gray50"

	# Borrow chromosome length function from the affy array annotations
	# Calculate offsets
	chrlens<-hgu95av2CHRLENGTHS[c(1:22,"X","Y","M")]
	offsets<-c(0,cumsum(as.numeric(chrlens)))

	names(offsets)=c(paste("chr", seq(1,22), sep=""),"chrX","chrY","chrM")
	dat$x<-dat$position + offsets[dat$chrfact]
	
	if(cell.type=="N"){
		dat$adj.ploidy<-cnv.means[dat[,y.axis2]+1]
		if(plot.type=="pdf"){
			pdf(file=plot.name, width=20, height=7)
		}else{
			png(file=plot.name, width=20, height=7)
		}
	
		#if a chromosome is given use this code
		if(!is.na(chr)){
			par(yaxt="s", mai=c(1,1,1,1))
			dat<-dat[dat$chromosome==chr,]
			y.max = max(dat[,y.axis1],na.rm=T) + .1
			x.max<-max(dat$position)+100
			x.min<-min(dat$position)-10
			plot(dat$position, dat[,y.axis1], col="grey50", pch=19, ylim=c(0,y.max),xlim=c(0,x.max),cex=dat$cex, xlab=chr, ylab=y.axis1, main=title)
			abline(h=cnv.means, lty=1, col="lightgrey")
			segments(x0=beg,y0=200,x1=end,y1=200, col="red", lwd=2)
			par(new=T, xaxt="n", yaxt="n")
			plot(dat$position, dat$adj.ploidy, col="blue", pch=19, ylim=c(0,y.max),xlim=c(0,x.max),cex=dat$cex, xlab="", ylab="")
			mtext("Ploidy", 4, adj=1, col="blue", at=mean(c(0,y.max)), line=3)
			#no chromosome given.  Plot the whole genome
		}else{
			y.max = max(dat[,y.axis1],na.rm=T) + .1
			x.min<-min(dat$x)-10
			x.max<-max(dat$x)+10
			print(x.min)
			par(xaxt="n",yaxt="s", mar=c(3,4,4,6)+0.1, xaxs="i")
			plot(dat$x, dat[,y.axis1], col=dat$cols, pch=19, ylim=c(0,y.max),xlim=c(0,x.max),cex=dat$cex, xlab="", ylab=y.axis1, main=title)
			abline(h=cnv.means, lty=1, col="lightgrey")
			par(new=T, xaxt="n", yaxt="n", xaxs="i")
			plot(dat$x, dat$adj.ploidy, col="blue", pch=19, ylim=c(0,y.max),xlim=c(0,x.max),cex=dat$cex, xlab="", ylab="")

			mid<-apply(cbind(offsets[1:25]+1,offsets[2:26]),1,mean)
			par(xaxt = "s")
			#create labels at a 45 degree angle
			text(x=mid,y=-27, labels=chrs, xpd=TRUE, srt=45, adj = 1, cex=1)
			#seperate the chromosomes by a grey vertical line
			for(i in 2:25){
      			abline(v=offsets[i],lty=1, col="lightgrey")
			}
			#add a label to the second axis
			mtext("Ploidy", 4, adj=1, col="blue", at=mean(c(0,y.max)), line=3.5)
		}
		#label the second y-axis with ploidy number and the average normalized coverage
		mtext(paste(names(cnv.means)," (",round(cnv.means),")", sep=""), 4, at=cnv.means, col = "blue", las=1, line=.25)
		dev.off()
	}else{
		#plot for tumor CNV
		#dat$adj.ploidy<-cnv.means[dat[,y.axis2]+1]
		pdf(file=plot.name, width=20, height=7)
		#if a chromosome is given use this code
		if(!is.na(chr)){
			par(yaxt="s", mai=c(1,1,1,1))
			dat<-dat[dat$chromosome==chr,]
			y.max = max(dat[,y.axis1],na.rm=T) + .1
			x.max<-max(dat$position)+100
			x.min<-min(dat$position)-10
			plot(dat$position, dat[,y.axis1], col="grey50", pch=19, ylim=c(0,y.max),xlim=c(0,x.max),cex=dat$cex, xlab=chr, ylab=y.axis1, main=title)
			abline(h=ploidy$relativeCvg, lty=1, col="lightgrey")
			if(!is.na(beg)){segments(x0=beg,y0=200,x1=end,y1=200, col="red", lwd=2)}
			par(new=T, xaxt="n", yaxt="n")
			plot(dat$position, dat[,y.axis2], col="blue", pch=19, ylim=c(0,y.max),xlim=c(0,x.max),cex=dat$cex, xlab="", ylab="")
			mtext("Level", 4, adj=1, col="blue", at=mean(c(0,y.max)), line=1.5)
			#no chromosome given.  Plot the whole genome
		}else{
			y.max = max(ploidy[,2])+.25
			#y.max = max(dat[,y.axis1],na.rm=T)
			print (y.max)
			x.min<-min(dat$x)-10
			x.max<-max(dat$x)+10
			par(xaxt="n",yaxt="s", mar=c(4,4,4,5)+0.1, xaxs="i")
			plot(dat$x, dat[,y.axis1], col=dat$cols, pch=19, ylim=c(0,y.max),xlim=c(0,x.max),cex=dat$cex, xlab="", ylab=y.axis1, main=title)
			print(as.character(ploidy$relativeCvg))
			abline(h=ploidy$relativeCvg, lty=1, col="lightgrey")
			par(new=T, xaxt="n", yaxt="n", xaxs="i")
			plot(dat$x, dat[,y.axis2], col="blue", pch=19, ylim=c(0,y.max),xlim=c(0,x.max),cex=dat$cex, xlab="", ylab="")
			mid<-apply(cbind(offsets[1:25]+1,offsets[2:26]),1,mean)
			par(xaxt = "s")
			#create labels at a 45 degree angle
			text(x=mid,y=-.3, labels=chrs, xpd=TRUE, srt=45, adj = 1, cex=1)
			#seperate the chromosomes by a grey vertical line
			for(i in 2:25){
      			abline(v=offsets[i],lty=1, col="lightgrey")
			}
			#add a label to the second axis
			mtext("Level", 4, adj=1, col="blue", at=mean(c(0,y.max)), line=1.5)
		}
		#label the second y-axis with ploidy number and the average normalized coverage
		mtext(ploidy$ploidy, 4, at=ploidy$relativeCvg, col = "blue", las=1, line=.25)
		dev.off()
	}
}
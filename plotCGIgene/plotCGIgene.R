library(gplots)

plotNoCoverage <- function(x, ...){
	if(as.numeric(x[6])<7){
		abline(v=x[1],...)
	}
}

plotCNV <- function(x, ...){
		segments(x[1],x[3],x[2],x[3],...)}

plotFeaturesClass <- function(x, ...){
	rect(x[19],x[20],x[21],x[22],col=x[15],...)}
plotFeaturesZygosity <- function(x, ...){
	rect(x[19],x[20],x[21],x[22],col=x[16],...)}
plotFeaturesComp <- function(x, ...){
	rect(x[19],x[20],x[21],x[22],col=x[17],...)}
plotFeaturesFunc <- function(x, ...){
	rect(x[19],x[20],x[21],x[22],col=x[18],...)}
lociText <- function(x,...){
	mtext(x[1],side=1,at=(as.numeric(x[19])+.5), ...)}
plotExons <- function(x, ...){
	rect(x[4],x[6],x[5],x[7],col="grey",lwd=.5,...)}
plotIntrons <- function(x, ...){
	segments(x[1],x[3],x[2],x[4],col="black", lwd=.25,...)}
plotIntrons <- function(x, ...){
	segments(x[1],x[3],x[2],x[4],col="black", lwd=.25,...)}

plotCGIgene<-function(file.name, gene.name, no.call.freq=0, half.call.num=0, coverage.file, allele.balance.file, cnv.file, variant.file, transcript.file, start=NA, end=NA){
	
	#nf <- layout(matrix(c(1,2,3,4,5,0),ncol=1,byrow=F),heights=c(1,1,1,1,1,.5), TRUE)
	pdf(file=file.name)
	nf<- layout(matrix(c(1,2,3,4,5,6,7,7,8,8,9,9),ncol=2,byrow=F),heights=c(.5,.5,.5,.5,1,1),widths=c(3.5,.5), TRUE)
	#layout.show(nf)
	par(mar=c(0,4,2,1))

	#coverage plot
	dat<-read.table(coverage.file, sep="\t", header=F)	
	myX<-dat[,1]
	myY<-dat[,6]
	min.x=min(dat[,1])
	max.x=max(dat[,1])
	eval<-(max(myX)-min(myX))+1
	myFit<-loess.smooth(myX,myY,evaluation=eval,span=0.06,degree=2)
	plot(myX,myY,xlim=c(min.x,max.x), ylab="Coverage", cex=.1, pch=20,xaxt="n", main=gene.name,cex.axis=.6,cex.lab=.7)
	lines(myFit,col="red")
	abline(h=7, col="blue")
	apply(dat,1,plotNoCoverage,lwd=.1,col="lightgrey")
	mtext(paste("Mean Coverage = ", round(mean(dat[,6]),2),sep=""),3,adj=1, cex=.50)
	mtext(paste("No-call Freq (bases) = ", no.call.freq, "%  Half-call allele # = ", half.call.num ,sep=""),3,adj=0, cex=.50)
	
	#allele balance
	if(file.info(allele.balance.file)$size > 0){
		dat<-read.table(allele.balance.file, sep="\t", header=F)
		par(mar=c(0,4,.5,1))
		plot(dat[,2],dat[,6],xlim=c(min.x,max.x), type="h", cex=.1,ylim=c(0,.5), xlab="", ylab="Het Ratio",xaxt="n",cex.axis=.6,cex.lab=.7)
	}

	#plot CNV data
	dat<-read.table(cnv.file, sep="\t", header=F, na.strings="N")
	dat$mean<-(dat[,1]+dat[,2])/2
	y.max<-max(dat[,3])+1
	dat$col="blue"
	dat[dat[,4]=="H","col"]<-"black"
	dat[dat[,4]=="I","col"]<-"black"
	par(mar=c(0,4,.5,1))
	plot(dat[,1],dat[,3],xlim=c(min.x,max.x), ylim=c(0,y.max), ylab="Rel. Coverage", type="n", xaxt="n",cex.axis=.6,cex.lab=.7)
	text(dat$mean,dat[,3],labels=dat[,4],col=dat$col)

	#plot variants scores
	dat<-read.table(variant.file, sep="\t", header=F)
	par(mar=c(0,4,0,0))
	barplot2(as.matrix(t(dat[,10:11])), beside = F, col = c("red", "blue"), ylab="TotalScore", cex.axis=.6,cex.lab=.7)
	legend("topleft",col=c("red","blue"),pch=19,cex=.5,legend=c("Allele1","Allele2"))

	#configure the colors
	dat$col.class<-"black"
	dat[dat[,6]=="snp","col.class"]<-"green"
	dat[dat[,6]=="ins","col.class"]<-"red"
	dat[dat[,6]=="del","col.class"]<-"blue"
	dat[dat[,6]=="sub","col.class"]<-"orange"
	dat[dat[,6]=="complex","col.class"]<-"black"

	dat$col.zygosity<-"black"
	dat[dat[,5]=="hap","col.zygosity"]<-"green"
	dat[dat[,5]=="hom","col.zygosity"]<-"red"
	dat[dat[,5]=="het-ref","col.zygosity"]<-"blue"
	dat[dat[,5]=="het-alt","col.zygosity"]<-"orange"

	dat$col.comp<-"white" #missing
	dat[dat[,12]=="CDS","col.comp"]<-"black"
	dat[dat[,12]=="INTRON","col.comp"]<-"green"
	dat[dat[,12]=="DONOR"|dat[,12]=="ACCEPTOR","col.comp"]<-"red"
	dat[dat[,12]=="TSS-UPSTREAM","col.comp"]<-"blue"
	dat[dat[,12]=="SPAN5"|dat[,12]=="SPAN3"|dat[,12]=="SPAN","col.comp"]<-"orange"
	dat[dat[,12]=="UTR5"|dat[,12]=="UTR3"|dat[,12]=="UTR","col.comp"]<-"grey"

	dat$col.func<-"tan" #missing
	dat[dat[,13]=="COMPATIBLE","col.func"]<-"black"
	dat[dat[,13]=="NO-CHANGE","col.func"]<-"black"
	dat[dat[,13]=="MISSENSE","col.func"]<-"green"
	dat[dat[,13]=="NONSENSE","col.func"]<-"red"
	dat[dat[,13]=="NONSTOP","col.func"]<-"blue"
	dat[dat[,13]=="DELETE"|dat[,13]=="DELETE+","col.func"]<-"orange"
	dat[dat[,13]=="INSERT"|dat[,13]=="INSERT+","col.func"]<-"purple"
	dat[dat[,13]=="FRAMESHIFT","col.func"]<-"brown"
	dat[dat[,13]=="MISSTART","col.func"]<-"yellow"
	dat[dat[,13]=="DISRUPT","col.func"]<-"pink"
	dat[dat[,13]=="UNKNOWN-VNC","col.func"]<-"lightgrey"
	dat[dat[,13]=="UNKNOWN-INC","col.func"]<-"darkgrey"
	dat[dat[,13]=="UNKNOWN-TR","col.func"]<-"greenyellow"
	#print(head(dat))
	#NO-CHANGE,COMPATIBLE,MISSENSE,NONSENSE,NONSTOP,DELETE,INSERT,DELETE+,INSERT+,FRAMESHIFT,MISSTART,DISRUPT,UNKNOWN-VNC,UNKNOWN-INC,UNKNOWN-TR

	par(mar=c(0,4,0,0))
	plot(0:nrow(dat),c(5,rep(0,nrow(dat))), type="n", yaxt="n", xaxt="n",bty="n", xlab="", ylab="", main="")

	dat$x1<-0:(nrow(dat)-1)
	dat$y1 = 0
	dat$x2 = dat$x1+1
	dat$y2 = 1
	apply(dat,1,plotFeaturesFunc)

	dat$x1<-0:(nrow(dat)-1)
	dat$y1 = 1
	dat$x2 = dat$x1+1
	dat$y2 = 2
	apply(dat,1,plotFeaturesComp)

	dat$x1<-0:(nrow(dat)-1)
	dat$y1 = 2
	dat$x2 = dat$x1+1
	dat$y2 = 3
	apply(dat,1,plotFeaturesZygosity)

	#class
	dat$x1<-0:(nrow(dat)-1)
	dat$y1 = 3
	dat$x2 = dat$x1+1
	dat$y2 = 4
	apply(dat,1,plotFeaturesClass)

	mtext("Class", side=2, at=3.5, las=1, cex=.5)
	mtext("Zygosity", side=2, at=2.5, las=1, cex=.5)
	mtext("Compartment", side=2, at=1.5, las=1, cex=.5)
	mtext("Function", side=2, at=0.5, las=1, cex=.5)
	mtext("Var Length", side=2, at=4.5, las=1, cex=.5)

	text((dat[,19]+.5),4.5,labels=dat[,14], cex=.4)

	dat$comp.het<-0
	dat[dat[,13]=="MISSENSE","comp.het"]<-1
	dat[dat[,13]=="NONSENSE","comp.het"]<-1
	dat[dat[,13]=="NONSTOP","comp.het"]<-1
	dat[dat[,13]=="DELETE"||dat[,13]=="DELETE+","comp.het"]<-1
	dat[dat[,13]=="INSERT"||dat[,13]=="INSERT+","comp.het"]<-1
	dat[dat[,13]=="FRAMESHIFT","comp.het"]<-1
	dat[dat[,13]=="MISSTART","comp.het"]<-1
	dat[dat[,13]=="DISRUPT","comp.het"]<-1
	dat[dat[,12]=="DONOR" || dat[,12]=="ACCEPTOR","comp.het"]<-1 #capture all potential splice
	#plot the transcripts
	dat.tx<-read.table(transcript.file, sep="\t", header=F)
	tx<-unique(dat.tx[,1])
	
	par(mar=c(5,4,0,0))
	plot(min.x:max.x,c(length(tx)+1,rep(0,length(min.x:(max.x-1)))),xlim=c(min.x,max.x), type="n", yaxt="n", bty="n", ylab="", xlab="Base Pair")
	for(i in 1:length(tx)){
		new.tx<-dat.tx[dat.tx[,1] == tx[i],]
		new.tx$y1<-i-.5
		new.tx$y2<-i
		apply(new.tx,1,plotExons)
		intron<-as.data.frame(cbind(new.tx[-nrow(new.tx),5],new.tx[2:nrow(new.tx),4]))
		intron$y1<-i-.25
		intron$y2<-i-.25
		apply(intron,1,plotIntrons)
		mtext(tx[i], side=2, at=i-.25, las=1, cex=.5)
	}
	
	#recalulate scale 
	diffscale <- (max.x-min.x)/(nrow(dat))
	new.lines<-NULL
	for(i in 1:nrow(dat)){
		for(j in 1:length(tx)){
			lines(rep(dat[i,3],2),c(j-0.3,j-0.2), col=dat[i,"col.func"])
		}
		segments(x0=as.numeric(dat[i,3]),y0=(length(tx)-.1),x1=as.numeric((min.x+(i)*diffscale-.5*diffscale)),y1=length(tx)+1, col=dat[i,"col.func"])
	}
	mtext(paste("Dmg. Vars = ", sum(dat$comp.het),sep=""), side=1, cex=.60, col="red")
	if(sum(dat$comp.het) >= 1){
		print(paste("ALERT: Found ",sum(dat$comp.het)," Damaging Variants", sep=""))
	} 
	#add figure legends
	par(mar=c(0,0,0,0))
	plot.new()
	legend("bottomleft",col=c("green","red","blue","orange","black"),pch=19,cex=.8,legend=c("snp","ins","del","sub","complex"), title="Class")

	plot.new()
	legend("left",col=c("green","red","blue","orange"),pch=19,,cex=.8,legend=c("Haploid","Homo","het-ref","het-alt"),title="Zygosity")

	plot.new()
	legend("topleft",col=c("black","green","red","blue","orange","grey","white"),pch=19,cex=.75,legend=c("CDS","Intron","SPLICE","UPSTREAM","SPAN","UTR","NA"), 	title="Component")

	legend("left",col=c("black","green","red","blue","orange","purple","brown","yellow","pink","lightgrey","darkgrey","greenyellow","tan"),pch=19,legend=c("Syn","Missense","Nonsense","Nonstop","Delete","Insert", "Frameshift","Misstart","Disrupt","Unknown-VNC","Unknown-INC","Unknown-TR","NA"), cex=.7,title="Function")

dev.off()
}
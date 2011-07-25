# graph_summary_data.r
# 18/07/2011

# Andrew R Wood
# Genetics of Complex Traits
# Peninsula Medical School
# University of Exeter
# Exeter EX1 2LU
# United Kingdom
#

#Jason Laramie
# jlaramie@completegenomics.com

# Instructions: load source and enter top directory or a sample.list (see below) containing ALL CGI directories to be scanned for summary files as the argument of the function
##Where sample list is a space delinated file that must have the column heading "dir" like and "DNA":
#proj     DNA    cust                 ver      GB  PB   dir
#XXXX   Sample1  XXXX                 1.11     37  a    /path/to/genome/root4sample1
#XXXX   Sample2  XXXX                 1.11     37  a    /path/to/genome/root4sample2


graph.summary.data <- function(top_dir, sample.list=NULL, assembly.version=1.12, test.pairwise=F) {
    require(lattice)
    require(grid)
    file_list=NULL
    # Pull out listing of summary files within sub-directories. 
    # Final list should be sorted by summary-* name so files have the same index across different file systems and directory structures
    print("Searching for summary-* files...")
    if(!is.null(sample.list)){
    	urls<-read.table(sample.list, header=T)
    	for(i in 1:nrow(urls)){
    		#print(paste("Working on ", as.character(urls[i,"DNA"], sep="")))
			new.url<-paste(urls[i,"dir"],"/ASM/", sep="")
			summary.file<-list.files(new.url,pattern="summary-")
			file_list[i] = paste(new.url,summary.file,sep="")
    	}
    	#print (file_list)
    	file_list <- funct_sort_files(file_list)	
    }else{
    	file_list <- funct_get_files(top_dir)
    	#print(file_list)
    }
    file_count <- length(file_list)
   	cat(file_count, "files found\n", sep=" ")

    if (file_count > 0) {

        # Generate the data frame
        print("Extracting data from summary files")
        summary_dataframe <- funct_fill_data_frame(file_list, assembly.version)
        attach(summary_dataframe)
    
        # Plot the data
        # Define X axis range:
        xlim <- c(0,max(File_Index))
        # Define lattice plot grid line markers:
        d <- seq(0, max(File_Index), by=5)
        # Define x labels dependent on sample size
        if(max(File_Index) <= 500) {
             xaxis=list(at = seq(0, max(File_Index), by=25))
             xaxis_ticks   = seq(0, max(File_Index), by=25)
        }
        else {
            xaxis=list(at = seq(0, max(File_Index), by=50))
            xaxis_ticks   = seq(0, max(File_Index), by=50)
        }
        print("Plotting data")

        # Plot data regarding fractions called for genome/exome
        funct_plot_genome_fractions(xlim,xaxis,d)
        funct_plot_genome_fractions_merged(xaxis_ticks)
        funct_plot_genome_fractions_hists()
        funct_plot_genome_fractions_box()
        # Plot data regarding coverage across the genome
        funct_plot_genome_coverage(xlim,xaxis,d)
        funct_plot_genome_coverage_merged(xaxis_ticks)
        funct_plot_genome_coverage_hists()
        funct_plot_genome_coverage_box()
        # Plot data regarding yields and library mate means
        funct_plot_genome_yields_cov_var(xlim,xaxis,d)
        funct_plot_genome_yields_cov_var_hists()
        funct_plot_genome_yields_cov_var_box()
        # Plot library mate mean data
        funct_plot_library_mate_data(xaxis_ticks)
        # Plot SNP stats
        funct_plot_genome_snp_stats(xlim,xaxis,d)
        funct_plot_genome_snp_stats_hists()
        funct_plot_genome_snp_stats_box()
        funct_plot_exome_snp_stats(xlim,xaxis,d)
        funct_plot_exome_snp_stats_hists()
        funct_plot_exome_snp_stats_box()
        # Plot INS stats
        funct_plot_ins_stats(xlim,xaxis,d)
        funct_plot_ins_stats_hists()
        funct_plot_ins_stats_box()
        # Plot DEL stats
        funct_plot_del_stats(xlim,xaxis,d)
        funct_plot_del_stats_hists()
        funct_plot_del_stats_box()
        # Plot SUB stats
        funct_plot_sub_stats(xlim,xaxis,d)
        funct_plot_sub_stats_hists()
        funct_plot_sub_stats_box()
        # Plot predicted function stats
        funct_plot_func_impact_snps(xlim,xaxis,d)
        funct_plot_func_impact_snps_hists()
        funct_plot_func_impact_snps_box()
        funct_plot_func_impact_framing(xlim,xaxis,d)
        funct_plot_func_impact_framing_hists()
        funct_plot_func_impact_framing_box()
        # Plot CNV stats
        funct_plot_cnvs(xlim,xaxis,d)
        funct_plot_cnvs_hists()
        funct_plot_cnvs_box()
        # Plot SV Ssts
        funct_plot_sv(xlim,xaxis,d)
        funct_plot_sv_hists()
        funct_plot_sv_box()
        #plot MEI data
        if(assembly.version >= 1.12){
        funct_plot_genome_fractions(xlim,xaxis,d)
        	funct_plot_mei(xlim,xaxis,d)
        }
        
        #find correlations among the data
        if(test.pairwise){
        	funct_find_corr(summary_dataframe)
        }
        dev.off()
        detach(summary_dataframe)
    }
    
    else {
        print("No summary files found")
    }
    
    rm(list=ls(all=TRUE))
    
}

    
funct_get_files <- function(top_dir) {
    # Get full path of summary files
    file_list <- list.files(path=top_dir, full.name=TRUE, recursive=TRUE, pattern="summary-")
    if (length(file_list)>0) {
        # Call function to sort the summary files by summary-* name.
        file_list <- funct_sort_files(file_list)
        # Return the file list ordered
    }
    return(file_list)
}    



funct_sort_files <- function(file_list) {

    # Generate vector to hold summary files in order
    summary_file_names = NULL
    for(i in 1:length(file_list)) {
        startindex = NULL
        this_file = file_list[i]       
        for(j in nchar(this_file):1) {
            if ((substr(this_file,start=j, stop=j)) == "/") {
                startindex = j+1
                break
            }
        }
        summary_file_names = c(summary_file_names,substr(this_file,start=startindex,stop=nchar(this_file))) 
    }
    # sort!
    summary_file_names = sort(summary_file_names)

    # Now grep the file_list vector in order of the summary_file_names to assign new indicies
    file_list_ordered = NULL
    for (i in 1:length(summary_file_names)) {
        this_file = summary_file_names[i]
        for (j in 1:length(file_list)) {
            if (grepl(this_file,file_list[j],fixed=TRUE)) {
                file_list_ordered[i] = file_list[j]
            }
        }
    }
    # return the file list
    return(file_list_ordered)
}


funct_fill_data_frame <- function(file_list, assembly.version) {

    summary_dataframe = NULL

    for (i in 1:length(file_list)) {
        summary_file_as_df <- read.table(file_list[i], header=TRUE, sep="\t")
        if (i == 1) {
            summary_dataframe = summary_file_as_df
        }
        else {
            summary_dataframe = cbind(summary_dataframe, summary_file_as_df$Value)
        }
        names(summary_dataframe)[i+2] = file_list[i]    
    }
    
    file_indicies=seq(1:length(file_list))
    trans_summary_dataframe = as.data.frame(t(summary_dataframe),stringsAsFactors=FALSE)
    trans_summary_dataframe = trans_summary_dataframe[-c(1,2),]
    trans_summary_dataframe = cbind(file_indicies,file_list,trans_summary_dataframe)
    row.names(trans_summary_dataframe) <- NULL 
	
	if(assembly.version==1.11){
    	names(trans_summary_dataframe) = c(
                                     "File_Index",
                                     "File",
                                     "Gender",
                                     "Genome_Fully_Called_Fraction",  
                                     "Genome_Partially_Called_Fraction",  
                                     "Genome_No_Called_Fraction",
                                     "Genome_Gross_Mapping_Yield_Gb", 
                                     "Genome_Both_Mates_Mapped_Yield_Gb", 
                                     "Genome_100k_Normalised_Coverage_Variability",
                                     "Genome_Fraction_WeightSumSequenceCoverage_gt_5x", 
                                     "Genome_Fraction_WeightSumSequenceCoverage_gt_10x",
                                     "Genome_Fraction_WeightSumSequenceCoverage_gt_20x",
                                     "Genome_Fraction_WeightSumSequenceCoverage_gt_30x",
                                     "Genome_Fraction_WeightSumSequenceCoverage_gt_40x",
                                     "Exome_Fully_Called_Fraction", 
                                     "Exome_Partially_Called_Fraction", 
                                     "Exome_No_Called_Fraction",
                                     "Exome_Fraction_WeightSumSequenceCoverage_gt_5x",
                                     "Exome_Fraction_WeightSumSequenceCoverage_gt_10x",
                                     "Exome_Fraction_WeightSumSequenceCoverage_gt_20x",
                                     "Exome_Fraction_WeightSumSequenceCoverage_gt_30x",
                                     "Exome_Fraction_WeightSumSequenceCoverage_gt_40x",
                                     "Library_Mate_Distribution_Mean", 
                                     "Library_Mate_Distribution_Range_95CI_Min", 
                                     "Library_Mate_Distribution_Range_95CI_Max",
                                     "Genome_var_SNP_Total_Count",
                                     "Genome_var_SNP_Homozygous_Count",
                                     "Genome_var_SNP_Heterozygous_Count",
                                     "Genome_var_SNP_Novel_Fraction",
                                     "Genome_var_SNP_Homozygous_Novel_Fraction",
                                     "Genome_var_SNP_Heterozygous_Novel_Fraction",
                                     "Genome_var_SNP_Het_Hom_Ratio",
                                     "Genome_var_SNP_Transition_Transversion_Ratio",
                                     "Genome_var_INS_Total_Count",
                                     "Genome_var_INS_Novel_Fraction",
                                     "Genome_var_INS_Het_Hom_Ratio",
                                     "Genome_var_DEL_Total_Count",
                                     "Genome_var_DEL_Novel_Fraction",
                                     "Genome_var_DEL_Het_Hom_Ratio",
                                     "Genome_var_SUB_Total_Count",
                                     "Genome_var_SUB_Novel_Fraction",
                                     "Genome_var_SUB_Het_Hom_Ratio",
                                     "Exome_var_SNP_Total_Count",
                                     "Exome_var_SNP_Homozygous_Count",
                                     "Exome_var_SNP_Heterozygous_Count",
                                     "Exome_var_SNP_Novel_Fraction",
                                     "Exome_var_SNP_Homozygous_Novel_Fraction",
                                     "Exome_var_SNP_Heterozygous_Novel_Fraction",
                                     "Exome_var_SNP_Het_Hom_Ratio",
                                     "Exome_var_SNP_Transition_Transversion_Ratio",
                                     "Exome_var_INS_Total_Count",
                                     "Exome_var_INS_Novel_Fraction",
                                     "Exome_var_INS_Het_Hom_Ratio",
                                     "Exome_var_DEL_Total_Count",
                                     "Exome_var_DEL_Novel_Fraction",
                                     "Exome_var_DEL_Het_Hom_Ratio",
                                     "Exome_var_SUB_Total_Count",
                                     "Exome_var_SUB_Novel_Fraction",
                                     "Exome_var_SUB_Het_Hom_Ratio",
                                     "Functional_impact_Synonymous_SNP_Loci",
                                     "Functional_impact_Non_synonymous_SNP_Loci",
                                     "Functional_impact_Missense_SNP_Loci",
                                     "Functional_impact_Nonsense_SNP_Loci",
                                     "Functional_impact_Nonstop_SNP_Loci",
                                     "Functional_impact_Misstart_SNP_Loci",
                                     "Functional_impact_Disrupt_SNP_Loci",
                                     "Functional_impact_Frame_Shifting_INS_Loci", 
                                     "Functional_impact_Frame_Shifting_DEL_Loci",
                                     "Functional_impact_Frame_shifting_SUB_Loci", 
                                     "Functional_impact_Frame_Preserving_INS_Loci",
                                     "Functional_impact_Frame_Preserving_DEL_Loci", 
                                     "Functional_impact_Frame_Preserving_SUB_Loci",
                                     "CNV_Total_CNV_Segment_Count", 
                                     "CNV_Total_Bases_in_CNV_Segments", 
                                     "CNV_Fraction_Novel_by_Segment_Count",
                                     "CNV_Fraction_Novel_by_Base_Count", 
                                     "SV_Total_Junction_Count", 
                                     "SV_High_Confidence_Junction_Count"
                                    )
	}else if(assembly.version == 1.12){
		names(trans_summary_dataframe) = c(
                                     "File_Index",
                                     "File",
                                     "Gender",
                                     "Genome_Fully_Called_Fraction",  
                                     "Genome_Partially_Called_Fraction",  
                                     "Genome_No_Called_Fraction",
                                     "Genome_Gross_Mapping_Yield_Gb", 
                                     "Genome_Both_Mates_Mapped_Yield_Gb", 
                                     "Genome_100k_Normalised_Coverage_Variability",
                                     "Genome_Fraction_WeightSumSequenceCoverage_gt_5x", 
                                     "Genome_Fraction_WeightSumSequenceCoverage_gt_10x",
                                     "Genome_Fraction_WeightSumSequenceCoverage_gt_20x",
                                     "Genome_Fraction_WeightSumSequenceCoverage_gt_30x",
                                     "Genome_Fraction_WeightSumSequenceCoverage_gt_40x",
                                     "Exome_Fully_Called_Fraction", 
                                     "Exome_Partially_Called_Fraction", 
                                     "Exome_No_Called_Fraction",
                                     "Exome_Fraction_WeightSumSequenceCoverage_gt_5x",
                                     "Exome_Fraction_WeightSumSequenceCoverage_gt_10x",
                                     "Exome_Fraction_WeightSumSequenceCoverage_gt_20x",
                                     "Exome_Fraction_WeightSumSequenceCoverage_gt_30x",
                                     "Exome_Fraction_WeightSumSequenceCoverage_gt_40x",
                                     "Library_Mate_Distribution_Mean", 
                                     "Library_Mate_Distribution_Range_95CI_Min", 
                                     "Library_Mate_Distribution_Range_95CI_Max",
                                     "Genome_var_SNP_Total_Count",
                                     "Genome_var_SNP_Homozygous_Count",
                                     "Genome_var_SNP_Heterozygous_Count",
                                     "Genome_var_SNP_Novel_Fraction",
                                     "Genome_var_SNP_Homozygous_Novel_Fraction",
                                     "Genome_var_SNP_Heterozygous_Novel_Fraction",
                                     "Genome_var_SNP_Het_Hom_Ratio",
                                     "Genome_var_SNP_Transition_Transversion_Ratio",
                                     "Genome_var_INS_Total_Count",
                                     "Genome_var_INS_Novel_Fraction",
                                     "Genome_var_INS_Het_Hom_Ratio",
                                     "Genome_var_DEL_Total_Count",
                                     "Genome_var_DEL_Novel_Fraction",
                                     "Genome_var_DEL_Het_Hom_Ratio",
                                     "Genome_var_SUB_Total_Count",
                                     "Genome_var_SUB_Novel_Fraction",
                                     "Genome_var_SUB_Het_Hom_Ratio",
                                     "Exome_var_SNP_Total_Count",
                                     "Exome_var_SNP_Homozygous_Count",
                                     "Exome_var_SNP_Heterozygous_Count",
                                     "Exome_var_SNP_Novel_Fraction",
                                     "Exome_var_SNP_Homozygous_Novel_Fraction",
                                     "Exome_var_SNP_Heterozygous_Novel_Fraction",
                                     "Exome_var_SNP_Het_Hom_Ratio",
                                     "Exome_var_SNP_Transition_Transversion_Ratio",
                                     "Exome_var_INS_Total_Count",
                                     "Exome_var_INS_Novel_Fraction",
                                     "Exome_var_INS_Het_Hom_Ratio",
                                     "Exome_var_DEL_Total_Count",
                                     "Exome_var_DEL_Novel_Fraction",
                                     "Exome_var_DEL_Het_Hom_Ratio",
                                     "Exome_var_SUB_Total_Count",
                                     "Exome_var_SUB_Novel_Fraction",
                                     "Exome_var_SUB_Het_Hom_Ratio",
                                     "Functional_impact_Synonymous_SNP_Loci",
                                     "Functional_impact_Non_synonymous_SNP_Loci",
                                     "Functional_impact_Missense_SNP_Loci",
                                     "Functional_impact_Nonsense_SNP_Loci",
                                     "Functional_impact_Nonstop_SNP_Loci",
                                     "Functional_impact_Misstart_SNP_Loci",
                                     "Functional_impact_Disrupt_SNP_Loci",
                                     "Functional_impact_Frame_Shifting_INS_Loci", 
                                     "Functional_impact_Frame_Shifting_DEL_Loci",
                                     "Functional_impact_Frame_shifting_SUB_Loci", 
                                     "Functional_impact_Frame_Preserving_INS_Loci",
                                     "Functional_impact_Frame_Preserving_DEL_Loci", 
                                     "Functional_impact_Frame_Preserving_SUB_Loci",
                                     "CNV_Total_CNV_Segment_Count", 
                                     "CNV_Total_Bases_in_CNV_Segments", 
                                     "CNV_Fraction_Novel_by_Segment_Count",
                                     "CNV_Fraction_Novel_by_Base_Count", 
                                     "SV_Total_Junction_Count", 
                                     "SV_High_Confidence_Junction_Count",
                                     "Mobile_element_insertion_count",
                                     "Fraction_of_novel_MEI"
                                    )
	}else{
		stop(paste("Assembly version ", assembly.version, " not known",sep=""))
	}
    trans_summary_dataframe[, c(4:ncol(trans_summary_dataframe))] <- sapply(trans_summary_dataframe[, c(4:ncol(trans_summary_dataframe))], as.numeric)
    write.table(trans_summary_dataframe,file="summary_data_plotted.tsv",sep="\t",row.names=FALSE)
    return(trans_summary_dataframe)
}

funct_plot_genome_fractions <- function(xlim,xaxis,d) {
 
    trellis.device(device="pdf", file="summary_plots.pdf", height=11, width=7.5, color=TRUE)

    #--Define plot titles:
    lab.gfc <- "Genome: Fully Called"
    lab.gpc <- "Genome: Partially Called"
    lab.gnc <- "Genome: Not Called"
    lab.efc <- "Exome: Fully Called"
    lab.epc <- "Exome: Partially Called"
    lab.enc <- "Exome: Not Called"
    
    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.enc, lab.epc, lab.efc, lab.gnc, lab.gpc, lab.gfc)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
        
    plot<-xyplot(Exome_No_Called_Fraction  + Exome_Partially_Called_Fraction  + Exome_Fully_Called_Fraction +
                 Genome_No_Called_Fraction + Genome_Partially_Called_Fraction + Genome_Fully_Called_Fraction ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 6, 1), ylab="Fraction", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("Fraction of Genome/Exome Called Across", length(File_Index), "Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          )
    print(plot)
}


funct_plot_genome_fractions_merged <- function(xaxis_ticks) {

    par(mfrow=c(3,1))
    xgrid<-c(seq(0,max(File_Index)+5,by=5))

    ygrid<-c(seq(0.9,1.0,by=0.02))
    plot(NULL,NULL,type="l",xlim=c(0,max(File_Index)+5),ylim=c(0.9,1), col="blue", main=paste("Fraction of Genome/Exome Fully Called Across",length(File_Index),"Subjects"), xlab="Subject Index", ylab="Fraction", xaxt="n")
    axis(1, at=xaxis_ticks)
    abline(h = ygrid, v = xgrid, col = "grey90", lty=1, lwd=0.5)
    lines(File_Index,Genome_Fully_Called_Fraction, col="blue")
    lines(File_Index,Exome_Fully_Called_Fraction, col="red")
    legend("bottomright", c("Genome","Exome"), cex=0.8, col=c("blue", "red"), lty=1, bty="n")

    ygrid<-c(seq(0,0.02,by=0.005))
    plot(NULL,NULL,type="l",xlim=c(0,max(File_Index)+5),ylim=c(0,0.02), col="blue", main=paste("Fraction of Genome/Exome Partially Called Across",length(File_Index),"Subjects"), xlab="Subject Index", ylab="Fraction", xaxt="n")
    axis(1, at=xaxis_ticks)
    abline(h = ygrid, v = xgrid, col = "grey90", lty=1, lwd=0.5)
    lines(File_Index,Genome_Partially_Called_Fraction,col="blue")
    lines(File_Index,Exome_Partially_Called_Fraction, col="red")
    legend("bottomright", c("Genome","Exome"), cex=0.8, col=c("blue", "red"), lty=1, bty="n")

    ygrid<-c(seq(0,0.06,by=0.01)) 
    plot(NULL,NULL,type="l",xlim=c(0,max(File_Index)+5),ylim=c(0,0.06), col="blue", main=paste("Fraction of Genome/Exome Not Called Across",length(File_Index),"Subjects"), xlab="Subject Index", ylab="Fraction", xaxt="n")
    axis(1, at=xaxis_ticks)
    abline(h = ygrid, v = xgrid, col = "grey90", lty=1, lwd=0.5)
    lines(File_Index,Genome_No_Called_Fraction, col="blue")
    lines(File_Index,Exome_No_Called_Fraction, col="red")
    legend("bottomright", c("Genome","Exome"), cex=0.8, col=c("blue", "red"), lty=1, bty="n")

}

funct_plot_genome_fractions_hists <- function() {

    par(mfrow=c(3,2))   
    hist(Genome_Fully_Called_Fraction,     col="grey", prob=TRUE, main="Histogram: Fraction of Genome Fully Called", xlab="Fraction of Genome")
    hist(Exome_Fully_Called_Fraction,      col="grey", prob=TRUE, main="Histogram: Fraction of Exome Fully Called", xlab="Fraction of Exome")
    hist(Genome_Partially_Called_Fraction, col="grey", prob=TRUE, main="Histogram: Fraction of Genome Partially Called", xlab="Fraction of Genome")
    hist(Exome_Partially_Called_Fraction,  col="grey", prob=TRUE, main="Histogram: Fraction of Exome Partially Called", xlab="Fraction of Exome")
    hist(Genome_No_Called_Fraction,        col="grey", prob=TRUE, main="Histogram: Fraction of Genome Not Called", xlab="Fraction of Genome")
    hist(Exome_No_Called_Fraction,         col="grey", prob=TRUE, main="Histogram: Fraction of Exome Not Called", xlab="Fraction of Exome")

}

funct_plot_genome_fractions_box <- function() {

    par(mfrow=c(3,2))
    boxplot(Genome_Fully_Called_Fraction,     main="Boxplot: Fraction of Genome Fully Called", xlab="Fraction of Genome", ylab="Fraction")
    boxplot(Exome_Fully_Called_Fraction,      main="Boxplot: Fraction of Exome Fully Called", xlab="Fraction of Exome", ylab="Fraction")
    boxplot(Genome_Partially_Called_Fraction, main="Boxplot: Fraction of Genome Partially Called", xlab="Fraction of Genome", ylab="Fraction")
    boxplot(Exome_Partially_Called_Fraction,  main="Boxplot: Fraction of Exome Partially Called", xlab="Fraction of Exome", ylab="Fraction")
    boxplot(Genome_No_Called_Fraction,        main="Boxplot: Fraction of Genome Not Called", xlab="Fraction of Genome", ylab="Fraction")
    boxplot(Exome_No_Called_Fraction,         main="Boxplot: Fraction of Exome Not Called", xlab="Fraction of Exome", ylab="Fraction")

}

funct_plot_genome_coverage <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.gac <- "Genome: >5X"
    lab.gbc <- "Genome: >10X"
    lab.gcc <- "Genome: >20X"
    lab.gdc <- "Genome: >30X"
    lab.gec <- "Genome: >40X"

    lab.eac <- "Exome: >5X"
    lab.ebc <- "Exome: >10X"
    lab.ecc <- "Exome: >20X"
    lab.edc <- "Exome: >30X"
    lab.eec <- "Exome: >40X"

    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:    
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.eec, lab.edc, lab.ecc, lab.ebc, lab.eac, lab.gec, lab.gdc, lab.gcc, lab.gbc, lab.gac)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
    
    plot<-xyplot(Exome_Fraction_WeightSumSequenceCoverage_gt_40x + Exome_Fraction_WeightSumSequenceCoverage_gt_30x + Exome_Fraction_WeightSumSequenceCoverage_gt_20x +
                 Exome_Fraction_WeightSumSequenceCoverage_gt_10x + Exome_Fraction_WeightSumSequenceCoverage_gt_5x +
                 Genome_Fraction_WeightSumSequenceCoverage_gt_40x + Genome_Fraction_WeightSumSequenceCoverage_gt_30x + Genome_Fraction_WeightSumSequenceCoverage_gt_20x +
                 Genome_Fraction_WeightSumSequenceCoverage_gt_10x + Genome_Fraction_WeightSumSequenceCoverage_gt_5x ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 10, 1), ylab="Fraction", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("Fraction of G/E with Coverage > X Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          ) 
    print(plot)
}

funct_plot_genome_coverage_merged <- function(xaxis_ticks) {

    par(mfrow=c(3,1))
    xgrid<-c(seq(0,max(File_Index)+5,by=5))

    ygrid<-c(seq(as.numeric(substr(min(Genome_Fraction_WeightSumSequenceCoverage_gt_40x),1,3))-0.1:1,by=0.1))
    plot(NULL,NULL,type="l",xlim=c(0,max(File_Index)),ylim=c(as.numeric(substr(min(Genome_Fraction_WeightSumSequenceCoverage_gt_40x)-0.1,1,3)),1.0), col="blue", main=paste("Genome X-fold Coverage Thresholds Across",length(File_Index),"Subjects"), xlab="Subject Index", ylab="Fraction of Genome", xaxt="n")
    axis(1, at=xaxis_ticks)
    abline(h = ygrid, v = xgrid, col = "grey90", lty=1, lwd=0.5)
    lines(File_Index,Genome_Fraction_WeightSumSequenceCoverage_gt_5x, col="red")
    lines(File_Index,Genome_Fraction_WeightSumSequenceCoverage_gt_10x, col="green")
    lines(File_Index,Genome_Fraction_WeightSumSequenceCoverage_gt_20x, col="darkgreen")
    lines(File_Index,Genome_Fraction_WeightSumSequenceCoverage_gt_30x, col="blue")
    lines(File_Index,Genome_Fraction_WeightSumSequenceCoverage_gt_40x, col="purple")
    legend("bottomleft",c(">5x",">10x",">20x",">30x",">40x"), cex=0.5, col=c("red", "green", "darkgreen", "blue", "purple"), lty=1, bty="n")

    ygrid<-c(seq(as.numeric(substr(min(Exome_Fraction_WeightSumSequenceCoverage_gt_40x),1,3))-0.1:1,by=0.1))
    plot(NULL,NULL,type="l",xlim=c(0,max(File_Index)),ylim=c(as.numeric(substr(min(Exome_Fraction_WeightSumSequenceCoverage_gt_40x)-0.1,1,3)),1.0), col="blue", main=paste("Exome X-fold Coverage Thresholds Across",length(File_Index),"Subjects"), xlab="Subject Index", ylab="Fraction of Exome", xaxt="n")
    axis(1, at=xaxis_ticks)
    abline(h = ygrid, v = xgrid, col = "grey90", lty=1, lwd=0.5)
    lines(File_Index,Exome_Fraction_WeightSumSequenceCoverage_gt_5x, col="red")
    lines(File_Index,Exome_Fraction_WeightSumSequenceCoverage_gt_10x, col="green")
    lines(File_Index,Exome_Fraction_WeightSumSequenceCoverage_gt_20x, col="darkgreen")
    lines(File_Index,Exome_Fraction_WeightSumSequenceCoverage_gt_30x, col="blue")
    lines(File_Index,Exome_Fraction_WeightSumSequenceCoverage_gt_40x, col="purple")
    legend("bottomleft",c(">5x",">10x",">20x",">30x",">40x"), cex=0.5, col=c("red", "green", "darkgreen", "blue", "purple"), lty=1, bty="n")

}

funct_plot_genome_coverage_hists <- function() {

    par(mfrow=c(5,2))
    hist(Genome_Fraction_WeightSumSequenceCoverage_gt_5x,  col="grey", prob=TRUE, main="Genome: Fraction WSSCov > 5X", xlab="Genome Fraction")
    hist(Exome_Fraction_WeightSumSequenceCoverage_gt_5x,   col="grey", prob=TRUE, main="Exome: Fraction WSSCov > 5X", xlab="Exome Fraction")
    hist(Genome_Fraction_WeightSumSequenceCoverage_gt_10x, col="grey", prob=TRUE, main="Genome: Fraction WSSCov > 10X", xlab="Genome Fraction")
    hist(Exome_Fraction_WeightSumSequenceCoverage_gt_10x,  col="grey", prob=TRUE, main="Exome: Fraction WSSCov > 10X", xlab="Exome Fraction")
    hist(Genome_Fraction_WeightSumSequenceCoverage_gt_20x, col="grey", prob=TRUE, main="Genome: Fraction WSSCov > 20X", xlab="Genome Fraction")
    hist(Exome_Fraction_WeightSumSequenceCoverage_gt_20x,  col="grey", prob=TRUE, main="Exome: Fraction WSSCov > 20X", xlab="Exome Fraction")
    hist(Genome_Fraction_WeightSumSequenceCoverage_gt_30x, col="grey", prob=TRUE, main="Genome: Fraction WSSCov > 30X", xlab="Genome Fraction")
    hist(Exome_Fraction_WeightSumSequenceCoverage_gt_30x,  col="grey", prob=TRUE, main="Exome: Fraction WSSCov > 30X", xlab="Exome Fraction")
    hist(Genome_Fraction_WeightSumSequenceCoverage_gt_40x, col="grey", prob=TRUE, main="Genome: Fraction WSSCov > 40X", xlab="Genome Fraction")
    hist(Exome_Fraction_WeightSumSequenceCoverage_gt_40x,  col="grey", prob=TRUE, main="Exome: Fraction WSSCov > 40X", xlab="Exome Fraction")

}

funct_plot_genome_coverage_box <- function() {

    par(mfrow=c(5,2))
    boxplot(Genome_Fraction_WeightSumSequenceCoverage_gt_5x,  main="Genome: Fraction WSSCov > 5X", xlab="Genome Fraction", ylab="Fraction")
    boxplot(Exome_Fraction_WeightSumSequenceCoverage_gt_5x,   main="Exome: Fraction WSSCov > 5X", xlab="Exome Fraction", ylab="Fraction")
    boxplot(Genome_Fraction_WeightSumSequenceCoverage_gt_10x, main="Genome: Fraction WSSCov > 10X", xlab="Genome Fraction", ylab="Fraction")
    boxplot(Exome_Fraction_WeightSumSequenceCoverage_gt_10x,  main="Exome: Fraction WSSCov > 10X", xlab="Exome Fraction", ylab="Fraction")
    boxplot(Genome_Fraction_WeightSumSequenceCoverage_gt_20x, main="Genome: Fraction WSSCov > 20X", xlab="Genome Fraction", ylab="Fraction")
    boxplot(Exome_Fraction_WeightSumSequenceCoverage_gt_20x,  main="Exome: Fraction WSSCov > 20X", xlab="Exome Fraction", ylab="Fraction")
    boxplot(Genome_Fraction_WeightSumSequenceCoverage_gt_30x, main="Genome: Fraction WSSCov > 30X", xlab="Genome Fraction", ylab="Fraction")
    boxplot(Exome_Fraction_WeightSumSequenceCoverage_gt_30x,  main="Exome: Fraction WSSCov > 30X", xlab="Exome Fraction", ylab="Fraction")
    boxplot(Genome_Fraction_WeightSumSequenceCoverage_gt_40x, main="Genome: Fraction WSSCov > 40X", xlab="Genome Fraction", ylab="Fraction")
    boxplot(Exome_Fraction_WeightSumSequenceCoverage_gt_40x,  main="Exome: Fraction WSSCov > 40X", xlab="Exome Fraction", ylab="Fraction")

}

funct_plot_genome_yields_cov_var <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.lab1 <- "Genome Gross Mapping Yield (Gb)"
    lab.lab2 <- "Genome Both Mates Mapped Yield (Gb)"
    lab.lab3 <- "Genome 100k Normalised Coverage Variability"
    
    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:   
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.lab3, lab.lab2, lab.lab1)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
        
    plot<-xyplot(Genome_100k_Normalised_Coverage_Variability + Genome_Both_Mates_Mapped_Yield_Gb + Genome_Gross_Mapping_Yield_Gb ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 3, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("Mapping Yield and Coverage Variability Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          )
    print(plot)
}

funct_plot_genome_yields_cov_var_hists <- function(files_indicies) {

    par(mfrow=c(3,2))
    hist(Genome_Gross_Mapping_Yield_Gb,               col="grey", prob=TRUE, main="Genome Gross Mapping Yield (Gb)",             xlab="Mapping Yield (Gb)")
    hist(Genome_Both_Mates_Mapped_Yield_Gb,           col="grey", prob=TRUE, main="Genome Both Mates Mapped Yield (Gb)",         xlab="Mapping Yield (Gb)")
    hist(Genome_100k_Normalised_Coverage_Variability, col="grey", prob=TRUE, main="Genome 100k Normalised Coverage Variability", xlab="Variability")

}

funct_plot_genome_yields_cov_var_box <- function() {

    par(mfrow=c(3,2))
    boxplot(Genome_Gross_Mapping_Yield_Gb,               main="Mapping Yield (Gb)", ylab="Yield (Gb)")
    boxplot(Genome_Both_Mates_Mapped_Yield_Gb,           main="Mapping Yield (Gb)", ylab="Yield (Gb)")
    boxplot(Genome_100k_Normalised_Coverage_Variability, main="Genome 100k Normalised Coverage Variability", ylab="Variability")

}

funct_plot_library_mate_data <- function(xaxis_ticks) {

    par(mfrow=c(3,1))
    xgrid<-c(seq(0,max(File_Index),by=5))
    ygrid<-c(seq(0,600,by=050))

    plot(NULL,NULL,type="l",xlim=c(0,max(File_Index)),ylim=c(0,600), col="blue", main=paste("Library Mate Mean with 95CI Min/Max Across",length(File_Index),"Subjects"), xlab="Subject Index", ylab="", xaxt="n")
    axis(1, at=xaxis_ticks)
    abline(h = ygrid, v = xgrid, col = "grey90", lty=1, lwd=0.5)
    lines(File_Index, Library_Mate_Distribution_Mean, col="red")
    lines(File_Index, Library_Mate_Distribution_Range_95CI_Min, col="green")
    lines(File_Index, Library_Mate_Distribution_Range_95CI_Max, col="blue")
    legend("bottomleft",c("Mate Distribution Mean","Mate Distribution 95CI Min","Mate Distribution 95CI Max"), cex=0.5, col=c("red", "green", "blue"), lty=1, bty="n")

}

funct_plot_genome_snp_stats <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.gasnp <- "Genome: SNP Total Count"
    lab.gbsnp <- "Genome: Homozygous Count"
    lab.gcsnp <- "Genome: Heterozygous Count"
    lab.gdsnp <- "Genome: Novel Fraction"
    lab.gesnp <- "Genome: Homozygous Novel Fraction"
    lab.gfsnp <- "Genome: Heterozygous Novel Fraction"
    lab.ggsnp <- "Genome: Heterozygous/Homozygous Ratio"
    lab.ghsnp <- "Genome: Transition/Transversion Ratio"
   
    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:   
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.ghsnp, lab.ggsnp, lab.gfsnp, lab.gesnp, lab.gdsnp, lab.gcsnp, lab.gbsnp, lab.gasnp)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
        
    plot<-xyplot(Genome_var_SNP_Transition_Transversion_Ratio + Genome_var_SNP_Het_Hom_Ratio + Genome_var_SNP_Heterozygous_Novel_Fraction +
                 Genome_var_SNP_Homozygous_Novel_Fraction + Genome_var_SNP_Novel_Fraction + Genome_var_SNP_Heterozygous_Count + Genome_var_SNP_Homozygous_Count +
                 Genome_var_SNP_Total_Count ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 8, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("Genome: SNP Data Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          )   
    print(plot)
}

funct_plot_genome_snp_stats_hists <- function() {

    par(mfrow=c(4,2))
    hist(Genome_var_SNP_Total_Count,                   col="grey", prob=TRUE, main="Genome: SNP Total Count",                   xlab="SNP Count")
    hist(Genome_var_SNP_Homozygous_Count,              col="grey", prob=TRUE, main="Genome: SNP Homozygous Count",              xlab="SNP Count")
    hist(Genome_var_SNP_Heterozygous_Count,            col="grey", prob=TRUE, main="Genome: SNP Heterozygous Count",            xlab="SNP Count")
    hist(Genome_var_SNP_Novel_Fraction,                col="grey", prob=TRUE, main="Genome: SNP Novel Fraction",                xlab="Novel Fraction")
    hist(Genome_var_SNP_Homozygous_Novel_Fraction,     col="grey", prob=TRUE, main="Genome: SNP Homozygous Novel Fraction",     xlab="Homozygous Novel Fraction")
    hist(Genome_var_SNP_Heterozygous_Novel_Fraction,   col="grey", prob=TRUE, main="Genome: SNP Heterozygous Novel Fraction",   xlab="Heterozygous Novel Fraction")
    hist(Genome_var_SNP_Het_Hom_Ratio,                 col="grey", prob=TRUE, main="Genome: SNP Het/Hom Ratio",                 xlab="Het/Hom Ratio")
    hist(Genome_var_SNP_Transition_Transversion_Ratio, col="grey", prob=TRUE, main="Genome: SNP Transition/Transversion Ratio", xlab="Transition/Transversion Ratio")
    
}

funct_plot_genome_snp_stats_box <- function() {

    par(mfrow=c(4,2))
    boxplot(Genome_var_SNP_Total_Count,                   main="Genome: SNP Total Count",                   ylab="SNP Count")
    boxplot(Genome_var_SNP_Homozygous_Count,              main="Genome: SNP Homozygous Count",              ylab="SNP Count")
    boxplot(Genome_var_SNP_Heterozygous_Count,            main="Genome: SNP Heterozygous Count",            ylab="SNP Count")
    boxplot(Genome_var_SNP_Novel_Fraction,                main="Genome: SNP Novel Fraction",                ylab="Novel Fraction")
    boxplot(Genome_var_SNP_Homozygous_Novel_Fraction,     main="Genome: SNP Homozygous Novel Fraction",     ylab="Homozygous Novel Fraction")
    boxplot(Genome_var_SNP_Heterozygous_Novel_Fraction,   main="Genome: SNP Heterozygous Novel Fraction",   ylab="Heterozygous Novel Fraction")
    boxplot(Genome_var_SNP_Het_Hom_Ratio,                 main="Genome: SNP Het/Hom Ratio",                 ylab="Het/Hom Ratio")
    boxplot(Genome_var_SNP_Transition_Transversion_Ratio, main="Genome: SNP Transition/Transversion Ratio", ylab="Transition/Transversion Ratio")

}

funct_plot_exome_snp_stats <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.easnp <- "Exome: SNP Total Count"
    lab.ebsnp <- "Exome: Homozygous Count"
    lab.ecsnp <- "Exome: Heterozygous Count"
    lab.edsnp <- "Exome: Novel Fraction"
    lab.eesnp <- "Exome: Homozygous Novel Fraction"
    lab.efsnp <- "Exome: Heterozygous Novel Fraction"
    lab.egsnp <- "Exome: Heterozygous/Homozygous Ratio"
    lab.ehsnp <- "Exome: Transition/Transversion Ratio"

    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:    
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.ehsnp, lab.egsnp, lab.efsnp, lab.eesnp, lab.edsnp, lab.ecsnp, lab.ebsnp, lab.easnp)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
        
    plot<-xyplot(Exome_var_SNP_Transition_Transversion_Ratio + Exome_var_SNP_Het_Hom_Ratio + Exome_var_SNP_Heterozygous_Novel_Fraction +
                 Exome_var_SNP_Homozygous_Novel_Fraction + Exome_var_SNP_Novel_Fraction + Exome_var_SNP_Heterozygous_Count + Exome_var_SNP_Homozygous_Count +
                 Exome_var_SNP_Total_Count ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 8, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("Exome: SNP Data Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          )  
    print(plot)
}

funct_plot_exome_snp_stats_hists <- function() {

    par(mfrow=c(4,2))
    hist(Exome_var_SNP_Total_Count,                   col="grey", prob=TRUE, main="Exome: SNP Total Count",                   xlab="SNP Count")
    hist(Exome_var_SNP_Homozygous_Count,              col="grey", prob=TRUE, main="Exome: SNP Homozygous Count",              xlab="SNP Count")
    hist(Exome_var_SNP_Heterozygous_Count,            col="grey", prob=TRUE, main="Exome: SNP Heterozygous Count",            xlab="SNP Count")
    hist(Exome_var_SNP_Novel_Fraction,                col="grey", prob=TRUE, main="Exome: SNP Novel Fraction",                xlab="Novel Fraction")
    hist(Exome_var_SNP_Homozygous_Novel_Fraction,     col="grey", prob=TRUE, main="Exome: SNP Homozygous Novel Fraction",     xlab="Homozygous Novel Fraction")
    hist(Exome_var_SNP_Heterozygous_Novel_Fraction,   col="grey", prob=TRUE, main="Exome: SNP Heterozygous Novel Fraction",   xlab="Heterozygous Novel Fraction")
    hist(Exome_var_SNP_Het_Hom_Ratio,                 col="grey", prob=TRUE, main="Exome: SNP Het/Hom Ratio",                 xlab="Het/Hom Ratio")
    hist(Exome_var_SNP_Transition_Transversion_Ratio, col="grey", prob=TRUE, main="Exome: SNP Transition/Transversion Ratio", xlab="Transition/Transversion Ratio")

}

funct_plot_exome_snp_stats_box <- function() {

    par(mfrow=c(4,2))
    boxplot(Exome_var_SNP_Total_Count,                   main="Exome: SNP Total Count",                   ylab="SNP Count")
    boxplot(Exome_var_SNP_Homozygous_Count,              main="Exome: SNP Homozygous Count",              ylab="SNP Count")
    boxplot(Exome_var_SNP_Heterozygous_Count,            main="Exome: SNP Heterozygous Count",            ylab="SNP Count")
    boxplot(Exome_var_SNP_Novel_Fraction,                main="Exome: SNP Novel Fraction",                ylab="Novel Fraction")
    boxplot(Exome_var_SNP_Homozygous_Novel_Fraction,     main="Exome: SNP Homozygous Novel Fraction",     ylab="Homozygous Novel Fraction")
    boxplot(Exome_var_SNP_Heterozygous_Novel_Fraction,   main="Exome: SNP Heterozygous Novel Fraction",   ylab="Heterozygous Novel Fraction")
    boxplot(Exome_var_SNP_Het_Hom_Ratio,                 main="Exome: SNP Het/Hom Ratio",                 ylab="Het/Hom Ratio")
    boxplot(Exome_var_SNP_Transition_Transversion_Ratio, main="Exome: SNP Transition/Transversion Ratio", ylab="Transition/Transversion Ratio")

}


funct_plot_ins_stats <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.gains <- "Genome: INS Total Count"
    lab.gbins <- "Genome: INS Novel Fraction"
    lab.gcins <- "Genome: INS Heterozygosity/Homozygosity Ratio"
    lab.eains <- "Exome: INS Total Count"
    lab.ebins <- "Exome: INS Novel Fraction"
    lab.ecins <- "Exome: INS Heterozygosity/Homozygosity Ratio"

    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:    
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.ecins, lab.ebins, lab.eains, lab.gcins, lab.gbins, lab.gains)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
    
    plot<-xyplot(Exome_var_INS_Het_Hom_Ratio  + Exome_var_INS_Novel_Fraction  + Exome_var_INS_Total_Count + 
                 Genome_var_INS_Het_Hom_Ratio + Genome_var_INS_Novel_Fraction + Genome_var_INS_Total_Count ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 6, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("Genome/Exome INS Data Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          )
    print(plot)
}

funct_plot_ins_stats_hists <- function() {

    par(mfrow=c(3,2))
    hist(Genome_var_INS_Total_Count,    col="grey", prob=TRUE, main="Genome: INS Total Count",    xlab="INS Count")
    hist(Exome_var_INS_Total_Count,     col="grey", prob=TRUE, main="Exome: INS Total Count",     xlab="INS Count")
    hist(Genome_var_INS_Novel_Fraction, col="grey", prob=TRUE, main="Genome: INS Novel Fraction", xlab="Novel Fraction")
    hist(Exome_var_INS_Novel_Fraction,  col="grey", prob=TRUE, main="Exome: INS Novel Fraction",  xlab="Novel Fraction")
    hist(Genome_var_INS_Het_Hom_Ratio,  col="grey", prob=TRUE, main="Genome: INS Het/Hom Ratio",  xlab="Het/Hom Ratio")
    hist(Exome_var_INS_Het_Hom_Ratio,   col="grey", prob=TRUE, main="Exome: INS Het/Hom Ratio",   xlab="Het/Hom Ratio")
    
}

funct_plot_ins_stats_box <- function() {

    par(mfrow=c(3,2))
    boxplot(Genome_var_INS_Total_Count,    main="Genome: INS Total Count",    ylab="INS Count")
    boxplot(Exome_var_INS_Total_Count,     main="Exome: INS Total Count",     ylab="INS Count")
    boxplot(Genome_var_INS_Novel_Fraction, main="Genome: INS Novel Fraction", ylab="Novel Fraction")
    boxplot(Exome_var_INS_Novel_Fraction,  main="Exome: INS Novel Fraction",  ylab="Novel Fraction")
    boxplot(Genome_var_INS_Het_Hom_Ratio,  main="Genome: INS Het/Hom Ratio",  ylab="Het/Hom Ratio")
    boxplot(Exome_var_INS_Het_Hom_Ratio,   main="Exome: INS Het/Hom Ratio",   ylab="Het/Hom Ratio")
    
}

funct_plot_del_stats <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.gadel <- "Genome: DEL Total Count"
    lab.gbdel <- "Genome: DEL Novel Fraction"
    lab.gcdel <- "Genome: DEL Heterozygosity/Homozygosity Ratio"
    lab.eadel <- "Exome: DEL Total Count"
    lab.ebdel <- "Exome: DEL Novel Fraction"
    lab.ecdel <- "Exome: DEL Heterozygosity/Homozygosity Ratio"

    
    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:   
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.ecdel, lab.ebdel, lab.eadel, lab.gcdel, lab.gbdel, lab.gadel)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
        
    plot<-xyplot(Exome_var_DEL_Het_Hom_Ratio  + Exome_var_DEL_Novel_Fraction  + Exome_var_DEL_Total_Count + 
                 Genome_var_DEL_Het_Hom_Ratio + Genome_var_DEL_Novel_Fraction + Genome_var_DEL_Total_Count ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 6, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("Genome/Exome DEL Data Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          )
    print(plot)
}

funct_plot_del_stats_hists <- function() {

    par(mfrow=c(3,2))
    hist(Genome_var_DEL_Total_Count,    col="grey", prob=TRUE, main="Genome: DEL Total Count",    xlab="DEL Count")
    hist(Exome_var_DEL_Total_Count,     col="grey", prob=TRUE, main="Exome: DEL Total Count",     xlab="DEL Count")
    hist(Genome_var_DEL_Novel_Fraction, col="grey", prob=TRUE, main="Genome: DEL Novel Fraction", xlab="Novel Fraction")
    hist(Exome_var_DEL_Novel_Fraction,  col="grey", prob=TRUE, main="Exome: DEL Novel Fraction",  xlab="Novel Fraction")
    hist(Genome_var_DEL_Het_Hom_Ratio,  col="grey", prob=TRUE, main="Genome: DEL Het/Hom Ratio",  xlab="Het/Hom Ratio")
    hist(Exome_var_DEL_Het_Hom_Ratio,   col="grey", prob=TRUE, main="Exome: DEL Het/Hom Ratio",   xlab="Het/Hom Ratio")
    
}

funct_plot_del_stats_box <- function() {

    par(mfrow=c(3,2))
    boxplot(Genome_var_DEL_Total_Count,    main="Genome: DEL Total Count",    ylab="DEL Count")
    boxplot(Exome_var_DEL_Total_Count,     main="Exome: DEL Total Count",     ylab="DEL Count")
    boxplot(Genome_var_DEL_Novel_Fraction, main="Genome: DEL Novel Fraction", ylab="Novel Fraction")
    boxplot(Exome_var_DEL_Novel_Fraction,  main="Exome: DEL Novel Fraction",  ylab="Novel Fraction")
    boxplot(Genome_var_DEL_Het_Hom_Ratio,  main="Genome: DEL Het/Hom Ratio",  ylab="Het/Hom Ratio")
    boxplot(Exome_var_DEL_Het_Hom_Ratio,   main="Exome: DEL Het/Hom Ratio",   ylab="Het/Hom Ratio")
    
}


funct_plot_sub_stats <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.gasub <- "Genome: SUB Total Count"
    lab.gbsub <- "Genome: SUB Novel Fraction"
    lab.gcsub <- "Genome: SUB Heterozygosity/Homozygosity Ratio"
    lab.easub <- "Exome: SUB Total Count"
    lab.ebsub <- "Exome: SUB Novel Fraction"
    lab.ecsub <- "Exome: SUB Heterozygosity/Homozygosity Ratio"

    
    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:    
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.ecsub, lab.ebsub, lab.easub, lab.gcsub, lab.gbsub, lab.gasub)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
        
    plot<-xyplot(Exome_var_SUB_Het_Hom_Ratio  + Exome_var_SUB_Novel_Fraction  + Exome_var_SUB_Total_Count + 
                 Genome_var_SUB_Het_Hom_Ratio + Genome_var_SUB_Novel_Fraction + Genome_var_SUB_Total_Count ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 6, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("Genome/Exome SUB Data Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          )
    print(plot)
}

funct_plot_sub_stats_hists <- function() {

    par(mfrow=c(3,2))
    hist(Genome_var_SUB_Total_Count,    col="grey", prob=TRUE, main="Genome: SUB Total Count",    xlab="SUB Count")
    hist(Exome_var_SUB_Total_Count,     col="grey", prob=TRUE, main="Exome: SUB Total Count",     xlab="SUB Count")
    hist(Genome_var_SUB_Novel_Fraction, col="grey", prob=TRUE, main="Genome: SUB Novel Fraction", xlab="Novel Fraction")
    hist(Exome_var_SUB_Novel_Fraction,  col="grey", prob=TRUE, main="Exome: SUB Novel Fraction",  xlab="Novel Fraction")
    hist(Genome_var_SUB_Het_Hom_Ratio,  col="grey", prob=TRUE, main="Genome: SUB Het/Hom Ratio",  xlab="Het/Hom Ratio")
    hist(Exome_var_SUB_Het_Hom_Ratio,   col="grey", prob=TRUE, main="Exome: SUB Het/Hom Ratio",   xlab="Het/Hom Ratio")
    
}

funct_plot_sub_stats_box <- function() {

    par(mfrow=c(3,2))
    boxplot(Genome_var_SUB_Total_Count,    main="Genome: SUB Total Count",    ylab="SUB Count")
    boxplot(Exome_var_SUB_Total_Count,     main="Exome: SUB Total Count",     ylab="SUB Count")
    boxplot(Genome_var_SUB_Novel_Fraction, main="Genome: SUB Novel Fraction", ylab="Novel Fraction")
    boxplot(Exome_var_SUB_Novel_Fraction,  main="Exome: SUB Novel Fraction",  ylab="Novel Fraction")
    boxplot(Genome_var_SUB_Het_Hom_Ratio,  main="Genome: SUB Het/Hom Ratio",  ylab="Het/Hom Ratio")
    boxplot(Exome_var_SUB_Het_Hom_Ratio,   main="Exome: SUB Het/Hom Ratio",   ylab="Het/Hom Ratio")
    
}

funct_plot_func_impact_snps <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.fisnpa <- "SNPs: Synonymous Counts"
    lab.fisnpb <- "SNPs: Non-Synonymous"
    lab.fisnpc <- "SNPs: Missense"
    lab.fisnpd <- "SNPs: Nonsense"
    lab.fisnpe <- "SNPs: Nonstop"
    lab.fisnpf <- "SNPs: Misstart"
    lab.fisnpg <- "SNPs: Disrupt"

    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:   
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.fisnpg, lab.fisnpf, lab.fisnpe, lab.fisnpd, lab.fisnpc, lab.fisnpb, lab.fisnpa)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
        
    plot<-xyplot(Functional_impact_Disrupt_SNP_Loci  + Functional_impact_Misstart_SNP_Loci       + Functional_impact_Nonstop_SNP_Loci    + Functional_impact_Nonsense_SNP_Loci +
                 Functional_impact_Missense_SNP_Loci + Functional_impact_Non_synonymous_SNP_Loci + Functional_impact_Synonymous_SNP_Loci ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 6, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("SNP Functional Impact Counts Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          )
    print(plot)
}

funct_plot_func_impact_snps_hists <- function() {
    
    par(mfrow=c(4,2))
    hist(Functional_impact_Synonymous_SNP_Loci,     col="grey", prob=TRUE, main="Synonymous SNPs",     xlab="SNP Count")
    hist(Functional_impact_Non_synonymous_SNP_Loci, col="grey", prob=TRUE, main="Non-synonymous SNPs", xlab="SNP Count")
    hist(Functional_impact_Missense_SNP_Loci,       col="grey", prob=TRUE, main="Missense SNPs",       xlab="SNP Count")
    hist(Functional_impact_Nonsense_SNP_Loci,       col="grey", prob=TRUE, main="Nonsense SNPs",       xlab="SNP Count")
    hist(Functional_impact_Nonstop_SNP_Loci,        col="grey", prob=TRUE, main="Nonstop SNPs",        xlab="SNP Count")
    hist(Functional_impact_Misstart_SNP_Loci,       col="grey", prob=TRUE, main="Misstart SNPs",       xlab="SNP Count")
    hist(Functional_impact_Disrupt_SNP_Loci,        col="grey", prob=TRUE, main="Disrupt SNPs",        xlab="SNP Count")

}

funct_plot_func_impact_snps_box <- function() {
    
    par(mfrow=c(4,2))
    boxplot(Functional_impact_Synonymous_SNP_Loci,     main="Synonymous SNPs",     ylab="SNP Count")
    boxplot(Functional_impact_Non_synonymous_SNP_Loci, main="Non-synonymous SNPs", ylab="SNP Count")
    boxplot(Functional_impact_Missense_SNP_Loci,       main="Missense SNPs",       ylab="SNP Count")
    boxplot(Functional_impact_Nonsense_SNP_Loci,       main="Nonsense SNPs",       ylab="SNP Count")
    boxplot(Functional_impact_Nonstop_SNP_Loci,        main="Nonstop SNPs",        ylab="SNP Count")
    boxplot(Functional_impact_Misstart_SNP_Loci,       main="Misstart SNPs",       ylab="SNP Count")
    boxplot(Functional_impact_Disrupt_SNP_Loci,        main="Disrupt SNPs",        ylab="SNP Count")

}

funct_plot_func_impact_framing <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.fiframea <- "Frame-Shifting: INS"
    lab.fiframeb <- "Frame-Shifting: DEL"
    lab.fiframec <- "Frame-Shifting: SUB"
    lab.fiframed <- "Frame-Preserving: INS"
    lab.fiframee <- "Frame-Preserving: DEL"
    lab.fiframef <- "Frame-Preserving: SUB"
 
    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.fiframef, lab.fiframee, lab.fiframed, lab.fiframec, lab.fiframeb, lab.fiframea)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
    
    plot<-xyplot(Functional_impact_Frame_Preserving_SUB_Loci + Functional_impact_Frame_Preserving_DEL_Loci + Functional_impact_Frame_Preserving_INS_Loci + 
                 Functional_impact_Frame_shifting_SUB_Loci   + Functional_impact_Frame_Shifting_DEL_Loci   + Functional_impact_Frame_Shifting_INS_Loci ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 6, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("Frame-* Functional Impact Counts Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          )
    print(plot)
}

funct_plot_func_impact_framing_hists <- function() {

    par(mfrow=c(3,2))
    hist(Functional_impact_Frame_Shifting_INS_Loci,   col="grey", prob=TRUE, main="INS: Frame-shifting",   xlab="INS Count")
    hist(Functional_impact_Frame_Preserving_INS_Loci, col="grey", prob=TRUE, main="INS: Frame-preserving", xlab="INS Count")
    hist(Functional_impact_Frame_Shifting_DEL_Loci,   col="grey", prob=TRUE, main="DEL: Frame-shifting",   xlab="DEL Count")
    hist(Functional_impact_Frame_Preserving_DEL_Loci, col="grey", prob=TRUE, main="DEL: Frame-preserving", xlab="DEL Count")
    hist(Functional_impact_Frame_shifting_SUB_Loci,   col="grey", prob=TRUE, main="SUB: Frame-shifting",   xlab="SUB Count")
    hist(Functional_impact_Frame_Preserving_SUB_Loci, col="grey", prob=TRUE, main="SUB: Frame-preserving", xlab="SUB Count")

}

funct_plot_func_impact_framing_box <- function() {

    par(mfrow=c(3,2))
    boxplot(Functional_impact_Frame_Shifting_INS_Loci,   main="INS: Frame-shifting",   ylab="INS Count")
    boxplot(Functional_impact_Frame_Preserving_INS_Loci, main="INS: Frame-preserving", ylab="INS Count")
    boxplot(Functional_impact_Frame_Shifting_DEL_Loci,   main="DEL: Frame-shifting",   ylab="DEL Count")
    boxplot(Functional_impact_Frame_Preserving_DEL_Loci, main="DEL: Frame-preserving", ylab="DEL Count")
    boxplot(Functional_impact_Frame_shifting_SUB_Loci,   main="SUB: Frame-shifting",   ylab="SUB Count")
    boxplot(Functional_impact_Frame_Preserving_SUB_Loci, main="SUB: Frame-preserving", ylab="SUB Count")

}

funct_plot_cnvs <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.cnva <- "CNV: Total Count"
    lab.cnvb <- "CNV: Total Bases in CNV Segments"
    lab.cnvc <- "CNV: Fraction Novel (by CNV count)"
    lab.cnvd <- "CNV: Fraction Novel (by base count)"

 
    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.cnvd, lab.cnvc, lab.cnvb, lab.cnva)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
    
    plot<-xyplot(CNV_Fraction_Novel_by_Base_Count + CNV_Fraction_Novel_by_Segment_Count + CNV_Total_Bases_in_CNV_Segments + CNV_Total_CNV_Segment_Count ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 4, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("CNV Data Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          )
    print(plot)
}

funct_plot_cnvs_hists <- function() {
    
    par(mfrow = c(2,2))
    hist(CNV_Total_CNV_Segment_Count,         col="grey", prob=TRUE, main="CNV Total Segment Count",             xlab="Count")
    hist(CNV_Total_Bases_in_CNV_Segments,     col="grey", prob=TRUE, main="CNV Total Bases in CNV Segments",     xlab="Count")
    hist(CNV_Fraction_Novel_by_Segment_Count, col="grey", prob=TRUE, main="CNV Fraction Novel by Segment Count", xlab="Fraction")
    hist(CNV_Fraction_Novel_by_Base_Count,    col="grey", prob=TRUE, main="CNV Fraction Novel by Base Count",    xlab="Fraction")

}

funct_plot_cnvs_box <- function() {
    
    par(mfrow = c(2,2))
    boxplot(CNV_Total_CNV_Segment_Count,         main="CNV Total Segment Count",             ylab="Count")
    boxplot(CNV_Total_Bases_in_CNV_Segments,     main="CNV Total Bases in CNV Segments",     ylab="Count")
    boxplot(CNV_Fraction_Novel_by_Segment_Count, main="CNV Fraction Novel by Segment Count", ylab="Fraction")
    boxplot(CNV_Fraction_Novel_by_Base_Count,    main="CNV Fraction Novel by Base Count",    ylab="Fraction")

}

funct_plot_sv <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.sva <- "SV: Total Junction Count"
    lab.svb <- "SV: High Confidence Junction Count"

    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.svb, lab.sva)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
        
    plot<-xyplot(SV_High_Confidence_Junction_Count + SV_Total_Junction_Count ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 2, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("SV Data Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          ) 
    print(plot)
}

funct_plot_sv_hists <- function() {

    par(mfrow = c(3,2))
    hist(SV_Total_Junction_Count,           col="grey", prob=TRUE, main="SV Total Junction Count",            xlab="Count")
    hist(SV_High_Confidence_Junction_Count, col="grey", prob=TRUE, main="SV High Confidence Junction Count",  xlab="Count")
}

funct_plot_sv_box <- function() {

    par(mfrow = c(3,2))
    boxplot(SV_Total_Junction_Count,           main="SV Total Junction Count",            ylab="Count")
    boxplot(SV_High_Confidence_Junction_Count, main="SV High Confidence Junction Count",  ylab="Count")

}

funct_find_corr <- function(summary_dataframe){
	plot.num = 0
	for(i in 4:(ncol(summary_dataframe)-1)){
		for(j in (i+1):ncol(summary_dataframe)){
			#test.cor<-cor.test(summary_dataframe[,i],summary_dataframe[,j], method="pearson")
			#test.cor<-summary(lm(summary_dataframe[,i] ~ summary_dataframe[,j]))
			#print(summary.lm(lm(summary_dataframe[,i] ~ summary_dataframe[,j])))
			lm.sum<-lm(summary_dataframe[,j] ~ summary_dataframe[,i])
			p.value<-summary.lm(lm.sum)$coefficients[2,4]
			#print(summary.lm(lm(summary_dataframe[,j] ~ summary_dataframe[,i])))
			#print(p.value)
			#p.value=1
			#print(names(test.cor))
			#if(!is.na(test.cor$p.value)){
				if(p.value <= .05){
					if((plot.num %% 6) == 0){
						par(mfrow = c(3,2))
					}
					plot(summary_dataframe[,i],summary_dataframe[,j], xlab=names(summary_dataframe)[i], ylab=names(summary_dataframe)[j])
					abline(lm.sum)
					if(p.value<2e-16){
						mtext(paste("Regression pvalue = <2e-16", sep=""))
					}else{
						mtext(paste("Regression pvalue = ", round(p.value,5), sep=""))
					} 
					plot.num = plot.num+1
				}
			#}		
		}	
	}
}

funct_plot_mei <- function(xlim,xaxis,d) {

    #--Define plot titles:
    lab.meia <- "MEI: Mobile element insertion count"
    lab.meib <- "MEI: Fraction of novel MEI"

    #--Define colours for raw data & median:
    col.raw <- "#377EB8"
    col.lm <- "grey20"
    
    #--Custom strip function:
    my.strip <- function(which.given, which.panel, ...) {
        strip.labels <- c(lab.meia, lab.meib)
        panel.rect(0, 0, 1, 1, col="Grey", border=1)
        panel.text(x=0.5, y=0.5, adj=c(0.5, 0.55), cex=0.95,lab=strip.labels[which.panel[which.given]],col="Black")
    }
        
    plot<-xyplot(SV_High_Confidence_Junction_Count + SV_Total_Junction_Count ~ File_Index,
                 scales=list(y="free", x=xaxis, rot=0), xlim=xlim,
                 strip=my.strip, outer=TRUE, layout=c(1, 2, 1), ylab="", xlab="Subject Index",
                 panel=function(x, y, ...) {
                     panel.grid(h=-1, v=0)	# plot default horizontal gridlines
                     panel.abline(v=d, col="grey90") # custom vertical gridlines
                     panel.xyplot(x, y, ..., type="l", col=col.raw, lwd=0.5) # raw data
                     panel.abline(h=median(y, na.rm=TRUE), lty=2, col=col.lm, lwd=1) # median value
                 },
                 key=list(text=list(c("Raw Data", "Median")),
                          title=paste("MEI Data Across",length(File_Index),"Subjects"),
                          col=c(col.raw, col.lm), lty=c(1, 2),
                          columns=2, cex=0.95,
                          lines=TRUE
                 ),
          ) 
    print(plot)
}

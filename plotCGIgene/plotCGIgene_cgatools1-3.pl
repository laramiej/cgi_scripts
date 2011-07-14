#!/usr/bin/perl -w

use strict;
use Getopt::Std;
my %options=();

getopts("i:s:g:n:d:r:x:", \%options);

#Jason M Laramie
#jlaramie@completegenomics.com
#program needs to be run with plotCGIgene.R.  Will create a pretty plot for a gene and a number of CGI data sources.
#will only work with Mastervars generated by cgatools 1.3

#damage rankings
my %func = ("NO-CHANGE" => 0, "COMPATIBLE" => 0, "MISSENSE"=>4,
			"NONSENSE"=>11, "NONSTOP"=>12, "DELETE"=>5,
			"DELETE+"=>7, "INSERT"=>6, "INSERT+"=>8,
			"FRAMESHIFT"=>13,"MISSTART"=>10,"DISRUPT"=>9,
			"UNKNOWN-VNC"=>1,"UNKNOWN-INC"=>3,"UNKNOWN-TR"=>2);
#componet rankings
my %comp = ("CDS"=>5,"INTRON"=>1,"DONOR"=>4,"ACCEPTOR"=>4,"TSS-UPSTREAM"=>0,"SPAN5"=>3,"SPAN3"=>3,"SPAN"=>3,"UTR5"=>2,"UTR3"=>2,"UTR"=>2);

#check for input
if (!defined $options{i}){
	if(!defined $options{s}){
		usage("option -s or option -i need to be defined\n");
	}
}else{
	$options{i} .= "/ASM";
	if (!(-d $options{i})){
		usage("Input directory $options{i} doesn't exist\n");
	}
}
if(!defined $options{g}){
	usage("option -g not defined\n");
}elsif(!defined $options{d}){
	usage("option -d not defined\n");
	if(!(-d $options{d})){
		usage("output directory $options{d} doesn't exist\n");
	}
}elsif(!defined $options{r}){
	usage("option -r not defined\n");
}elsif(!defined $options{x}){
	$options{x}=0;
}elsif(!defined $options{n}){
	$options{n}=999999999999999999;
}
#need to set these according to your personal setup
#Location of UCSC tracks for gene locations
my $pathToKnownGene = "/home/jlaramie/refs/$options{r}/knownGene.txt";
my $pathTokgXref = "/home/jlaramie/refs/$options{r}/kgXref.txt";
my $pathToKnown2locusId = "/home/jlaramie/refs/$options{r}/knownToLocusLink.txt";
my $pathTorefLink = "/home/jlaramie/refs/$options{r}/refLink.txt";

#path to plotCGIgene.R
my $pathToPlotter = "/home/jlaramie/biogen/scripts/plotCGIgene.R";

#Prefixes of specific files
my $masterVarPrefix = "^mastervar";
my $cnvdetailsPrefix = "^cnvDetailsBeta";


my $pathToASM = $options{i};
my $geneList = $options{g};
my %geneLookUp;

	
sub getCoverage ($$$$){
	my ($chr, $start, $end, $geneName) = @_;
	my $sum;
	
	my $outRefScore = $geneName . "-refscore.txt";	
	my $OUT_FH = open_file_for_writing($options{d}, $outRefScore, 1);
	my ($IN_FH, $doc_file) = open_file_for_reading("$pathToASM/REF/", qr{^coverageRefScore-$chr-.*});
	
	my $lineCount = 0;
	my $ntCount = 0;
	while(<$IN_FH>){
		chomp;
		$lineCount++;
		if($lineCount <= 10){
			next;
		}
		if($lineCount % 10000000 == 0){
			print "Working on $lineCount line\n";
		}	
		my @line = split("\t", $_);
		if($line[0] <= $end && $line[0] >= $start){
			print $OUT_FH "$_\n";
			$ntCount++;
			$sum += $line[5];
		}
		if($line[0] >= $end){
			#print "End of gene region. Exiting\n";
			last;
		}
	}
	my $avg = int($sum/$ntCount);
	print "$geneName has the average gross coverage of $avg" . "X\n";
	close($OUT_FH);
	close($IN_FH);
}

sub calcPloidy($$$$){
	my ($chr,$start,$end, $geneName) = @_;
	
	my ($IN_FH, $doc_file) = open_file_for_reading("$pathToASM/CNV/", qr{^$cnvdetailsPrefix-.*});
	my $outCNV = $geneName . "-cnv.txt";
	my $OUT_FH = open_file_for_writing($options{d}, $outCNV, 1);
	
	while(<$IN_FH>){
		last if(/^>/);
	}
	while(<$IN_FH>){
		chomp;
		my @line=split("\t",$_);
		if(($line[0] eq $chr) && ($line[1] <= $end && $line[1] >= $start)){
			my $segment = $line[1]+2000;
			if($line[7] eq "hypervariable"){
				print $OUT_FH "$line[1]\t$segment\t$line[5]\tH\n";
			}elsif($line[7] eq "invariant"){
				print $OUT_FH "$line[1]\t$segment\t$line[5]\tI\n";
			}else{
				print $OUT_FH "$line[1]\t$segment\t$line[5]\t$line[6]\n";
			}
		}
	}
	close($OUT_FH);
	close($IN_FH);
}

sub calcAllelicImbalance ($$$$){
	my ($chr,$start,$end, $geneName) = @_;
	
	my ($IN_FH, $doc_file) = open_file_for_reading("$pathToASM/", qr{^$masterVarPrefix-.*});

	my $outAlleleBalance = $geneName . "-alleleBalance.txt";
	my $outVariants = $geneName . "-alleleVariants.txt";
	
	my $OUT_FH = open_file_for_writing($options{d}, $outAlleleBalance, 1);
	my $VAR_FH = open_file_for_writing($options{d}, $outVariants, 1);
	
	my $lineCount = 0;
	my $noCallBases = 0;
	my $noCallPercent;
	my $halfCallCount = 0;
	while(<$IN_FH>){
		chomp;
		$lineCount++;
		if($lineCount <= 20){
			next;
		}
		if($lineCount % 10000000 == 0){
			print "Working on $lineCount line\n";
		}
		#look for het calls
		my @line = split("\t", $_);
		
		#Am I on the right chromosome and within the start & stop of the canonical gene?
		if(($line[2] eq $chr) && ($line[3] <= $end && $line[3] >= $start)){
			#in the gene Area
			my $count=0;
			#foreach my $elem (@line){
			#	print "$count\t$elem\n";
			#	$count++
			#}
			#die;
			if($line[6] eq "snp" && $line[5] eq "het-ref"){
				print $OUT_FH "$line[2]\t$line[3]\t$line[6]\t$line[16]\t$line[17]\t";
				my $ratio;
				if($line[16] < $line[17]){
					$ratio = $line[16]/($line[16]+$line[17]);
				}else{
					$ratio = $line[17]/($line[17]+$line[16]);
				}
				print $OUT_FH "$ratio\n";
			}
			#find variation 
			if($line[5] ne "no-call" && $line[5] ne "half"){
				if($line[6] ne "ref"){
					#print "gene: $line[21]\n";
					my @transFunction = split("\;",$line[21]);
					#print "@transFunction\n";
					my $highImpact = -1;
					my $highComp = -1;
					my $bestComp = "missing";
					my $bestFunc= "missing";
					foreach my $transFunc (@transFunction){
						my ($geneId,$mrnaAcc,$symbol,$pfam,$pfam2,$pfam3,$component,$impact) = split("\:", $transFunc);
						#print "$symbol\t$impact\n";
						if($symbol eq $geneName){
							#add rules for impact
							if($highComp<=$comp{$component}){
								$highComp = $comp{$component};
								$bestComp = $component;
							}
							if($highImpact<=$func{$impact}){
								$highImpact = $func{$impact};
								$bestFunc = $impact;
							}
						}
					}
					my $numbp = 0;
					if($line[6] eq "snp"){
						$numbp=1;
					}elsif($line[6] eq "ins"){
						$numbp = length($line[8]);
					}elsif($line[6] eq "del"){
						$numbp = length($line[7]);
					}elsif($line[6] eq "sub"){
						$numbp = $line[4]-$line[3];
					}else{
						$numbp="NA";
					}
					print $VAR_FH "$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\t$line[11]\t$bestComp\t$bestFunc\t$numbp\n";
				}
			}elsif($line[5] eq "no-call"){
				$noCallBases += ($line[4] - $line[3]);
				$noCallPercent = int(($noCallBases / ($end-$start))*100);
			}elsif($line[5] eq "half"){
				$halfCallCount++;	
			}
		}
		if(($line[2] eq $chr) && ($line[3] > $end)){
			#past where we need to be
			last;
		}	 	
	}
	print "No-call percent: $noCallPercent" . "% with $noCallBases bases not called\n";
	print "The number of half-called alleles: $halfCallCount\n";
	close($IN_FH);
	close($VAR_FH);
	close($OUT_FH);
	return($noCallPercent,$halfCallCount);
}	

sub open_file_for_reading($$){
    my ($dir, $fmask) = @_;
    opendir DIR, $dir or die "Cannot open directory $dir: $!, stopped ";
    my @files = grep { /$fmask/ } readdir DIR;
    (scalar @files == 1)
        or die (((scalar @files == 0) ? 'No' : 'Multiple') .
            " files matching file mask $fmask exist in directory $dir\n" .
            join("\n", @files) .
            ", stopped ");
    my $fname = $files[0];
    print "Reading $fname...\n";
    my $fhandle;   
    if ($fname=~/\.gz$/) { # if it ends with .gz
    	open ($fhandle,"zcat $dir/$fname |")  || die "can't open 'zcat $dir/$fname'";
	} elsif ($fname=~/\.bz2$/) { # if it ends with .bz2
    	open ($fhandle,"bzcat $dir/$fname |") || die "can't open 'bzcat $dir/$fname'";
	} else {
    	open ($fhandle,"<$dir/$fname") || die "can't open '$dir/$fname'";
	}
    return ($fhandle, $fname);
}

sub open_file_for_writing($$$){
    my ($dir, $fname, $clobber) = @_;
    if (-e "$dir/$fname" and !$clobber) {
        die "Output file $dir/$fname already exists, stopped ";
    }
    open (my $fhandle, '>', "$dir/$fname")
        or die "Cannot open output file $dir/$fname for writing: $!, stopped ";
    return $fhandle;
}

sub createPlot($$$$){
	 my ($geneName, $noCall, $halfCall, $dir) = @_;
	 my $plotname = $geneName . ".pdf";
	 my $r_cmd = <<"R_CMD";
        R --no-save --silent <<R_SCRIPT
            source("$pathToPlotter")
            setwd("$dir")
            plotCGIgene(file.name="$plotname", gene.name="$geneName", no.call.freq=$noCall, half.call.num=$halfCall, coverage.file="$geneName-refscore.txt", allele.balance.file="$geneName-alleleBalance.txt", cnv.file="$geneName-cnv.txt", variant.file="$geneName-alleleVariants.txt", transcript.file="$geneName-transcripts.txt")
R_SCRIPT
R_CMD

    print "About to execute internal R script for plotting...\n";
    system ($r_cmd) == 0
        or print "Failed to execute R script: $!";
    print "Successfully executed R script\n";
}	

sub usage{
	#print the error message
	my $errorString = shift;
	print "$errorString";
	print "perl plotCGIgene.pl -i -g -d -n\n";
	print "\t-i\t/path/to/genome/data/ASM directory\n";
	print "\t-g\tinput genelist\n";
	print "\t-d\t/path/to/output/directory\n";
	print "\t-r\treference genome version [hg18|hg19]\n";
	print "\t-x\textra bp to add to start/stop of transcript\n";
	print "\t-n\tinteger number basepairs to restrict plotting\n\n";
	print "\nWritten by Jason M Laramie\n\n";
	die;
}

sub getTime{
	#get a timestamp
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	print "$theTime\n";
}

getTime();
#Create a database to look up gene positions
open(IN,$pathToKnownGene) or die "Can't read: $pathToKnownGene: $!\n";
my %ucscGene;
my %alias;
print "Reading $pathToKnownGene...\n";
my $count = 0;
while(<IN>){
	chomp;
	my @line = split("\t", $_);
	#remove alternate assemblies
	if($line[1] eq "chr1" || $line[1] eq "chr2" || $line[1] eq "chr3" || $line[1] eq "chr4" || $line[1] eq "chr5" || $line[1] eq "chr6" || $line[1] eq "chr7" || $line[1] eq "chr8" || $line[1] eq "chr9" || $line[1] eq "chr10" || $line[1] eq "chr11" || $line[1] eq "chr12" || $line[1] eq "chr13" || $line[1] eq "chr14" || $line[1] eq "chr15" || $line[1] eq "chr16" || $line[1] eq "chr17" || $line[1] eq "chr18" || $line[1] eq "chr19" || $line[1] eq "chr20" || $line[1] eq "chr21" || $line[1] eq "chr22" || $line[1] eq "chrY" || $line[1] eq "chrX" || $line[1] eq "chrM"){
		$count++;
		$ucscGene{$line[0]} = $_;
	}
}
print "Found $count transcripts\n";
close(IN);

my %ucscid2name;
open(IN, $pathTokgXref) or die "Can't read: $pathTokgXref: $!\n";
while(<IN>){
	chomp;
	my @line = split("\t", $_);
	$ucscid2name{$line[0]} = $line[4];
}
close(IN);


my %entrez2ucscid;
my %entrez2name;
open(IN, $pathToKnown2locusId) or die "Can't read: $pathToKnown2locusId: $!\n";
while(<IN>){
	chomp;
	my @line = split("\t", $_);
	if(defined($ucscGene{$line[0]})){
		if(!defined($entrez2ucscid{$line[1]})){
			$entrez2ucscid{$line[1]} = $line[0];
			$entrez2name{$line[1]} = $ucscid2name{$line[0]};
		}else{
			$entrez2ucscid{$line[1]} .= "|$line[0]";
		}
	}
}
close(IN);

#gene list
sub processTranscripts($$$){
	my ($transcriptIds, $geneName, $FH) = @_;
	my @transcripts = split("\\|",$transcriptIds);
	my $start = 999999999999; 
	my $end = -1;
	my $chr;
	my $genelength = 0;
	foreach my $id (@transcripts){
		my @line = split("\t", $ucscGene{$id});
		$chr = $line[1];
		if($start > $line[3]){
			$start = $line[3];
		}
		if($end < $line[4]){
			$end = $line[4];
		}
		my @starts=split("\,",$line[8]);
		my @ends=split("\,",$line[9]);
		for(my $i=0; $i<=$#starts; $i++){
			print $FH "$line[0]\t$line[1]\t$line[2]\t$starts[$i]\t$ends[$i]\n";
		}
	}
	print "Using genome coordinates of $chr:$start-$end\n";
	#adjust distance
	$start = $start - $options{x};
	$end = $options{x} + $end;
	
	if($options{x} != 0){
		print "Adjusting genome coordinates to $chr:$start-$end\n";
	}
	$genelength = $end - $start;
	print "Calculating Coverage\n";
	getCoverage($chr,$start,$end, $geneName);
	print "Finding variations and calculating allelic balance\n";
	my ($noCall, $halfCall)=calcAllelicImbalance($chr,$start,$end, $geneName);
	print "Finding ploidy\n";
	calcPloidy($chr,$start,$end,$geneName);
	print "Attempting to plot $geneName...\n";
	if(defined $options{n}){
		if($genelength <= $options{n}){
			createPlot($geneName,$noCall,$halfCall, $options{d});
		}else{
			print "Gene length is $genelength and is greater than the threshold $options{n}\n";
			print "Therefore not plotting $geneName\n";		
		}
	}else{
		createPlot($geneName,$noCall,$halfCall,$options{d});
	}	
}

open(IN, $geneList) or die "Can't open $geneList: $!\n";


my $counter=0;
while(<IN>){
	chomp;
	if($counter == 0){
		$counter++;
		next;
	}
	my @line = split(" ",$_);
	print "Working on entrezID $line[0]\n";
	my $OUT_FH;
	#look up the entrezID
	if(defined($entrez2ucscid{$line[0]})){
		print "Found the transcript(s) $entrez2ucscid{$line[0]} ";
		print "for gene $entrez2name{$line[0]}\n";
		my $txs = $entrez2name{$line[0]} . "-transcripts.txt";
		$OUT_FH = open_file_for_writing($options{d}, $txs, 1);
		processTranscripts($entrez2ucscid{$line[0]}, $entrez2name{$line[0]}, $OUT_FH);
		close($OUT_FH);
	}else{
		print "Transcript IDs not found.  Moving on\n";
	}
}
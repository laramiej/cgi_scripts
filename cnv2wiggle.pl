#!/usr/bin/perl
use strict;

# A Very Simple CNV Variants to wig Converter and BED file generator

#This script uses the cnvDetailsBeta file to create a wiggle track and the cnvSegments
#file to create a BED track.
#update to create the full genome track - must break wiggle into chromosomes
my $usage=  "usage: perl cnv2wiggle.pl cnvDetailsBeta-*-ASM cnvSegmentsBeta-*-ASM chromosome\n";
my $parse=0;
# check that we got 3 args
die $usage if (!defined($ARGV[0]) || !defined($ARGV[1]) || !defined($ARGV[2]));

# arg #1 & #2 ($ARGV[0]) are the file names

# arg #3 is the chromosome (only single tracks are supported
die "arg 3 '$ARGV[2]' must be a chromosome name or the word 'all'\n$usage"
  unless ($ARGV[2]=~/^chr\d+$/ || $ARGV[2]=~/^chr[MXY]$/ || $ARGV[2] eq "all");

# set values from args
my $fileDetails=    $ARGV[0]; #for the wiggle plots
my $fileSegments=    $ARGV[1]; #for the bed track
my $chr=     $ARGV[2];

# open file
if ($fileDetails=~/\.gz$/) { # if it ends with .gz
    open (FILE,"zcat $fileDetails |")  || die "can't open 'zcat $fileDetails'";
} elsif ($fileDetails=~/\.bz2$/) { # if it ends with .bz2
    open (FILE,"bzcat $fileDetails |") || die "can't open 'bzcat $fileDetails'";
} else {
    open (FILE,"<$fileDetails") || die "can't open '$fileDetails'";
}



my $bedFile = "cnv-ALL.bed";
my $bedFileGene = "cnv-overlapgene-ALL.bed";
#open(WIG, ">$wigFile") or die "Can't open $wigFile\n";
open(BED, ">$bedFile") or die "Can't open $bedFile\n";
open(BEDGENE, ">$bedFileGene") or die "Can't open $bedFileGene\n";
# read var file until we get to a line which starts with >
while (<FILE>) {
    last if (/^\>/);
}

# print header line for BED file
#print WIG "track type=wiggle_0\n";
#print WIG "variableStep  chrom=$chr span=2000\n";

my $name = $chr . "CNV";
print BED "track name=$name itemRgb=On\n";

# Here we parse the body of the cnv details file, which looks like:
#
#    0    1      	    2     			   3      		     4      	    5         		6                 7            8            9
# >chr  position   avgNormalizedCvg    gcCorrectedCvg  fractionUnique  relativeCvg     calledPloidy    calledCNVType   ploidyScore  CNVTypeScore
#chr1    1000    83.5    47.2    0.05    1.25    N       hypervariable   0       0
#chr1    3000    103.6   143.0   0.05    1.55    N       hypervariable   0       0
#chr1    5000    112.8   125.6   0.10    1.68    N       hypervariable   0       0
#chr1    7000    122.8   124.9   0.14    1.83    N       hypervariable   0       0
#
print "Parsing $fileDetails to make the wiggle track...\n";
my $WIG;
if($chr eq "all"){
	my $previousChr = "NA";
	while (<FILE>) {
		chomp;
		my @line= split(/\t/);
		if($previousChr eq "NA"){
			#this is the first line in the file
			print "Working on $line[0]\n";
			$WIG = openWigFile($line[0]);
		}elsif($previousChr ne $line[0]){
			#new chr close previous file and open new
			close $WIG;
			print "Working on $line[0]\n";
			$WIG = openWigFile($line[0]);
		}
		if($line[5] ne "N"){ #is this the right chromosome and not a hypervar or invar region
			print $WIG "$line[1]\t$line[5]\n"; 
		}
		$previousChr = $line[0];
	}
}else{
	$WIG = openWigFile($chr);
	while (<FILE>) {
		chomp;
		my @line= split(/\t/);
		if($line[0] eq $chr && $line[5] ne "N"){ #is this the right chromosome and not a hypervar or invar region
			print $WIG "$line[1]\t$line[5]\n"; 
		}
	}
}
close FILE;
close $WIG;

#segments
# 0       1       2             3                     4               5                 6            7                8    
#>chr    begin   end     avgNormalizedCvg        relativeCvg     calledPloidy    calledCNVType   ploidyScore     CNVTypeScore    overlappingGene knownCNV        repeats
#chr1    0       167280  78.4    1.17    N       hypervariable   0       0                       
#chr1    217280  257582  78.4    1.17    N       hypervariable   0       0                       
#chr1    307582  461231  78.4    1.17    N       hypervariable   0       0 

if ($fileSegments=~/\.gz$/) { # if it ends with .gz
    open (FILE,"zcat $fileSegments |")  || die "can't open 'zcat $fileSegments'";
} elsif ($fileSegments=~/\.bz2$/) { # if it ends with .bz2
    open (FILE,"bzcat $fileSegments |") || die "can't open 'bzcat $fileSegments'";
} else {
    open (FILE,"<$fileSegments") || die "can't open '$fileSegments'";
}

while (<FILE>) {
    last if (/^\>/);
}

# read loci
print "Parsing $fileSegments to make the BED track\n";
while (<FILE>) {
	chomp;
	my @line= split(/\t/);
		if($line[5] ne "N" && $line[6] ne "="){
			print BED "$line[0]\t$line[1]\t$line[2]\t$line[5]\t$line[8]\t+\t$line[1]\t$line[2]\t";
   			if($line[6] eq "+"){
    	 		print BED "0,255,0\n";
    		}elsif($line[6] eq "-"){
    			print BED "255,0,0\n";
    		}else{
    			print BED "0,0,0\n";
    		}
		}
		if($line[5] ne "N" && $line[6] ne "=" && $line[9]){ #only want those segements that overlap a gene
			print BEDGENE "$line[0]\t$line[1]\t$line[2]\t$line[5]\t$line[8]\t+\t$line[1]\t$line[2]\t";
   			if($line[6] eq "+"){
    	 		print BEDGENE "0,255,0\n";
    		}elsif($line[6] eq "-"){
    			print BEDGENE "255,0,0\n";
    	}else{
    		print BEDGENE "0,0,0\n";
    	}
	}		
}
close BED;
close BEDGENE;
close FILE;
exit 0;

sub openWigFile{
	my $chr = shift;
	my $wigFile = "cnv-" . $chr . ".wig";
	open(WIG, ">$wigFile") or die "Can't open $wigFile\n";
	print WIG "track type=wiggle_0\n";
	print WIG "variableStep  chrom=$chr span=2000\n";
	my $FH = *WIG;
	return($FH)
}


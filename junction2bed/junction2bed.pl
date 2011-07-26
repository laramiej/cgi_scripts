#!/usr/bin/perl -w
use strict;

my $usage = "perl junction2bed.pl <path/to/evidenceDNB> <chr> <browser start> <browser end> <junctionID 1> <junctionID 2>\n";
die $usage if ($#ARGV <= 3);


# set values from args
my $fileDetails=    $ARGV[0];
my @juncs;
if($ARGV[5]){@juncs=($ARGV[4],$ARGV[5]);}else{@juncs=($ARGV[4])}
 
my $chr=     $ARGV[1];
my $begin=	$ARGV[2];
my $end=	$ARGV[3];



print "browser position $chr:$begin-$end\n";
print "browser hide all\n";


foreach my $j (@juncs){
	if ($fileDetails=~/\.gz$/) { # if it ends with .gz
    	open (FILE,"zcat $fileDetails |")  || die "can't open 'zcat $fileDetails'";
	} elsif ($fileDetails=~/\.bz2$/) { # if it ends with .bz2
    	open (FILE,"bzcat $fileDetails |") || die "can't open 'bzcat $fileDetails'";
	} else {
    	open (FILE,"<$fileDetails") || die "can't open '$fileDetails'";
	}

	while (<FILE>) {
    	last if (/^\>/);
	}
	
	print "track name=\"Junction $j\" description=\"mate pairs in evidence of junction $j\" visibility=squish itemRgb=on\n";
	my $lineCount=0;
	while (<FILE>) {
		chomp;
		$lineCount++;
		my $color;
		my @line = split("\t", $_);
		if($line[0] != $j){next;}
		if($line[5] eq "L" && $line[6] eq "-"){
			$color = "0,255,0";
		}elsif($line[5] eq "R" && $line[6] eq "+"){
			$color = "100,100,0";
		}elsif($line[5] eq "L" && $line[6] eq "+"){
			$color = "255,0,0";
		}elsif($line[5] eq "R" && $line[6] eq "-"){
			$color = "0,0,255";
		}
		print "$line[7]\t$line[8]\t$line[14]\tmatePair$lineCount\t1\t$line[6]\t$line[8]\t$line[14]\t$color\t2\t42,42\t0,";
		my $diff = $line[14] - 42 - $line[8];
		print "$diff\n";
	}
	close(FILE);
}
	
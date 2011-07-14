#!/usr/bin/perl -w

use strict;use Getopt::Std;
my %options=();

getopts("i:d:m", \%options);

# Will create a set of plink files from a testvariants output for only SNPs.
#This assumes unrelated people.  Change the output'd tped file if this is not the case.


#check for input
if (!defined $options{i}){
	usage("option -i not defined\n");
}elsif(!defined $options{d}){
	usage("option -d not defined\n");
}

open(IN, $options{i}) or die "Can't open $options{i}: $!\n";
my $tped = $options{d} . "/plink.tped";
open(OUT, ">$tped") or die "Can't open $options{d}/plink.tped: $!\n";;
my $count = 0;
while(<IN>){
	chomp;
	if($count == 0){
		processHeader($_);
		$count++;
		next;
	}
	if($count % 1000000 == 0){
		print "Working on line $count\n";
	}
	my @line = split("\t", $_);
	if($line[4] ne "snp"){
		next;
	}elsif($line[1] =~ "chrX|chrY|chrM"){
		next;
	}
	my $ref = $line[5];
	my $alleleSeq = $line[6];
	my $missing_flag = 0;
	my $genotypes = ();
	#print "$#line\n";
	for(my $i = 8; $i<= $#line; $i++){
		if($line[$i] =~ "N"){
			#exclude missing by default
			if(!$options{m}){
				$missing_flag=1;
				last;
			}else{
				print "$line[$i]\n";
				$genotypes .= " 0 0";
			}
		}elsif($line[$i] =~ "00"){
			$genotypes .= " $ref $ref";
		}elsif($line[$i] =~ "01|10"){
			$genotypes .= " $ref $alleleSeq";
		}else{
			$genotypes .= " $alleleSeq $alleleSeq";
		}
	}
	if(!$missing_flag){
		$line[1] =~ s/chr//;
		if($line[7] ne ""){
			my @rsNums = split("\;",$line[7]);
			#print "@rsNums\n";
			my ($dbVer,$rsId) = split("\:",$rsNums[0]);
				if($rsId =~ "rs"){
					print OUT "$line[1] $rsId 0 $line[4]";
				}else{
					print OUT "$line[1] $line[0] 0 $line[4]";
				}
		}else{
			print OUT "$line[1] $line[0] 0 $line[4]";
		}
		print OUT "$genotypes\n";	
	}
	$count++;	
}
close(OUT);


sub processHeader{
	my @header=split("\t",$_);
	my @ids = splice(@header,8);
	my $tfam = 	$options{d} . "/plink.tfam";
	open(HEAD, ">$tfam") or die "Can't open plink.ped: $!\n";
	my $countId=0;
	foreach my $id (@ids){
		$countId++;
		print HEAD "$countId\t$id\t0\t0\t?\t-9\n";
	}
	close(HEAD);	
}

sub usage {
	#print the error message
	my $errorString = shift;
	print "$errorString";
	print "perl createPlink.pl -i -d\n";
	print "\t-i\tA test variant file\n";
	print "\t-d\toutput directory\n";
	print "\t-m\tflag to include missing data (these will be encoded as NN even is it is a half-call\n";
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

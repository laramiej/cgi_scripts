#!/usr/bin/perl -w

use strict;use Getopt::Std;
my %options=();
use File::Basename;

# Options to run callDiff (These are set at the default).  
my $scoreThresholdA = 10;
my $scoreThresholdB = 3;
my $distance = 500;
my $minlength = 1000;


#using a sample file that contains tumor / normal pairs will run junction diff between the tumor / normal pair matches.  The output files are then compared against each other to see what junctions are in common between the tumors.
getopts("s:d:b:htk", \%options);

#check for input
if (!defined $options{s}){
	usage("option -s not defined\n");
}elsif(!defined $options{d}){
	usage("option -d not defined\n");
}

#Must include full path to cgatools here if cgatools is part of your $PATH then leave this alone.
my $path2cgatools = "cgatools";

#must include full path to reference datasets
my ($path2repeatMask,$path2genes,$path2reference);
if($options{b} eq "hg18"){
	$path2repeatMask = "/rnd/home/jbaccash/TestData/cgatools/0.0.0/ucsc/rmsk36.tsv.gz";
	$path2genes = "/rnd/home/jbaccash/TestData/cgatools/0.0.0/ucsc/gene36.tsv.gz";
	$path2reference = "/sup/sv-01/data/ref/hg18.crr";
}elsif($options{b} eq "hg19"){
	$path2repeatMask = "/rnd/home/jbaccash/TestData/cgatools/0.0.0/ucsc/rmsk37.tsv.gz";
	$path2genes = "/rnd/home/jbaccash/TestData/cgatools/0.0.0/ucsc/gene37.tsv.gz";
	$path2reference = "/sup/sv-01/data/ref/build37.crr";
}else{
	usage("genome build -b must be hg18 or hg19\n");
}

my $allJunctionsPrefix = "allJunctionsBeta-";
my $highConfidencePrefix = "highConfidenceJunctionsBeta-";
my $varPrefix = "var-";


sub find_file($$){
    my ($dir, $fmask) = @_;
    opendir DIR, "$dir" or die "Cannot open directory $dir: $!, stopped ";
    my @files = grep { /$fmask/ } readdir DIR;
    (scalar @files == 1)
        or die (((scalar @files == 0) ? 'No' : 'Multiple') .
            " files matching file mask $fmask exist in directory $dir\n" .
            join("\n", @files) .
            ", stopped ");
    my $fname = $dir . "/" . $files[0];
    return ($fname);
}


#Step 1 run through the samples file and run junctionDiff
open(IN, $options{s}) or die "Can't open $options{s}: $!\n";
my @junctionDiffFiles;
my $count = 0;
if($options{k}){
	print "Skipping Step 1: Not Comparing Tumor to matched normal\n";
}else{
	print "Executing Step 1: Comparing Tumor to match normal\n";
	my $header_ref;
	while(<IN>){
		chomp;
		if($count == 0){
			#skip the header. Could check this if need be.
			my @header = split(" ",$_);
			$header_ref = findCol(\@header);
			$count++;
			next;
		}
		my @line = split(" ",$_);
	
		my $normaljunctions = find_file("$line[$$header_ref{'normal'}]/ASM/SV", qr{^$allJunctionsPrefix.*});
		my $tumorjunctions = find_file("$line[$$header_ref{'tumor'}]/ASM/SV", qr{^$highConfidencePrefix.*});
	
		my $outputPrefix = $options{d} . "/" . $line[$$header_ref{'DNAid'}] . "-notTumor-";
		my $outPrefix = $line[$$header_ref{'DNAid'}] . "-notTumor-diff-";
		push(@junctionDiffFiles, $outPrefix);
		#SAMPLE A IS TUMOR 
		#sample B is matched normal
	
		my $command = $path2cgatools . " junctiondiff --beta --junctionsA " . $tumorjunctions . " --junctionsB " . 	$normaljunctions . " --scoreThresholdA " . $scoreThresholdA . " --scoreThresholdB " .  $scoreThresholdB . " --distance " . $distance . " --minlength " . $minlength . " --statout --output-prefix " . $outputPrefix . " --reference " .$path2reference;
		if($options{t}){
			print "In Debug Mode:\n$command\n";
			die;
		}else{
			print "Going to execute full command\n";
			print "$command\n";
			system ($command) == 0
        	or die "Failed to execute command: $!, stopped ";
    		print "Successfully executed command\n";
		}
	}
}

#step 2 is to buffer the junctionIDs so that we can see what is removed.
#buffer first file
my %junctionIds; #<- this will be the superset
$count = 0;
if($options{k}){
	#since we skipped step 1, the sample file will provide the files to compare
	my $header_ref;
	while(<IN>){
		chomp;
		if($count == 0){
			#skip the header. Could check this if need be.
			my @header = split(" ",$_);
			$header_ref = findCol(\@header);
			$count++;
			next;
		}
		my @line = split(" ",$_);
		my $normalJunctions = find_file("$line[$$header_ref{'normal'}]/ASM/SV", qr{^$allJunctionsPrefix.*});
		my $normalHighJunctions = find_file("$line[$$header_ref{'normal'}]/ASM/SV", qr{^$highConfidencePrefix.*});
		if($options{h}){
			#only use the high confidence junctions
			push(@junctionDiffFiles, $normalHighJunctions);
		}else{
			#use alljunctions
			push(@junctionDiffFiles, $normalJunctions);
		}
	}
	close(IN);
}


	my $firstFilemask = shift(@junctionDiffFiles);
	my $firstfullFile;
	if($options{k}){
		$firstfullFile = $firstFilemask;
	}else{
		$firstfullFile = find_file($options{d},qr{^$firstFilemask});
	}
		
	print "\n\nStarting Step 2: collecting the same junctions\nUsing $firstfullFile as the superset file\n\n";

	open(IN,$firstfullFile) or die "Can't open $firstfullFile: $!\n";
	open(OUT, ">$options{d}/mergedSV.txt") or die "Can't open outfile: $!\n";
	$count=0;
	while(<IN>){
		chomp;
		$count++;
		if($count<12){
			print OUT "$_\n";
			next;
		}
		
		my @line = split("\t",$_);
		#hash the junction ID
		$junctionIds{$line[0]} = $_;
	}
	print "Found " . ($count-12) . " junctions in $firstfullFile\n\n";
	close(IN);

	my @comparisonFiles;
	for(my $i = 0; $i<=$#junctionDiffFiles; $i++){
		#this is the difference file for the next
		my $fullFile;
		if($options{k}){
			$fullFile = $junctionDiffFiles[$i];
		}else{
			$fullFile = find_file($options{d},qr{^$junctionDiffFiles[$i]});
		} 
		
		open(IN,$fullFile) or die "Can't open $fullFile: $!\n";
		my($filename, $directories,$outputPrefix,$outputPre);
		if($options{k}){
			($filename, $directories) = fileparse($fullFile);
			$outputPrefix = "comparison-" . $filename;
			$outputPre = "comparison-" . $filename . "diff-";
		}else{
			$outputPrefix = "comparison-" . $junctionDiffFiles[$i];
			$outputPre = "comparison-" . $junctionDiffFiles[$i] . "diff-";
		}
		
		push(@comparisonFiles,$outputPre);
		my $command = $path2cgatools . " junctiondiff --beta --junctionsA " . $firstfullFile . " --junctionsB " . $fullFile . " --scoreThresholdA " . $scoreThresholdA . " --scoreThresholdB " . $scoreThresholdB . " --distance " . $distance . " --minlength " . $minlength . " --statout --output-prefix " . $outputPrefix . " --reference " .$path2reference;
		print "Comparing $firstfullFile to $fullFile\n";
		print "Going to execute full command\n";
		print "$command\n";
		system ($command) == 0
        	or die "Failed to execute command: $!, stopped ";
    	print "Successfully executed command\n";
	}
	print "Finished comparing the samples\n\n\n";

	foreach my $file (@comparisonFiles){
		my $fullFile = find_file($options{d},qr{^$file});
		print "Opening $fullFile\n";
		$count = 0;
		open(IN,$fullFile) or die "Can't open $fullFile: $!\n";
		while(<IN>){
			$count++;
			if($count<12){
				next;
			}
			my @line = split("\t",$_);
			if(defined($junctionIds{$line[0]})){
			#if we see a junction in a junctiondiff file then it is different.  Remove it.
				delete($junctionIds{$line[0]});
			}
		}
		print "Found " . ($count-12) . " junctions that differ between files\n";
		close(IN);
	}

	my $junctCount = 0;
	foreach my $id (sort {$a <=> $b} keys %junctionIds){
		$junctCount++;
		print OUT "$junctionIds{$id}\n";
	}
	close(OUT);
	print "\n\nFound $junctCount overlapping junctions.  These are outputted to $options{d}/mergedSV.txt\n";	


sub findCol {
	my @header = split(" ",$_);
	my $count = 0;
	my %cols;
	foreach my $elm (@header){
		$cols{$elm} = $count;
		$count++;
	}
	return(\%cols);
} 

sub usage {
	#print the error message
	my $errorString = shift;
	print "$errorString";
	print "perl findOverlappingJunctions.pl -s -d -h -f -t\n";
	print "\t-s\tA tab-delim sample manifest\n";
	print "\t-d\tfull path of the output directory\n";
	print "\t-b\tgenome build [hg18 or hg19]\n";
	print "\t-h\tflag to use high confidence calls\n";
	print "\t-k\tflag to skip step one which is to compare tumor to matched normal\n";
	print "\t-t\t(optional) for debugging\n";
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

#!/usr/bin/perl
use strict;
$| = 1;

# Rick_ListVariants
# Takes as input list of var files
# Extracts a non-redundant set of variations
# Adds a count of how many times each variation occurs in the var files [optional]
# Unlike cgatools listvariants, does not attempt to canonicalise the variants
# So the same variation in different forms across genomes will remain entered as different variations

# perl prog file(s) [flag] filename
# ie
# perl Rick_ListVariants \
# --Input_File varfile_1 \
# --Input_File varfile_2 \
# ...
# --Input_File varfile_n \
# --Export_Counts \ [optional, no value to be entered]
# --Output_File filename
# eg
# perl /Users/rtearle/Documents/Programming/Perl/Scripts/Dev/Rick_ListVariants_0_0_1.pl \
# --Input_File /Yoruban_Trio_1100_37/GS19238-1100-37/GS00028-DNA_A01/ASM/var-GS19238-1100-37-ASM.tsv.bz2 \
# --Input_File /Yoruban_Trio_1100_37/GS19239-1100-37/GS00028-DNA_B01/ASM/var-GS19239-1100-37-ASM.tsv.bz2 \
# --Input_File /Yoruban_Trio_1100_37/GS19240-1100-37/GS00028-DNA_C01/ASM/var-GS19240-1100-37-ASM.tsv.bz2 \
# --Export_Counts \
# --Output_File YRI_Trio_ListVariants.tsv

# var file fields
# >locus ploidy haplotype chromosome begin end varType reference alleleSeq totalScore hapLink xRef

# Parsing and storing input parameters
# Only input_file fields can be repeated
# input paramaters are case insensitive

print "\nProcessing input parameters\n";
my $NrParams;
my %ExpectedParams =  GetExpectedParams (); # list of expected parms
my %EnteredParams = GetEnteredParams (); # list of entered params
TestParameters (%EnteredParams); # compares expectd and entered params

# Input Files
my $FilesIn = $EnteredParams{input_file}; # ptr to list of input files
print "input files:\n";
my $NrInputFiles = int(@$FilesIn);
foreach my $File (@$FilesIn) {unless (-f $File) {die "Testvariants input file $File not found\n";}} # requires existing files
for my $n (0.. $NrInputFiles-2) # look for duplicates in list
{
	for my $m ($n+1.. $NrInputFiles-1) {if ($$FilesIn[$n] eq $$FilesIn[$m] ) {die "File $$FilesIn[$n] is repeated in input file list\n";}}
}

# Ouput File
my $FileOut = $EnteredParams{output_file}; # output file
if (-f $FileOut) # ouput file exists, create a new one based on the name
{
	print "Output file $FileOut exists, modifying to unique file name ";
	$FileOut =~ /^(.+?)\./; # find name without extensions
	my $Stub = $1; # set stub to name without extensions
	$FileOut =~ /(\..+)?$/; # get extension(s)
	my $Ext = $1; # set ext to extensions
	my $n = 1; # n will increment to find a unique name
	my $Suff = ""; # suff tracks n
	while (-f $Stub.$Suff.$Ext) {$Suff = "-$n"; $n++;} # loop till we have a new unique filename
	$FileOut = $Stub.$Suff.$Ext; # file out now has same name, same extensions, but also -n at the end of the name, making it unique
	print "$FileOut\n";
}

my $ExportCounts = $EnteredParams{export_counts}; # flag to export counts

#print "Files\n",join("\n",@$FilesIn),"\n\n";
#print "Ouput Dir\n$DirectoryOut\n";
#print "Ouput File\n$FileOut\n";

# Loading chr nrs, setting up var hash
my @ChrNames = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
				'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
				'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'); # using this array forces the output order of chrs into the correct order
my %Vars; # hash to store var records in an array for each chr
foreach my $Chr (@ChrNames) {$Vars{$Chr} = {};} # print "$Chr\t";}  # set up hash of hashes, one for each chr
#print "\n"; #exit;

# Extract Column Header
my $IN = OpenFile ($$FilesIn[0]); # open the first file with correct file format
my $Header = GetHeader ($IN); # get header
unless ($Header) {close $IN; $IN = OpenFile ($$FilesIn[0]);} # if there is no header, close and reopen file, ie start file again
my $ColHeader = <$IN>; # get col header, first remaining line
#$ColHeader =~ s/^>\w+\t\w+\t\w+\t/>/; # remove first three fields: locus ploidy and allele
chomp $ColHeader;
my @x = split("\t",$ColHeader);
$ColHeader = join "\t", @x[3..8,11];
$ColHeader .= "\tCount" if $ExportCounts; # add count to end of col header, ie count of how many times var is found in the gene files
# print "$ColHeader\n"; exit;

# Process Files
my $RecCount = 0; # total nr recs
foreach my $File (@$FilesIn)
{
	my $Count = 0; # count recs for this file
	print "Processing file $File\n";
	my $IN = OpenFile ($File); # open the file with correct file format
	my $Header = GetHeader ($IN); # get  header
	unless ($Header) {close $IN; $IN = OpenFile ($$FilesIn[0]);} # if no header, close and reopen file, ie start file again
	<$IN>; # get col header ie first remaining line, not using so not keeping
	while (<$IN>) # loop through remainder of file ie data
	{
		my ($Chr, $Rec) = ExtractFields ($_); # sub extracts wanted fields in rec as array
		#print $Rec,"\n";
		#print "Chr: $Chr\tRec: ",$Rec,"*\n" if $Chr;
		#if ($Chr) {print "HERE: ",$Chr,"\t",$Rec,"\n"; exit;}
		if ($Chr) # sub above returns null if not a var
		{
			$Vars{$Chr}->{$Rec}++; # hash for each chr, containing a hash with rec as key, incrementing count for this key
			$Count++; # increment count of vars
		}
	}
	print "Nr matched records for this file: $Count\n";
	$RecCount += $Count; # add this file's count to total count
	close $IN;
}
print "Nr matched records across all files:\t $RecCount\n"; # total count

print "Saving to file $FileOut ...\n";
open my $OUT, ">", $FileOut;
print $OUT "$ColHeader\n";

$RecCount = 0; # reuse total count for nr of non-reduntant vars

foreach my $Chr (@ChrNames) # sort records in each chr array and print with count
{
	if ($ExportCounts)
	{
		foreach my $Rec (sort {SortStringsasArrays ($a, $b)} keys %{$Vars{$Chr}}) # using sub to sort on begin, end fields
		{
			print $OUT "$Rec\t$Vars{$Chr}->{$Rec}\n"; # printing rec and count
			$RecCount++; # increment count of vars
		}
	}
	else
	{
		foreach my $Rec (sort {SortStringsasArrays ($a, $b)} keys %{$Vars{$Chr}}) # using sub to sort on begin, end fields
		{
			print $OUT "$Rec\n"; # printing rec but not count
			$RecCount++; # increment count of vars
		}
	}
}
print "Nr saved (non-reduntant) records:\t $RecCount\n"; # count of non-redundant vars, c/f all vars above


###########################################################################
#                                   SUBS                                  #
###########################################################################

sub GetExpectedParams
{
	my %Hash = # hash to store expected params
	(
		"input_file" => [],
		"output_file" => "",
		"export_counts" => 0,
	);
	$NrParams = int keys %Hash;
	return %Hash;
}

sub GetEnteredParams
{
	# Processing @ARGV
	my %Hash;
	my @ARGVs = split /--/, join (" ",@ARGV); # split args on --, into array
	for my $n (1..$#ARGVs) # parse each [nb arg 0 is empty so ignored]
	{
		$ARGVs[$n] =~ s/\s+$//; # remove any trailing spaces
		my ($Key, $Val) = split / /, $ARGVs[$n], 2; # put first element into key, any other elements into val
		$Key = lc $Key; # making all lower case so will be case insensitive
		if ($Key eq "input_file") # multiple entries expected, setting up array
		{
			push @{$Hash{$Key}}, $Val; # add input to input hash
		}
		else
		{
			$Hash{$Key} = $Val; # make a hash entry out of key and val
		}
	}
	return %Hash; # hash now has each --entry param, with associated values
}

sub TestParameters
{
	# Test if parameters are filled, all are present, no extras present - not fully tested so may be redundant/incomplete
	my $ArgErrors = 0;
	foreach my $ExpArg (keys %ExpectedParams)
	{
		if ($ExpArg eq "export_counts") {$EnteredParams{$ExpArg} = 1;}
		elsif (!defined $EnteredParams{$ExpArg}) #  if there is no entry for this expected param
		{print "Parameter --$ExpArg appears to be missing from arguments\n"; $ArgErrors++;}
		elsif ($EnteredParams{$ExpArg} == -1) # if no entry has been entered from the input list
		{print "No value for parameter --$ExpArg in arguments"; $ArgErrors++;}
	}
	if (int keys %EnteredParams > $NrParams) {die "There are too many parameter in the arguments";} # die if too many params
	if ($ArgErrors) {die "$ArgErrors detected in input parameters";} # die if a param not entered
}

sub OpenFile
{
    my $File = shift;
    my $FH;

    if ($File =~ /.bz2$/)
    {
		open ($FH, "bzcat $File |") or die ("$!: can't open file $File");
    }
	elsif ($File =~ /.gz$/)
    {
		open ($FH, "gunzip -c $File |") or die ("$!: can't open file $File");
    }
    elsif ($File =~ /.tsv$/ or $File =~ /.txt$/)
    {
		open ($FH, "cat $File |") or die ("$!: can't open file $File");
    }
    else
    {
		print ("Do not recognise file type for file $File.\nOpening as text file\n");
		open ($FH, "cat $File |") or die ("$!: can't open file $File");
    }
    return $FH;
}

sub GetHeader
{
    my $FH = shift;
    my @Header;
    my $Count = 0;
    while (<$FH>) # loop until a line is empty
    {
		if ($_ eq "\n") # exit when empty line
		{
			chomp @Header; # removing all terminated returns
			return \@Header ; # return ref to array
		}
		else
		{
			push @Header, $_;
		}
		return my @Array if $Count++ > 50; # too many lines for a header, must be no header, return empty array
    }
}

sub ExtractFields # expects a gene file rec, strips out file specific fields, gets chr
{
	# >locus ploidy haplotype chromosome begin end varType reference alleleSeq totalScore hapLink xRef
	my $Rec = shift;
	my $Wanted;
	my @Fields;

	if ($Rec =~ /\tsnp|ins|del|sub\t/) # only take recs which contain variants
	{
		chomp $Rec; # remove return
		@Fields = split "\t", $Rec;
		$Wanted = join ("\t",@Fields[3..8,11]);
		#print "$Rec\n";
		#print $Wanted,"\n";
	}
	return ($Fields[3], $Wanted);
}

sub SortStringsasArrays # sorts based on begin and end of two recs
{
	my $String1 = shift; # first string
	my $String2 = shift; # second string

	my @Array1 = split "\t", $String1; # put fields into array
	my @Array2 = split "\t", $String2;

	# array[1] is begin, array[2] is end, returning order based on these fields
	if ($Array1[1] < $Array2[1]) # begin of 1 < begin of 2
	{
		return -1;
	}
	elsif ($Array1[1] == $Array2[1]) # begin of 1 == begin of 2
	{
		if ($Array1[2] < $Array2[2]) # end of 1 < end of 2
		{
			return -1;
		}
		elsif ($Array1[2] == $Array2[2]) # end of 1 == end of 2
		{
			return 0;
		}
		else # end of 1 > end of 2
		{
			return 1;
		}
	}
	else # begin of 1 > begin of 2
	{
		return 1;
	}

}
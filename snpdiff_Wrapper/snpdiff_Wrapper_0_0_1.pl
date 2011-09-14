#!/usr/bin/perl
use strict;
use File::Basename;
$| = 1;

# snpdiff_Wrapper_0_0_1
# Script to pass snpdiff a modified copy of a genotype file, avoiding modifying the original file
# snpdiff requires genotype data with specific column headers and 0 based positions
# But users may not want to modify genotype files to fit this format
# This wrapper takes a genotype file, creates a temporary copy, modifies and passes it to snpdiff
# When snpdiff finishes, the temporary copy is deleted

# Required field names are: Chromosome, Offset0Based, GenotypesStrand Genotypes
# User can provide current names or field nrs (0 based) for these fields using parameters below
# User can provide offset of positions, either 0 or 1 based

# perl prog file dir filename
# ie
# perl snpdiff_Wrapper \
# --Reference reference_file \
# --Variants var_file_ \
# --Genotypes genotype_file_ \
# --Output_Prefix output_filestub \
# --Offset integer [0|1]
# --Reports [output,stats,verbose]
# --Delete_TempFile [optional, yes|no, default is no]
# --Chromosome_Field name [or --Chromosome_Field_Nr integer] \
# --Offset0Based_Field name [or --Offset0Based_Field_Nr integer] \
# --GenotypesStrand_Field name [or --GenotypesStrand_Field_Nr integer] \
# --Genotypes_Field name [or --Genotypes_Field_Nr integer \

# eg
# perl /Users/rtearle/Documents/Programming/Perl/Scripts/Dev/snpdiff_Wrapper_0_0_1.pl \
# --Reference /Users/rtearle/Documents/TBF/Reference_Files/Reference_Genomes/build37.crr \
# --Variants /Users/rtearle/Documents/TBF/snpdiff_Wrapper-Test/var-GS19240-1120-37-21-ASM.tsv.bz2 \
# --Genotypes /Users/rtearle/Documents/TBF/snpdiff_Wrapper-Test/NA19240_HapMap_Infinium_Genotypes_37_chr21-Modified.tsv \
# --Output_Prefix /Users/rtearle/Documents/TBF/Test-snpdiff_Wrapper \
# --Offset 1 \
# --Delete_TempFile no
# --Reports Output,Stats \
# --Chromosome_Field Chr \
# --Offset0Based_Field Position \
# --GenotypesStrand_Field Orient \
# --Genotypes_Field Call \

# or

# perl /Users/rtearle/Documents/Programming/Perl/Scripts/Dev/snpdiff_Wrapper_0_0_1.pl \
# --Reference /Users/rtearle/Documents/TBF/Reference_Files/Reference_Genomes/build37.crr \
# --Variants /Users/rtearle/Documents/TBF/snpdiff_Wrapper-Test/var-GS19240-1120-37-21-ASM.tsv.bz2 \
# --Genotypes /Users/rtearle/Documents/TBF/snpdiff_Wrapper-Test/NA19240_HapMap_Infinium_Genotypes_37_chr21-Modified.tsv \
# --Output_Prefix /Users/rtearle/Documents/TBF/Test-snpdiff_Wrapper \
# --Offset 1 \
# --Reports Output,Stats \
# --Chromosome_Field_Nr 2 \
# --Offset0Based_Field_Nr 3 \
# --GenotypesStrand_Field_Nr 4 \
# --Genotypes_Field_Nr 5


print "\nProcessing input parameters\n";
my $NrParams;
my %ExpectedParams =  GetExpectedParams (); # list of expected parms
my %EnteredParams = GetEnteredParams (); # list of entered params
TestParameters (%EnteredParams); # compares expected and entered params

# Input Files
my $ReferenceFile = $EnteredParams{reference}; # reference file
unless (-f $ReferenceFile) {die "Reference file $ReferenceFile not found\n";} # requires existing file
my $VarFile = $EnteredParams{variants}; # variant file
unless (-f $VarFile) {die "Variant file $VarFile not found\n";} # requires existing file
my $GenotypeFile = $EnteredParams{genotypes}; # genotype file
unless (-f $GenotypeFile) {die "Genotype file $GenotypeFile not found\n";} # requires existing file

# Check if Output Dir exists
my $FileOutPrefix = $EnteredParams{output_prefix}; # output dir
$FileOutPrefix =~ /(.+\/)/g; # matching string ending in '/',
my $DirOut = $1; # not really necessary, keeping in case we use dir out in a later version
if ($DirOut) # if string exists, test dir exists
{
	unless (-d $DirOut) {die "Output directory $1 not found\n";} # requires existing dir
}

# Temp Output File
my $Seconds = time;
my $TempGenotypeFile = $FileOutPrefix."-Genotypes$Seconds.tmp"; # temp genotype file
open my $OUT, ">", $TempGenotypeFile; # open temp genotype file for output
#print "$ReferenceFile\n$VarFile\n$GenotypeFile\n";

# Input Params
my $Reports = $EnteredParams{reports} || $ExpectedParams{reports}; # get reports list, otherwise use all reports
my $DeleteTempFile = lc $EnteredParams{delete_tempfile} || "yes"; # whether to delete temp file, default is yes
my $OffSet = $EnteredParams{offset}; # get offset, either 0 or 1
#print "$Reports\n$OffSet\n";

my %Defaults = qw (chromosome Chromosome offset0based Offset0Based genotypesstrand GenotypesStrand genotypes Genotypes); # hash of default field names
my %Fields = qw (chromosome "" offset0based "" genotypesstrand "" genotypes ""); # hash to store current field names or nrs
foreach my $Field (keys %Fields) # loop through the fields hash
{
	if (defined $EnteredParams{$Field."_field"}) # field name exists
	{
		$Fields{$Field} = $Field."_field"; # set value to field name
	}
	elsif (defined $EnteredParams{$Field."_field_nr"})  # field nr exists
	{
		$Fields{$Field} = $Field."_field_nr"; # set value to field nr
	}
	else
	{
		delete $Fields{$Field}; # no value, will use current field name, should be the default, will test when file opened
	}
}

print "\nPreparing Modified Genotype File\n";
# Assuming genotype file has a column header by not a header - can modify if necessary
my $IN = OpenFile ($GenotypeFile);
my $ColHeader = <$IN>;
chomp $ColHeader;
#print "Before:\n$ColHeader\n";

print "\nModifying Column Header\n" if int keys %Fields; # no need to do this if there are no name changes to make
foreach my $Field (keys %Fields) # loop through looking for field names to change
{
	if ($Fields{$Field} !~ /_nr/) # field name not nr, so process
	{
		$ColHeader =~ s/$EnteredParams{$Fields{$Field}}/$Defaults{$Field}/;
	}
}
#print "First Pass:\n$ColHeader\n";
my @ColHeader = split "\t", $ColHeader; # converting header in to array to look for specific names by field nr
foreach my $Field (keys %Fields) # loop through looking for field nrs to change
{
	if ($Fields{$Field} =~ /_nr/) # field nr not name, so process
	{
		$ColHeader[$EnteredParams{$Fields{$Field}}] = $Defaults{$Field};
	}
}
$ColHeader = join "\t", @ColHeader; # convert header array back to string
#print "Second Pass:\n$ColHeader\n";

# Check field names
my $ErrorCount = 0;
foreach my $Field (keys %Defaults) # loop through default field names
{
	unless ($ColHeader =~ /(?:^|\t)$Defaults{$Field}(?:\t|$)/) # if col header does not contain default field name
	{
		print "Fatal error: Column header must contain a field named $Defaults{$Field} and does not\n"; # tell the user
		$ErrorCount++; # increment error count
	}
}
if ($ErrorCount)
{
	exit; # exit if they are not all present
}
else
{
	print "Column Header contains all required field names: ";
	foreach my $Key (keys %Defaults) {print "$Defaults{$Key} ";}
	print "\n";
}

#$ColHeader =~ s/Chromosome/C/g; print "$ColHeader\n"; # error testing
print $OUT "$ColHeader\n"; # col header to file

# Processing data in genotype file
if ($OffSet) # offset is 1, need to make 0
{
	my $ColNr = -1; # column nr containing positions
	for my $n (0..int(@ColHeader)-1) # loop through col fields to find nr of offset col
	{
		if ($ColHeader[$n] eq "Offset0Based") # looking for field called Offset0Based
		{
			$ColNr = $n; # found, storing nr of field
			last;
		}
	}
	if ($ColNr == -1) {die "Cannot find a field named Offset0Based\n";} # in case field is not found

	while (<$IN>) # looping through records of genotype file
	{
		chomp; # remove return
		my @Fields = split "\t", $_; # split fields to array
		$Fields[$ColNr]--; # decrement pos field
		print $OUT join "\t", @Fields, "\n"; # output record to file
	}
}
else # no need to change offset, just output records
{
	while (<$IN>)
	{
		print $OUT $_; # output record to file
	}
}

print "\nRunning cgatools calldiff\n";
open my $Process, "cgatools snpdiff --reference $ReferenceFile --variants $VarFile --genotypes $TempGenotypeFile --output-prefix $FileOutPrefix --reports $Reports |";
while (<$Process>) {print $_;}

if ($DeleteTempFile eq "yes")
{
	print "\nDeleting temporary file $TempGenotypeFile\n";
	open my $Process, "rm -f $TempGenotypeFile |";
	while (<$Process>) {print $_;}
}

print "\n";

###########################################################################
#                                   SUBS                                  #
###########################################################################

sub GetExpectedParams
{
	my %Hash = # hash to store expected params
	(
		reference => "",
		variants => "",
		genotypes => "",
		output_prefix => "",
		reports => "output,stats,verbose",
		offset => 0,
		delete_tempfile => "Yes",
		chromosome_field => "Chromosome",
		offset0based_field => "Offset0Based",
		genotypesstrand_field => "GenotypesStrand",
		genotypes_field => "Genotypes",
		chromosome_field_nr => -1,
		offset0based_field_nr => -1,
		genotypesstrand_field_nr => -1,
		genotypes_field_nr => -1,
	);
	$NrParams = 11; # first 7 plus 4 of the next 8 in pairs
	#foreach my $Key (keys %Hash) {print "$Key   $Hash{$Key}\n";
	return %Hash;
}

sub GetEnteredParams
{
	# Processing @ARGV
	my %Hash;
	my @ARGVs = split /--/, join (" ",@ARGV); # split args on --, into array
	for my $n (1..$#ARGVs) # parse each [nb arg 0 is empty so ignored]
	{
		$ARGVs[$n] =~ s/^\s+//; # remove any leading spaces
		$ARGVs[$n] =~ s/\s+$//; # remove any trailing spaces
		my ($Key, $Val) = split / /, $ARGVs[$n], 2; # put first element into key, any other elements into val
		#$Key = lc $Key; # make lower case, ie case insensitive
		$Hash{lc $Key} = $Val; # make a hash entry out of key and val
	}
	#print "Ent: ", join " ", keys %Hash, "\n\n";
	return %Hash; # hash now has each --entry param, with associated values
}

sub TestParameters
{
	# Test if parameters are filled, all are present, no extras present - not fully tested so may be redundant/incomplete
	my $ArgErrors = 0;
	#foreach my $ExpArg (keys %ExpectedParams) {print "$ExpArg\t$ExpectedParams{$ExpArg}\n";}
	#foreach my $ExpArg (keys %EnteredParams) {print "$ExpArg\t$EnteredParams{$ExpArg}\n";}
	foreach my $ExpArg (keys %ExpectedParams)
	{
		#print "$ExpArg $EnteredParams{$ExpArg}\n";
		if ($ExpArg =~ /reports|chromosome_|offset0based_|genotypesstrand_|genotypes_/) # "optional" params, not tested
			{
				next;
			}
		elsif ($ExpArg eq "offset")
		{
			if (defined $EnteredParams{$ExpArg} and $EnteredParams{$ExpArg} != 0 and $EnteredParams{$ExpArg} != 1) # entry in arg list, must be 0 or 1
			{
				{print "--$ExpArg must be 0 or 1, not $EnteredParams{$ExpArg}\n"; $ArgErrors++;}
			}
		}
		elsif ($ExpArg eq "delete_tempfile")
		{
			$EnteredParams{$ExpArg} = lc $EnteredParams{$ExpArg};
			if (defined $EnteredParams{$ExpArg} and $EnteredParams{$ExpArg} ne "yes" and $EnteredParams{$ExpArg} ne "no") # entry in arg list, must be 0 or 1
			{
				{print "--$ExpArg must be Yes or No, not $EnteredParams{$ExpArg}\n"; $ArgErrors++;}
			}
		}
		else # $ExpArg = reference, variants, genotypes, output_prefix
		{
			if (! defined $EnteredParams{$ExpArg}) # no entry has been entered from the input list
			{print "No value for --$ExpArg in arguments\n"; $ArgErrors++;}
		}
	}

	if ($ArgErrors)
		{die "$ArgErrors error(s) detected in input parameters, please review\n";} # die if a param not entered
	if (int keys %EnteredParams > $NrParams)
		{print "There are too many parameters in the arguments, continuing anyway\n";} # die if too many params
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

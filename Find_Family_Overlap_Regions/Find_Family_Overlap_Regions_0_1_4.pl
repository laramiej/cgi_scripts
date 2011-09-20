#!/usr/bin/perl
use strict;
use File::Basename;
$| = 1;

# Find_Family_Overlap_Regions
#
# Determines which regions of a set of child genomes have been inherited in common.
# Takes in to account mode of inheritance and whether the children share (afffected) or don't share (unaffected).
#
# input from child comparison state files from Family_Inheritqnce_State_Analysis.pl
# eg *_Ch1_Ch2_Pat_States_GT_5.tsv
# Each child must be entered as affected or unaffected
# Regions shared by affected siblings but not shared by unaffected siblings will be outputted
# user can select whether to output a header.
# Header structure is the same as for CompleteGenomics headers.
#
# The mode of inheritance determines which files are used.
# Files which do not match mode of inheritance are ignored.
# So eg a set of files can be inputted and the mode of inheritance varied.
#
# Format:
# perl Family_Overlap_Region_Finder \
# --Input_Mat input_mat_state_file \ [multiple allowed]
# --Input_Pat input_mat_state_file \ [multiple allowed]
# --output_dir output_dir \
# --Header [optional]
# --Affected nr [nr] \ [optional but one of Affected and Unaffected must be present]
# --Unaffected nr [nr] \ [optional as above]
# --Inheritance_Mode [maternal_dominant|paternal_dominant|recessive]
# eg
# perl /hoard/home/rtearle/perl/Find_Family_Overlap_Regions_0_1_1.pl \
# --Input_Mat /Comparisons/Quartet_testvariants-chr20-22_4_Ch1_Ch2_Mat_States_GT_5.tsv \
# --Input_Mat /Comparisons/Quartet_testvariants-chr20-22_4_Ch1_Ch3_Mat_States_GT_5.tsv \
# --Input_Pat /Comparisons/Quartet_testvariants-chr20-22_4_Ch1_Ch2_Pat_States_GT_5.tsv \
# --Input_Pat /Comparisons/Quartet_testvariants-chr20-22_4_Ch1_Ch3_Pat_States_GT_5.tsv \
# --Output_Dir /data/Family_Analysis/Comparisons/ \
# --Header
# --Inheritance_Mode recessive
# --Affected 1 2
# --Unaffected 3

# Rick Tearle 2011

# Input (Chx_Chy_Par_States) fields
# Chr Begin End State Length Records

my %Colours = SetColours (); # colours to use
my @ChrNames = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
				'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22');

# Parsing and storing input parameters
# Input Files, Affecteds and Unaffecteds may be repeated
print "\nProcessing input parameters\n";
my %ExpectedParams =  GetExpectedParams ();
my %EnteredParams = GetEnteredParams ();
TestParameters (); # no return val used, dies in subroutine if no correct

# Set up prog params from input params
my @Affecteds = map int $_, split / /, $EnteredParams{affected};
foreach my $Child (@Affecteds) {unless ($Child >0) {die "Affected children must be numbered greater than 0\n";}}
my $NrAffecteds = int @Affecteds;
my @Unaffecteds = map int $_, split / /, $EnteredParams{unaffected};
foreach my $Child (@Unaffecteds) {unless ($Child >0) {die "Unaffected children must be numbered greater than 0\n";}}
my $NrUnaffecteds = int @Unaffecteds;
my $NrChildren = $NrAffecteds + $NrUnaffecteds;
if ($NrChildren < 2) {die "Need at least two children to analyse\n";}

my @MatFilesIn = @{$EnteredParams{input_mat}}; # array of mat input files
my $NrMatFiles = int @MatFilesIn;
foreach my $File (@MatFilesIn) {unless (-f $File) {die "Maternal input file $File not found\n";}} # must exist
my @PatFilesIn = @{$EnteredParams{input_pat}}; # array of pat input files
my $NrPatFiles = int @PatFilesIn;
foreach my $File (@PatFilesIn) {unless (-f $File) {die "Paternal input file $File not found\n";}} # rmust exist

# Check correct files present for MOI
my $MOI = lc $EnteredParams{inheritance_mode};
if    ($MOI eq "recessive")         {unless ($NrMatFiles and $NrPatFiles) {die "Need both maternal and paternal state files for recessive and dominant modes of inheritance\n";}}
elsif ($MOI eq "maternal_dominant") {unless ($NrMatFiles) {die "Need maternal state files for maternal-dominant mode of inheritance\n";}}
elsif ($MOI eq "paternal_dominant") {unless ($NrMatFiles) {die "Need paternal state files for paternal-dominant mode of inheritance\n";}}
else								{die "Do not recognise $MOI inheritance mode\nAccepted modes are recessive, maternal_dominant and paternal_dominant\n";}

# Check Dir Out
my $DirectoryOut = $EnteredParams{output_dir}; #
$DirectoryOut =~ s/\/$//; # remove trailing slash if present
unless (-d $DirectoryOut) {die "Output directory $DirectoryOut not found\n";} # requires existing file

# Print processed params
print "Maternal Files:\n\t", join("\n\t",@MatFilesIn),
"\nPaternal Files:\n\t", join("\n\t",@PatFilesIn),
"\nOutput Dir:\n\t$DirectoryOut\n";
print "Header in Output" if $EnteredParams{header};
print "\nMode of Inheritance: $MOI",
"\nAffected Children:\t", join(" ",@Affecteds),
"\nUnaffected Children:\t", join(" ",@Unaffecteds),"\n\n";

# Get list of children IDs from input files, to compare to IDs from input args
my %ChildrenFromFiles; # hash to store nrs of children in input file names and nr occurrences
my (@AffectedMatFiles, @UnaffectedMatFiles, @AffectedPatFiles, @UnaffectedPatFiles);
foreach my $ParentFiles (\@MatFilesIn, \@PatFilesIn)
{
	foreach my $File (@{$ParentFiles})
	{
		#print "\t$File\n";
		my (@AffectedFiles, @UnaffectedFiles); # temp arrays for aff and unaff files
		my $Short = basename $File; # removing path
		my @Children =  $Short =~ /_Ch(\d+)/g; # get array of child nrs from file name
		my $AffectedsCount = my $UnaffectedsCount = 0; # count of children in aff, unaff
		foreach my $Child (@Children) # take each number
		{
			$ChildrenFromFiles{$Child}++; # add to a hash of numbers as key, nr occurrences as val
			foreach my $Child2 (@Affecteds) {if ($Child == $Child2) {$AffectedsCount++;}} # add child to aff count
			foreach my $Child2 (@Unaffecteds) {if ($Child == $Child2) {$UnaffectedsCount++;}} # add child to unaff count
			# any file for 2 unaffected children not used
		}

		if ($AffectedsCount == 2) 								{push @AffectedFiles, $File;} # both children in affecteds, add file to affected file list
		elsif ($AffectedsCount == 1 and $UnaffectedsCount == 1) {push @UnaffectedFiles, $File;} # else one in affected, one in unaffected, add file to unaffected file list
		# if both unaffected, will be processed vs an affected child separately, so not added here

		if ($ParentFiles eq \@MatFilesIn) # processing mat files, so add to mat file arrays
		{
			push @AffectedMatFiles, @AffectedFiles;
			push @UnaffectedMatFiles, @UnaffectedFiles;
		}
		else # processing pat files, so add to pat file arrays
		{
			push @AffectedPatFiles, @AffectedFiles;
			push @UnaffectedPatFiles, @UnaffectedFiles;
		}
	}
}

print "Affected Maternal Files:\n\t",join("\n\t",@AffectedMatFiles),"\n" if int(@AffectedMatFiles);
print "Affected/Unaffected Maternal Files:\n\t",join("\n\t",@UnaffectedMatFiles),"\n" if int(@UnaffectedMatFiles);
print "Affected Paternal Files:\n",join("\n\t",@AffectedPatFiles),"\n" if int(@AffectedPatFiles);
print "Affected/Unaffected Paternal Files:\n\t",join("\n\t",@UnaffectedPatFiles),"\n" if int(@UnaffectedPatFiles);

# Assign files to process, depends on mode of inheritance
my (@AffectedFiles2Use,  @UnaffectedFiles2Use);
if ($MOI eq "maternal_dominant") # mat dom so do not require pat data
{
	@AffectedFiles2Use = @AffectedMatFiles;
	@UnaffectedFiles2Use = @UnaffectedMatFiles;
}
elsif ($MOI eq "paternal_dominant") # pat dom so do not require mat data
{
	@AffectedFiles2Use = @AffectedPatFiles;
	@UnaffectedFiles2Use = @UnaffectedPatFiles;
}
else # recessive, require all
{
	@AffectedFiles2Use = (@AffectedMatFiles, @AffectedPatFiles);
	@UnaffectedFiles2Use = (@UnaffectedMatFiles, @UnaffectedPatFiles);
}

# Check that children listed in input params match children listed in files, and that we have affecteds
my $Match = my $NoMatch = 0; # counters
foreach my $NrIn (@Affecteds, @Unaffecteds) # loop through list of children
{
	if ($ChildrenFromFiles{$NrIn}) {$Match++;} # increment match counter
	else						   {$NoMatch++;} # increment no match counter
}
if ($NrChildren != $Match) {print "Number of children from input paramaters ($NrChildren) does not equal number of children from input files ($Match)\n";}
elsif ($NoMatch)           {print "There are children present in the input files, that are not listed in the input parameters\n";}
unless ($NrAffecteds > 0) {die "Must have at least one affected child\n"}; # require at least one affected

# Header
my $Header;
if ($EnteredParams{header}) # header in params
{
	my $ScriptName = basename ($0); # get script name
	my ($Sec,$Min,$Hour,$Mday,$Mon,$Year,$Wday,$Yday,$Isdst) = localtime(time);
	my @Months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
	$Year += 1900;
    my $DateString = "$Months[$Mon] $Mday $Year\t$Hour:$Min:$Sec";
	$Header .= "# Created By\t$ScriptName\n";
	$Header .= "# When\t$DateString\n";
	$Header .= "# Mode of Inheritance\t$MOI\n";
	$Header .= "# Affected Children\t$EnteredParams{affected}\n";
	$Header .= "# Unaffected Children\t$EnteredParams{unaffected}\n";
	foreach my $File (@AffectedFiles2Use) {$Header .= "# Affected Child File\t".basename($File)."\n";}
	foreach my $File (@UnaffectedFiles2Use) {$Header .= "# Affected/Unaffected Child File\t".basename($File)."\n";}
	$Header .= "\n";
}

# Get stub for output files from an input file
my $FileOut;
if (int(@AffectedFiles2Use)) {($FileOut, undef, undef) = fileparse ($AffectedFiles2Use[0], qr/\.[^.]*/);}
else						 {($FileOut, undef, undef) = fileparse ($UnaffectedFiles2Use[0], qr/\.[^.]*/);}
$FileOut =~ /(^.+?)_Ch?/;
$FileOut = $1;

# Create unique output filenames
print "\nOpening output files:\n";
my $n = 1; my $Suff = "";
my $Prefix = "$DirectoryOut/$FileOut-$MOI-$NrAffecteds"."_Aff";
$Prefix .= "_$NrUnaffecteds"."_Unaff" if $NrUnaffecteds > 0;
$Prefix .= "-Blocks";

while (-f $Prefix.$Suff.".tsv" and -f $Prefix.$Suff.".bed") # loop till name is found that does not exist, ensure we do not overwrite existing file
{$Suff = "_$n"; $n++;} # loop till we have a new unique filename

my $MOIRegionFileOut = $Prefix.$Suff.".tsv";
my $MOIBedFileOut = $Prefix.$Suff.".bed";
print "\t",$MOIRegionFileOut,"\n"; #exit;
print "\t",$MOIBedFileOut,"\n\n"; #exit;
open my $OUT1, ">", $MOIRegionFileOut;
open my $OUT2, ">", $MOIBedFileOut;

# Output Header/Col Headers
my $StateColHeader = ">Chr\tBegin\tEnd\tState\tLength\n";
print $OUT1 $Header; # header for tsv file, if it exists
print $OUT1 $StateColHeader; # col header for tsv file

if ($MOI =~ /paternal/)    {print $OUT2 "track name=Paternal_Overlaps description=Overlap_Blocks itemRgb=on\n";} # col header for bed file
elsif ($MOI =~ /maternal/) {print $OUT2 "track name=Maternal_Overlaps description=Overlap_Blocks itemRgb=on\n";}
else 					   {print $OUT2 "track name=All_Overlaps description=Overlap_Blocks itemRgb=on\n";}

# MAIN #
my $MOIBlocks = NewChrArrayHash (); # hash of chr arrays to store recs
if (int(@AffectedFiles2Use)) # if there are 2 affected children there will be at least 1 affected file
{
	print "Loading inital Affected file\n";
	print "\t$AffectedFiles2Use[0]\n";
	LoadBlocks ($AffectedFiles2Use[0],$MOIBlocks); # load all blocks
	$MOIBlocks = FilterBlocks ($MOIBlocks,"D"); # filters out D blocks, U blocks between D blocks
	print "Loading other Affected file(s)\n" if int(@AffectedFiles2Use) > 0; # only loaded if more than 1
	for my $n (1..int(@AffectedFiles2Use)-1)
	{
		print "\t$AffectedFiles2Use[$n]\n";
		my $MOIBlocks2 = NewChrArrayHash (); # hash of chr arrays to store recs
		LoadBlocks ($AffectedFiles2Use[$n],$MOIBlocks2); # load all blocks
		$MOIBlocks2 = FilterBlocks ($MOIBlocks2, "D"); # filters out D blocks, U blocks between D blocks
		print "\tFinding overlaps\n";
		$MOIBlocks = FindOverlaps ($MOIBlocks, $MOIBlocks2,"S"); # finds overlaps between O|S|U and S|U|G blocks
	}
	print "Loading Affected/Unaffected file(s)\n" if int(@UnaffectedFiles2Use) > 0;
	for my $n (0..int(@UnaffectedFiles2Use)-1)
	{
		print "\t$UnaffectedFiles2Use[$n]\n";
		my $MOIBlocks2 = NewChrArrayHash (); # hash of chr arrays to store recs
		LoadBlocks ($UnaffectedFiles2Use[$n], $MOIBlocks2); # load all blocks
		$MOIBlocks2 = FilterBlocks ($MOIBlocks2,"S"); # filters out S blocks, U blocks between S blocks
		print "\tFinding overlaps\n";
		$MOIBlocks = FindOverlaps ($MOIBlocks,$MOIBlocks2,"D"); # finds overlaps between O|U and D|U|G blocks
	}
	OutputMOIBlocks ($OUT1, $MOIBlocks); #
	OutputMOIBedBlocks ($OUT2, $MOIBlocks); #
}
else # only 1 affected child, so compare just to unaffecteds
{
	print "Loading inital Affected/Unaffected file\n";
	print "\t$UnaffectedFiles2Use[0]\n";
	my $MOIBlocks = NewChrArrayHash (); # hash of chr arrays to store recs
	LoadBlocks ($UnaffectedFiles2Use[0], $MOIBlocks); # load all blocks
	$MOIBlocks = FilterBlocks ($MOIBlocks, "S"); # filters out S blocks, U blocks between S blocks
	print "Loading other Affected/Unaffected file(s)\n" if int(@UnaffectedFiles2Use) > 1;
	for my $n (1..int(@UnaffectedFiles2Use)-1)
	{
		print "\t$UnaffectedFiles2Use[$n]\n";
		my $MOIBlocks2 = NewChrArrayHash (); # hash of chr arrays to store recs
		LoadBlocks ($UnaffectedFiles2Use[$n], $MOIBlocks2); # load all blocks
		$MOIBlocks2 = FilterBlocks ($MOIBlocks2,"S"); # filters out S blocks, U blocks between S blocks
		print "\tFinding overlaps\n";
		$MOIBlocks = FindOverlaps ($MOIBlocks,$MOIBlocks2,"D"); # finds overlaps between O|D|U and D|U|G blocks
	}
	OutputMOIBlocks ($OUT1, $MOIBlocks); # placed here because for some reason the output fails if put outside loops
	OutputMOIBedBlocks ($OUT2, $MOIBlocks); # placed here because for some reason the output fails if put outside loops
}

#OutputMOIBlocks ($OUT1, $MOIBlocks); #
#OutputMOIBedBlocks ($OUT2, $MOIBlocks); #


###########################################################################
#                                   SUBS                                  #
###########################################################################

sub GetExpectedParams # to store parameters and values
{
	my %Hash =
	(
	input_mat => [],
	input_pat => [],
	output_dir => "",
	header => 0,
	inheritance_mode => "",
	affected => -1,
	unaffected => -1,
	);
	return %Hash;
}

sub GetEnteredParams
{
	# Processing @ARGV
	my %Hash;
	my @ARGVs = split /--/, join (" ",@ARGV); # split args on --, into array
	for my $n (1..$#ARGVs) # parse each
	{
		$ARGVs[$n] =~ s/^\s+//; # remove any leading spaces
		$ARGVs[$n] =~ s/\s+$//; # remove any trailing spaces
		my ($Key, $Val) = split /\s+/, $ARGVs[$n], 2; # put first element into key, any other elements into val
		$Key = lc $Key;
		if ($Key eq "input_mat" or $Key eq "input_pat") # into arrays
		{
			push @{$Hash{$Key}}, $Val; # make a hash entry out of key and val
		}
		else
		{
			$Hash{$Key} = $Val; # make a hash entry out of key and val
		}
	}
	return %Hash; # hash now has each -- entry param, with associated values
}

sub TestParameters
{
	# Testing Parameters are filled, all are present, no extras present
	my $NrExpParams = int keys %ExpectedParams;
	my $ParamErrors = 0;
	foreach my $ExpArg (keys %ExpectedParams)
	{
		if ($ExpArg =~/unaffected/) # optional
		{
			next; # increment counter
		}
		if ($ExpArg =~/header/) # optional, give it a value
		{
			$EnteredParams{header} = 1; # increment counter
		}
		elsif ($ExpArg =~/input_/) # mat and pat - need to test later depending on MOI
		{
			next; # increment counter
		}
		elsif (!defined $EnteredParams{$ExpArg})
		{
			print "Parameter --$ExpArg appears to be missing from arguments\n"; $ParamErrors++;
		}
	}
	if ($ParamErrors) {die "$ParamErrors detected in input parameters";}
}

sub SetColours
{
	my %Colours;
	$Colours{Red} = "255,0,0"; # red
	$Colours{Blue} = "0,0,255"; # blue
	$Colours{Green} = "0,255,0"; # green
	$Colours{DarkPurple} = "125,0,125"; # dark purple
	$Colours{LightPurple} = "255,125,255"; # light purple
	$Colours{Purple} = "255,0,255"; # purple
	$Colours{PaleBlue} = "150,150,255"; # pale blue
	$Colours{Grey} = "150,150,150"; # grey
	$Colours{PaleGrey} = "50,50,50"; # pale grey

	$Colours{SameClr} = $Colours{DarkPurple}; # S = red
	$Colours{UncertClr} = $Colours{LightPurple}; # U = grey

	return %Colours;
}

sub LoadBlocks
{
	my $File = shift;
	my $HashPtr = shift; # ptr to hash of chr state rec arrays
	my @Blocks;
	open my $FH, "<", $File;
	<$FH>; # header
	while (<$FH>)
	{
		chomp;
		my @Array = split "\t", $_;
		{
			my $Rec = NewStateRecord ();
			$Rec->{Begin} = $Array[1];
			$Rec->{End} = $Array[2];
			$Rec->{State} = $Array[3];
			push @{$$HashPtr{$Array[0]}}, $Rec;
		}
	}
}

sub FilterBlocks # Removes S or D blocks and U blocks between S or D blocks
{
	my $Blocks = shift; # ptr to hash of chr state rec arrays
	my $NoKeep = shift; # type to discard - S|D

	my $NewBlocks = NewChrArrayHash (); # ptr to hash to store processed recs from rec comparisons

	foreach my $Chr (@ChrNames)
	{
		my $NrRecs = int(@{$$Blocks{$Chr}});
		my $Debug = 0;

		# First entry
		unless ($$Blocks{$Chr}->[0]->{State} eq $NoKeep or # this block is not Keep or
				($$Blocks{$Chr}->[0]->{State} eq "U" and $$Blocks{$Chr}->[1]->{State} eq $NoKeep)) # this block and next block are not U and then Keep
		{
			push @{$$NewBlocks{$Chr}}, $$Blocks{$Chr}->[0]; # so keep this block
		}

		# Entries 2 to n-1
		for my $n (1..$NrRecs-2)
		{
			unless ($$Blocks{$Chr}->[$n]->{State} eq $NoKeep or # this block is NoKeep or
					($$Blocks{$Chr}->[$n-1]->{State} eq $NoKeep and $$Blocks{$Chr}->[$n]->{State} eq "U" and $$Blocks{$Chr}->[$n+1]->{State} eq $NoKeep)) # prev block, this block and next block are not Keep then U then Keep
			{
				push @{$$NewBlocks{$Chr}}, $$Blocks{$Chr}->[$n]; # so keep this block
			}
		}
		# Last entry
		unless ($$Blocks{$Chr}->[$NrRecs-1]->{State} eq $NoKeep or # this block is not Keep or
				($$Blocks{$Chr}->[$NrRecs-2]->{State} eq $NoKeep and $$Blocks{$Chr}->[$NrRecs-1]->{State} eq "U")) # this block and next block are not U then Keep
		{
			push @{$$NewBlocks{$Chr}}, $$Blocks{$Chr}->[$NrRecs-1]; # so keep this block
		}
	}
	return $NewBlocks;
}

sub PrintOne # for debugging
{
	my $Blocks = shift; # ptr to hash of chr state rec arrays
	my $Tot = shift;
	my $Count = 0;
	foreach my $Chr (qw (chr2))
	{
		foreach my $Rec (@{$$Blocks{$Chr}}) # loop through current blocks
		{
			print "$Chr\t"; foreach my $Key (qw (Chr Begin End State)) {print "$Rec->{$Key}\t";}
			print "\n";
			return if $Count++ >= $Tot;
		}
	}
	return;
}

sub FindOverlaps # find recs that overlap, send to sub to find overlap co-ordinates
{
	my $Blocks1 = shift; # ptr to hash of chr state rec arrays
	my $Blocks2 = shift; # ptr to hash of chr state rec arrays for incoming file
	my $Keep = shift; # which type of rec to keep, S|D

	my $NewBlocks = NewChrArrayHash (); # ptr to hash to store processed recs from rec comparisons

	foreach my $Chr (@ChrNames)
	{
		my $Count = 0;
		foreach my $Rec1 (@{$$Blocks1{$Chr}}) # loop through current blocks
		{
			foreach my $Rec2 (@{$$Blocks2{$Chr}}) # loop through current blocks
			{
				unless ($Rec1->{End} < $Rec2->{Begin} or $Rec1->{Begin} > $Rec2->{End}) # compare if overlap
				{
					if ($Rec2->{State} =~ /$Keep|O|U|G/) # overlap, same, uncertain or gap, so find union
					{
						my $Rec = FindIntersection ($Rec1, $Rec2); # Find overlap, retun rec
						push @{$$NewBlocks{$Chr}}, $Rec if $Rec->{Length} > 0; # add new rec to new blocks array if length > 0
					}
				}
			}
		}

	}
	return $NewBlocks;
}

sub FindIntersection # returning range of rec1 that is also Keep|U|G in rec2, assumes there is overlap
{
	my $Rec1 = shift;
	my $Rec2 = shift;

	my $Rec = NewStateRecord ();

	if ($Rec1->{Begin} < $Rec2->{Begin}) {$Rec->{Begin} = $Rec2->{Begin};} # get larger of the Begins
	else								 {$Rec->{Begin} = $Rec1->{Begin};}
	if ($Rec1->{End} > $Rec2->{End})	 {$Rec->{End} = $Rec2->{End};} # get smaler of the Ends
	else								 {$Rec->{End} = $Rec1->{End};}

	$Rec->{Length} = $Rec->{End} - $Rec->{Begin};
	if ($Rec->{Length} > 0) # length > zero, finish and return
	{
		if ($Rec1->{State} =~ /U|G/ or $Rec2->{State} =~ /U|G/) {$Rec->{State} = "U";} # one or both uncertain, so mark as U
		else												  	{$Rec->{State} = "O";} # neithe uncertain, mark as O
	}
	return $Rec;
}

sub OutputMOIBedBlocks
{
	my $FH = shift;
	my $Blocks1 = shift; # ptr to hash of chr state rec arrays

	foreach my $Chr (@ChrNames)
	{
		foreach my $Rec (@{$$Blocks1{$Chr}})
		{
			print $FH  "$Chr\t$Rec->{Begin}\t$Rec->{End}\t$Rec->{State}\t",($Rec->{End}-$Rec->{Begin}),"\t+\t$Rec->{Begin}\t$Rec->{End}\t";
			if ($Rec->{State} eq "O" ) {print $FH  "$Colours{SameClr}\n";}
			else					   {print $FH "$Colours{UncertClr}\n";}
		}
	}
}

sub OutputMOIBlocks
{
	my $FH = shift;
	my $Blocks = shift; # ptr to hash of chr state rec arrays

	foreach my $Chr (@ChrNames)
	{
		foreach my $Rec (@{$$Blocks{$Chr}})
		{
			print $FH  "$Chr\t$Rec->{Begin}\t$Rec->{End}\t$Rec->{State}\t",($Rec->{End}-$Rec->{Begin}),"\n";
		}
	}
}


sub Max
{
	my $Val1 = shift;
	my $Val2 = shift;
	if ($Val1 > $Val2) {return $Val1;}
	else               {return $Val2;}
}

sub Min
{
	my $Val1 = shift;
	my $Val2 = shift;
	if ($Val1 < $Val2) {return $Val1;}
	else               {return $Val2;}
}

sub NewChrArrayHash
{
	my %Hash;
	foreach my $Chr (@ChrNames)
	{
		next if $Chr eq "chrX" or $Chr eq "chrY" or $Chr eq "chrM"; # not processing these yet
		$Hash{$Chr} = [];  # set up ptr to array for each chr, as hash
	}
	return \%Hash;
}

sub NewStateRecord
{
    my $Rec =
    {
		Chr => "",
        Begin => -1,
        End => -1,
		State => "",
        Records => 0,
		Length => -1,
    };
	return $Rec;
}

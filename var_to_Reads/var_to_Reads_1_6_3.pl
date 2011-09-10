#!/usr/bin/perl
use strict;
use File::Find;
$| = 1;

# var to reads file
# takes a nominated chr & range, extracts related information from assembly files
# assembles reads from evidence files if they exist

# Format:
# perl prog dir dir chr range file
# ie
# perl var_to_Reads.pl \
#  --Input_Dir assembly_dir \
#  --Output_Dir dir_path \
#  --Chromosome chr \
#  --Range n1-n2 \
#  --Reference ref.crr
# eg
# perl /rnd/home/rtearle/perl/var_to_Reads_1_6_3.pl \
#  --Input_Dir /GS19240-1100-37-ASM/GS00028-DNA_C01 \
#  --Output_Dir /rtearle/data/ \
#  --Chromosome chr1 \
#  --Range 104206533-104206534 \
#  --Reference /ref/build37.crr

# note that the range separator can be anything as long as it is not a number or a tab
# but it must be of the form integer separator integer

# now covers a range of contiguous positions
# these can encompass more than one interval, intervals aligned and printed separately
# now includes a count of each base type per position, for the two reads
# ruler and numbering should now align correctly
# added reference from cgatools decodecrr

# Flow:
# starts with chr and pos
# gets var entries for this chr and pos from var file
# gets CoverageRefScore entries for this chr and surrounding posns from CoverageRefScore file
# gets evInt entry covering this chr and pos from evInt file
# uses evInt ID to get evDNB entries from evDNB file
# aligns the evDNB reads covering the pos
# counts the bases for each pos

my $Debug = 0;

my $usage= "usage:  perl var_to_Reads_1_6_3.pl input_dir output_dir chr pos1-pos2 ref_crr_file\n";

# Variables from Args - expect input directory, output directory, chr, & pos
my %ExpectedParams =  GetExpectedParams ();
my %EnteredParams = GetEnteredParams ();
TestParameters (%EnteredParams); # no return val used, dies in subroutine if no correct

# Setting up prog paras from input paras
my $FileIn = $EnteredParams{input};


unless (@ARGV >= 4){die "Need main assembly directory, ouput directory, chromosome and position ie four parameters in that order\nA5th parameter, reference crr file, is optional";}
my $DirectoryIn = $EnteredParams{input_dir}; # input dir
unless (-d $DirectoryIn) {die "$!: Main directory $DirectoryIn not found\n";}
$DirectoryIn =~ s/\/$//; # remove trailing /
my $DirectoryOut = $EnteredParams{output_dir}; # output dir
unless (-d $DirectoryOut) {mkdir $DirectoryOut or die "$!: Cannot find/create output directory $DirectoryOut";} # make it if it does not exist or die trying
$DirectoryOut =~ s/\/$//; # remove trailing /
my $Chr = $EnteredParams{chromosome}; # chr
my $Range = $EnteredParams{range}; # position(s)
my ($Pos, $Pos2); # start and end
if ($Range =~ /\D+/) 		{($Pos, $Pos2) = ($Range =~ /^(\d+)\D+(\d+)/)} # if is a range, extract start and end pos
elsif ($Range =~ /^(\d+)/) 	{$Pos = $1;$Pos2 = $Pos + 1;} # if just one pos
else 						{die "Input range is not of the form 'integer' or 'integer separator integer'\n";} # else bad format, exit

my $RefFile = $EnteredParams{reference};
unless (-f $RefFile) {die "$!: Reference file $RefFile not found\n";} # we have a param for ref crr file

# Other Variables
my $IN; # file_in FH
my $Header; # ptr to header lines
my $ColHeader; # column header

# Open temp file in $DirOut
my $FileOut = "$DirectoryOut/var_cov_ev-$Chr-$Pos-$Pos2";
my $n = 1; my $Suff = ".tsv";
while (-f $FileOut.$Suff) {$Suff = "-$n.tsv";$n++;} # loop till name is found that does not exist
my $FileOut2 = $FileOut .= $Suff;
$FileOut2 =~ s/.tsv$/-tabbed.tsv/;
open my $OUT1, ">", $FileOut or die "$!: Cannot open file $FileOut to output data";
open my $OUT2, ">", $FileOut2 or die "$!: Cannot open file $FileOut2 to output data";

# Get paths to ASM and specific subfolders
print "\nFinding and collating data for var at position $Range on $Chr\n";
print "Finding Directory and Subdirectories\n";
my $DirectoryASM = "$DirectoryIn/ASM";
unless (-d $DirectoryASM) {die "$!: Cannot find ASM directory $DirectoryASM\n";}
my $DirectoryREF = "$DirectoryASM/REF";
unless (-d $DirectoryREF) {die "$!: Cannot find REF directory $DirectoryREF\n";}
my $DirectoryEV = "$DirectoryASM/EVIDENCE";
unless (-d $DirectoryEV) {die "$!: Cannot find EVIDENCE directory $DirectoryEV\n";}

# Find path and name of specific files
my ($dbSNPFile, $VarFile, $MasterVarFile, $EvIntFile, $EvDNBFile, $EvCorrFile, $CoverageFile, $GeneFile, $GeneSumFile);
# Find ASM files
opendir ASM, $DirectoryASM or die "$!: Cannot open dir $DirectoryASM";
my @TempFiles = readdir ASM; # get list of dirs and files in ASM
close ASM;
print "Finding Files\n";
foreach my $File (@TempFiles)
{
    if ($File =~ /^dbSNPAnnotated-/) {$dbSNPFile = $File;} # dbSNPFile (not used yet)
    elsif ($File =~ /^var-/)	    {$VarFile = $File;} # var file
    elsif ($File =~ /^masterVar/)	    {$MasterVarFile = $File;} # var file
    elsif ($File =~ /^geneVarSummary-/)    {$GeneSumFile = $File;} # geneVarSummary file
    elsif ($File =~ /^gene-/)      {$GeneFile = $File;} # gene file
}

# Write file Output Header
print $OUT1 "Data for variation(s) at position $Range on $Chr\n";
print $OUT1 "Data for a var entry from coverage, var and evidence files from $DirectoryIn\n";
print $OUT1 "Note that extra entries flanking the first and last desired entries are often returned to give context\n";

#### Get entries from coverageRefScore file ####
print "Finding coverageRefScore File for $Chr \n";
my $CoverageRefScoreFile = FindCoverageRefScoreFile ($DirectoryREF, $Chr);
unless ($CoverageRefScoreFile) {die "$!: cannot find a coverageRefScore File in $DirectoryREF for $Chr";}
unless (-f "$DirectoryREF/$CoverageRefScoreFile") {die "$!: cannot open CoverageRefScore File $DirectoryREF/$CoverageRefScoreFile";}

print "Searching coverageRefScore File $DirectoryREF/$CoverageRefScoreFile for Entries\n";
$IN = OpenFile ("$DirectoryREF/$CoverageRefScoreFile"); # open the file with correct file format
$Header = GetHeader ($IN); # extract header
unless (int @$Header) {die "File $DirectoryASM/$CoverageRefScoreFile does not appear to have a header";}
$ColHeader = GetColHeader ($IN); # extract column header
unless ($ColHeader) {die "File $DirectoryASM/$CoverageRefScoreFile does not appear to have a column header";}
my $Window = 5; # arbitrary choice of nr lines to display above and below range [change this to take into account width of interval(s)?]
my $CoverageRefScores = GetCoverageRefScoreEntries ($IN, $Pos, $Pos2, $Window); # pass fh, posns & window, return CoverageRefScore entries around the range
unless (int @$CoverageRefScores) {die "Cannot find an entry for $Chr position(s) $Range and surrounding positions in $DirectoryREF/$CoverageRefScoreFile\n";} # exit if cannot find entries
PrintHeaders1 ($CoverageRefScoreFile, "CoverageRefScore", $ColHeader); # print to file: col header
PrintData1 ($CoverageRefScores); # print to file: CoverageRefScore entries

#### Get records for entry from var file ####
print "Searching var File for Entries\n";
$IN = OpenFile ("$DirectoryASM/$VarFile"); # open the file with correct file format
$Header = GetHeader ($IN); # extract header
unless (int @$Header) {die "File $DirectoryASM/$VarFile does not appear to have a header";}
my $Version = GetVersion ($Header); # get software version because not all files the same columns
unless ($Version) {die "File $DirectoryASM/$VarFile does not appear to have a version";}
$ColHeader = GetColHeader ($IN); # extract col header
unless ($ColHeader) {die "File $DirectoryASM/$VarFile does not appear to have a column header";}
my $Vars = GetVarEntries ($IN, $Chr, $Pos, $Pos2); # pass fh, chr, posns, return var entries
close $IN;
unless (int @$Vars) {die "Cannot find an entry for $Chr position $Pos in $DirectoryASM/$VarFile\n";} # no entries
PrintHeaders1 ($VarFile, "Variation", $ColHeader); # standard stuff, col header
PrintData1 ($Vars); # leave formatting of output to sub
####################

#### Get records for entry from masterVar file ####
if ($MasterVarFile)
{
	print "Searching masterVar File for Entries\n";
	$IN = OpenFile ("$DirectoryASM/$MasterVarFile"); # open the file with correct file format
	$Header = GetHeader ($IN); # extract header
	unless (int @$Header) {die "File $DirectoryASM/$MasterVarFile does not appear to have a header";}
	$ColHeader = GetColHeader ($IN); # extract col header
	unless ($ColHeader) {die "File $DirectoryASM/$MasterVarFile does not appear to have a column header";}
	my $Vars = GetMasterVarEntries ($IN, $Chr, $Pos, $Pos2); # pass fh, chr, posns, return var entries
	close $IN;
	unless (int @$Vars) {die "Cannot find an entry for $Chr position $Pos in $DirectoryASM/$VarFile\n";} # no entries
	PrintHeaders1 ($VarFile, "Variation (masterVar)", $ColHeader); # standard stuff, col header
	PrintData1 ($Vars); # leave formatting of output to sub
}
else
{
	print "No masterVar File found for this assembly\n";

}

#### Get records for entry from var file ####
print "Searching var File for Entries\n";
$IN = OpenFile ("$DirectoryASM/$GeneFile"); # open the file with correct file format
$Header = GetHeader ($IN); # extract header
unless (int @$Header) {die "File $DirectoryASM/$GeneFile does not appear to have a header";}
$ColHeader = GetColHeader ($IN); # extract col header
unless ($ColHeader) {die "File $DirectoryASM/$GeneFile does not appear to have a column header";}
my $GeneVars = GetGeneVarEntries ($IN, $Chr, $Pos, $Pos2); # pass fh, chr, posns, return var entries
close $IN;
if (int @$GeneVars)
{
	PrintHeaders1 ($GeneFile, "Gene Variation", $ColHeader); # standard stuff, col header
	PrintData1 ($GeneVars); # leave formatting of output to sub
 }
else # no entries
{
	print "No gene entries for $Chr position $Pos in $DirectoryASM/$GeneFile\n";
}

#### Finding evInt, evDNB, correlation files for chr ####
print "Searching for evInt and evDNB Files for $Chr\n";
opendir EV, $DirectoryEV or die "Cannot open dir $DirectoryEV";
@TempFiles = readdir EV;
close EV;
foreach my $File (@TempFiles)
{
    if ($File =~ /^evidenceIntervals-$Chr-/) {$EvIntFile = $File;}
    elsif ($File =~ /^evidenceDnbs-$Chr-/)   {$EvDNBFile = $File;}
    elsif ($File =~ /^correlation-/)         {$EvCorrFile = $File;}
}
unless ($EvIntFile) {die "Cannot find a $Chr evidence interval file in $DirectoryEV\n";}
unless ($EvDNBFile) {die "Cannot find $Chr evidence DNB file in $DirectoryEV\n";}
unless ($EvCorrFile) {die "Cannot find an evidence correlation file in $DirectoryEV\n";}

#### Get record for entry from Evidence Int file ####
# IntervalId Chromosome OffsetInChromosome Length Ploidy AlleleIndexes Score
# Allele0 Allele1 Allele2 Allele1Alignment Allele2Alignment
print "Searching evInt File for Entry(s)\n";
$IN = OpenFile ("$DirectoryEV/$EvIntFile");
$Header = GetHeader ($IN);
unless (int @$Header) {die "File $DirectoryASM/$EvIntFile does not appear to have a header";}
my $EvIntColHeader = GetColHeader ($IN);
unless ($ColHeader) {die "File $DirectoryASM/$EvIntFile does not appear to have a column header";}
my $EvInts = GetEvIntEntries ($IN, $Pos, $Pos2); # pass fh, pos, return evint one entry as ref
close $IN;
unless (int keys %$EvInts) {die "Cannot find any evidence interval entries for $Chr position $Range in $DirectoryEV/$EvIntFile\n";}
PrintHeaders1 ($EvIntFile, "Evidence Interval", $EvIntColHeader); # standard stuff, col header etc
PrintEvInts ($EvInts);

#### Get records for entry from Evidence DNB file ####
print "Searching evDNB File for Entries\n";
$IN = OpenFile ("$DirectoryEV/$EvDNBFile");
$Header = GetHeader ($IN);
unless (int @$Header) {die "File $DirectoryASM/$EvDNBFile does not appear to have a header";}
$ColHeader = GetColHeader ($IN);
unless ($ColHeader) {die "File $DirectoryASM/$EvDNBFile does not appear to have a column header";}
my $Found = GetEvDNBs ($IN, $EvInts); # pass fh, intvl hash ref, dnb entries added to intvl hash
close $IN;
unless ($Found) {die "Cannot find entries for $Chr evidence intervals covering position $Range in $DirectoryEV/$EvDNBFile\n";}
PrintHeaders1 ($EvDNBFile, "Evidence DNBs", ""); # not passing column header, need to do it internally
PrintEvDNBs ($OUT1,$EvInts);

#### Get correlation records for entry from correlation file ####
print "Searching Correlation File for Entries\n";
$IN = OpenFile ("$DirectoryEV/$EvCorrFile");
$Header = GetHeader ($IN);
unless (int @$Header) {die "File $DirectoryASM/$EvCorrFile does not appear to have a header";}
$ColHeader = GetColHeader ($IN);
unless ($ColHeader) {die "File $DirectoryASM/$EvCorrFile does not appear to have a column header";}
my $EvCorrs = [];
$EvCorrs = GetCorrelations ($IN, $EvInts, $Chr); # pass fh, intvl hash ref, return corr hash ref
close $IN;
if (@$EvCorrs)
{
	PrintHeaders1 ($EvCorrFile, "Correlated Intervals", $ColHeader); # not passing column header, need to do it internally
	PrintEvCorrs ($OUT1, $EvCorrs);
}
else
{
	print "There are no correlation intervals for $Chr covering position $Range in $DirectoryEV/$EvCorrFile\n\n";
}

#### Aligning evidence reads for each interval ####
print "Aligning evidence interval alleles\n";
foreach my $ID (sort {$a <=> $b} keys %$EvInts) # sorting keys ie intvl ids
{
    AlignIntervalAlleles ($$EvInts{$ID}); # pass interval rec, align allele seqs
	BuildAlleleRuler ($$EvInts{$ID}); # pass interval rec, format allele seqs
	PrintAlleleRuler ($$EvInts{$ID});
}

#### Extract data to align evidence reads ####
print "Extracting reads from evidence DNB entries\n";
ExtractEvReads ($EvInts, $Version); # pass evintdnb entries, sw version nr, return extracted entries

#### Aligning evidence reads for each interval ####
print "Aligning Interval Reads\n";
foreach my $ID (sort {$a <=> $b} keys %$EvInts) # sorting keys ie intvl ids
{
	BuildReadsRuler ($$EvInts{$ID});
	AlignEvReads ($$EvInts{$ID}, $ID); # pass interval hash ref, aligned reads then in passed hash
	#GetRuler ($$EvInts{$ID}, $Chr); # pass ptr to ruler array, return formatted ruler
}
############################


#### Count Bases per Position ####
print "Counting number of bases per position across Evidence Intervals\n";
foreach my $ID (sort {$a <=> $b} keys %$EvInts) # sorting keys ie intvl ids
{
    GetEvReadBaseCounts ($EvInts, $ID);
}
############################

#### Count Bases per Position ####
print "Counting number of each base per position across Evidence Intervals\n";
foreach my $ID (sort {$a <=> $b} keys %$EvInts) # sorting keys ie intvl ids
{
    GetEvReadBases ($$EvInts{$ID});
}
############################

##### Printing alignments and counts for each interval #####
print "Printing Alignments and Counts\n";
foreach my $ID (sort {$a <=> $b} keys %$EvInts) # sorting keys ie intvl ids
{
    PrintEvAlignmentHeader($OUT1); # pass fh
    PrintEvAlignments ($OUT1, $ID, $$EvInts{$ID}); # pass fh, intvl id, ref to intvl, ruler
    PrintTabbedEvAlignments ($OUT2, $ID, $$EvInts{$ID}); # pass fh, intvl id, ref to intvl, ruler
    PrintEvReadBaseCountHeader ($OUT1); # pass fh
    PrintEvReadBaseCounts ($OUT1, $ID, $$EvInts{$ID}); # pass pass fh, intvl id, ref to intvl
}
############################

#### Clean Up ####
print "Completed\n";
############################

###########################
####        Subs       ####
###########################
sub GetExpectedParams
{
	my %Hash =
	(
	"input_dir" => -1,
	"output_dir" => -1,
	"chromosome" => -1,
	"range" => -1,
	"reference" => -1,
	); # to store parameters and values
	return %Hash;
}

sub GetEnteredParams
{
	# Processing @ARGV
	my %Hash;

	my @ARGVs = split /--/, join (" ",@ARGV); # split args on --, into array
	#print "Start\n", join ("\n",@ARGVs),"\n",int @ARGVs - 1,"\n\n" if $Debug;
	#print "Key\tVal\n" if $Debug; #exit;
	for my $n (1..$#ARGVs) # parse each
	{
		$ARGVs[$n] =~ s/\s+$//; # remove any trailing spaces
		my ($Key, $Val) = split / /, $ARGVs[$n], 2; # put first element into key, any other elements into val
		$Key = lc $Key;
		$Hash{$Key} = $Val; # make a hash entry out of key and val
	}
	#foreach my $Key (keys %Hash) {print "$Key\t$Hash{$Key}\n";} exit;
	return %Hash; # hash now has each -- entry param, with associated values
}

sub TestParameters
{
	# Testing Parameters are filled, all are present, no extras present
	my $ArgErrors = 0;
	foreach my $ExpArg (keys %ExpectedParams)
	{
		if (!defined $EnteredParams{$ExpArg})
			{print "Parameter --$ExpArg appears to be missing from arguments\n"; $ArgErrors++;}
		elsif ($EnteredParams{$ExpArg} == -1)
			{print "No value for parameter --$ExpArg in arguments"; $ArgErrors++;}
	}
	if (int keys %EnteredParams > 6) {die "There are too many parameter in the arguments";}
	if ($ArgErrors) {die "$ArgErrors detected in input parameters";}
}

sub GetdbSNPEntry # curr not used
{
    my $dbSNPRec; # record to store nominated dbSNP entry in
    my $FoundEntry = 0;

}

sub GetVarEntries
{
	# >locus ploidy haplotype chromosome begin end varType reference alleleSeq
	# totalScore hapLink xRef
    my $FH = shift;
    my $Chr = shift;
    my $Pos = shift;
    my $Pos2 = shift;
    my @Array = (); # record to store entries for nominated posn
    my $Leave = 0; # used to get out after all of contiguous pos entries have been found
    while (<$FH>)
    {
	if (/$Chr\t/)
	{
	    (undef, undef, undef, undef, my $P1, my $P2, undef) = split /\t/, $_, 7;
	    unless ($Pos > $P2+1 or $Pos2 < $P1-1) # searching for entries at the nom pos
	    {
		chomp $_;
		push @Array, $_;
		$Leave = 1; # found records
	    }
	}
	else
	{
	    last if $Leave;
	}
    }
    return \@Array; #

}

sub GetMasterVarEntries
{
	# >locus ploidy chromosome begin end zygosity varType reference allele1Seq allele2Seq allele1Score allele2Score
	# allele1HapLink allele2HapLink xRef evidenceIntervalId allele1ReadCount allele2ReadCount neitherAlleleReadCount totalReadCount
	# allele1Gene allele2Gene miRBaseId repeatMasker segDupOverlap relativeCoverag calledPloidy
    my $FH = shift;
    my $Chr = shift;
    my $Pos = shift;
    my $Pos2 = shift;
    my @Array = (); # record to store entries for nominated posn
    my $Leave = 0; # used to get out after all of contiguous pos entries have been found
    while (<$FH>)
    {
		if (/$Chr\t/)
		{
			(undef, undef, undef, my $P1, my $P2, undef) = split /\t/, $_, 6;
			unless ($Pos > $P2+1 or $Pos2 < $P1-1) # searching for entries at the nom pos
			{
				chomp $_;
				push @Array, $_;
				$Leave = 1; # found records
			}
		}
		else
		{
			last if $Leave;
		}
    }
    return \@Array; #

}

sub GetGeneVarEntries
{
	# >index locus allele chromosome begin end varType reference call xRef geneId mrnaAcc proteinAcc symbol orientation component componentIndex codingRegionKnown impact nucleotidePos proteinPos annotationRefSequence sampleSequence genomeRefSequence
    my $FH = shift;
    my $Chr = shift;
    my $Pos = shift;
    my $Pos2 = shift;
    my @Array = (); # record to store entries for nominated posn
    my $Leave = 0; # used to get out after all of contiguous pos entries have been found
    while (<$FH>)
    {
	if (/$Chr\t/)
	{
	    (undef, undef, undef, undef, my $P1, my $P2, undef) = split /\t/, $_, 7;
	    unless ($Pos > $P2+1 or $Pos2 < $P1-1) # searching for entries at the nom pos
	    {
		chomp $_;
		push @Array, $_;
		$Leave = 1; # found records
	    }
	}
	else
	{
	    last if $Leave;
	}
    }
    return \@Array; #

}

sub GetCoverageRefScoreEntries
{
    my $FH = shift;
    my $Start = shift;
    my $Stop = shift;
    my $Window = shift;
    my @Array;
    if ($Window < 0) # ignoring -ve values, reset to 1, honouring 0 or above (need to build in to user input?)
    {
	print "Window value of $Window for finding Reference Coverage entries will be ignored\n";
	$Window = 1;
    }
    $Start = Max2Ints($Start-$Window, 0);
    $Stop += $Window;
    #print "$Start\t$Stop\n";
    my $Pos;
    my $Leave = 0; # used to get out after all of contiguous matching entries have been found
    while (<$FH>) # now finding remainder of posns matching range from $Start to $Sop
    {
	/^(\d+)(?{$Pos=$^N})\t/; # seems to be a memory killer?
	if ($Pos >= $Start and $Pos <= $Stop)
	{
	    #print "$Start\t$Stop\t$Pos\n";
	    chomp;
	    push @Array, $_; # chomp it, add to array
	    $Leave = 1; # found records
	}
	else # out of the block, leave
	{
	    last if $Leave;
	}
    }
    return \@Array;
}

sub GetEvIntEntries
{
	# IntervalId Chromosome OffsetInChromosome Length Ploidy AlleleIndexes Score
	# Allele0 Allele1 Allele2 Allele1Alignment Allele2Alignment
    my $FH = shift;
    my $Pos = shift;
    my $Pos2 = shift;
    my %IntHash = (); # storing evint recs, number depends on range inputted
    my $Found = 0; # used to get out after all of contiguous matching entries have been found
    while (<$FH>)
    {
		my ($ID, undef, $P1, $P2, undef) = split /\t/, $_, 5; # using offset and length
		$P2 += $P1;
		unless ($Pos > $P2+1 or $Pos2 < $P1-1) # searching for entries at the nom pos
		{
			my $P = $IntHash{$ID} = NewIntervalRec();
			chomp;
			$P->{IntEntry} = $_;
			$P->{IntAlleleCIGARs}[0] = ($P2-$P1).'M';
			(undef, $P->{Chr}, undef, undef, undef, $P->{Alleles}, undef, $P->{IntAlleleSeqs}[0], $P->{IntAlleleSeqs}[1], $P->{IntAlleleSeqs}[2], $P->{IntAlleleCIGARs}[1], $P->{IntAlleleCIGARs}[2]) = split /\t/, $_;
			$P->{Begin} = $P1;
			$P->{End} = $P2;
			$Found = 1; # found records
			#print "$P->{IntEntry}\n";
			#print "$P->{Chr}\t$P->{Begin}\t$P->{End}\t$P->{Alleles}\t";
			#print join("\t",@{$P->{IntAlleleSeqs}}),"\t";
			#print join("\t",@{$P->{IntAlleleCIGARs}}),"\n";
		}
		else # out of the block, leave
		{
			last if $Found;
		}
    }
    return (\%IntHash);
}

sub GetCorrelations
{
    my $FH = shift;
    my $IntHash = shift; # ref to hash
    my $Chr = shift;
    my $Found = 0;
    my @EvCorrs = (); # storing corr recs
    my @Temp;
    foreach my $IntID (sort {$a <=> $b} keys %$IntHash) # sorting keys ie intvl ids
    {
		push @Temp, $IntID; # stored sorted intvl ids temporarily
    }
    #my $P1 = $Temp[0]; # get first ID in array
    #my $P2 = $Temp[$#Temp]; # get last ID in array, all IDs should be contiguous
    #print "$Temp[0]\t$Temp[$#Temp]\n";
    #print "$$IntHash{$Temp[0]}->{IntEntry}\n";
    #print "$$IntHash{$Temp[$#Temp]}->{IntEntry}\n";
    my (undef, undef, $First, undef) = split /\t/, $$IntHash{$Temp[0]}->{IntEntry},4;
    my (undef, undef, $Last, undef) = split /\t/, $$IntHash{$Temp[$#Temp]}->{IntEntry},4;
    #print "$First\t$Last\n"; #exit;

    while (<$FH>)
    {
		my ($Chr0, $Start0, undef, $Chr1, $Start1, undef) = split /\t/, $_, 6; # using offset and length
		#if ($Pos > $P2+1 or $Pos2 < $P1-1) # searching for entries at the nom pos
		if ($Chr eq $Chr0 or $Chr eq $Chr1)
		{
			if (($First >= $Start0 and $Last <= $Start0) or ($First >= $Start1 and $Last <= $Start1))
			{
				chomp;
				push @EvCorrs, $_;
				$Found = 1; # found records
			}
		}
		else # out of the block, leave
		{
			last if $Found;
		}
    }
    return (\@EvCorrs);
}

sub PrintEvInts
{
    my $HRef = shift;
    foreach my $ID (sort {$a <=> $b} keys %$HRef) # sorting keys ie intvl ids
    {
		print $OUT1 "$$HRef{$ID}->{IntEntry}\n";
    }
}


#sub PrintAlignedIntervalAlleles
#{
#    my $HRef = shift;
#    foreach my $ID (sort {$a <=> $b} keys %$HRef) # sorting keys ie intvl ids
#    {
#		my $HP = $$HRef{$ID};
#		AlignIntervalAlleles ($HP);
#		GetAlleleRulers ($HP, $Chr);
#		for my $a (0..2) # for each allele
#		{
#			next unless $HP->{Alleles} =~ /$a/; # don't process if no data for this allele
#			print $OUT1 "Allele $a of $HP->{Alleles} for EvInt $ID\n";
#			print $OUT1 $HP->{Ruler}[$n]->{ScaleNrs},"\n"; # ptr to scale nrs array
#			print $OUT1 $HP->{Ruler}[$n]->{ScaleLine},"\n"; # ptr to scale nrs array
#			print $OUT1 $HP->{Ruler}[$n]->{ScaleSeq},"\n"; # ptr to scale nrs array
#			#print join("\n",@{$HP->{Ruler}[$a]->{RefNrs}}),"\n";
#		}
#    }
#}

sub AlignIntervalAlleles
{
    my $HP = shift; # pointer to ID rec ref
	my @Posns; # array of alleles (1 & 2) to store posns
	for my $a (0..2) {$Posns[$a] = [];} # array to store insertion posns and lengths for allele
	my @CIGAR; # array of CIGAR values from CIGAR string
	#(undef, undef, my $Start, my $Len, undef, undef, undef, undef, undef,
	# undef, $CIGAR[1], $CIGAR[2]) = split "\t",  $HP->{IntEntry};
	#$CIGAR[0] = $Len.'M';

	# Set up ruler records for each allele
	for my $n (0..2) {$HP->{AlleleRuler}[$n] = NewRulerRecord ();} # recs to store ruler

	# Interval
	for my $a (0..2) # for each allele
	{
		next unless ($HP->{Alleles} =~ /$a/ or $a == 0); # don't process if no data for this allele
		#print "allele: $a\talleles: $HP->{Alleles}\n";
		my @CigCom = ($HP->{IntAlleleCIGARs}[$a] =~ m/(\D+)/g); # CIGAR commands
		my @CigLen = ($HP->{IntAlleleCIGARs}[$a] =~ m/(\d+)/g); # CIGAR lengths used for offsets for base counting
		$HP->{AlleleRuler}[$a]->{CIGAR} = NewCIGARRecord (); # attach com array to intvl rec
		$HP->{AlleleRuler}[$a]->{CIGAR}->{CigCom} = \@CigCom; # attach com array to intvl rec
		$HP->{AlleleRuler}[$a]->{CIGAR}->{CigLen} = \@CigLen; # attach len array to intvl rec
		my @Seq = split //, $HP->{IntAlleleSeqs}[$a]; # allele sequence as array
		my $s = 0;
		#print "$a\tPosns",$Posns[$a],"\tRuler",$HP->{AlleleRuler}[$a]->{RefNrs},"\n";
		#print "CIGAR for $a\t",$HP->{IntAlleleCIGARs}[$a],"\n";
		my $OffSet = $HP->{Begin};
		my $p = 0; # position wrt reference
		for my $c (0..int(@CigCom)-1) # loop through CIGARs for this allele
		{
			if ($CigCom[$c] =~ "M")
			{
				for my $i (1..$CigLen[$c])
				{
					$Posns[$a][$p] = $OffSet;
					$HP->{AlleleRuler}[$a]->{AlleleSeq} .= $Seq[$s];
					$p++; $OffSet++;$s++;
				}
			}
			if ($CigCom[$c] =~ "N|P")
			{
				for my $i (1..$CigLen[$c])
				{
					$Posns[$a][$p] = $OffSet;
					$HP->{AlleleRuler}[$a]->{AlleleSeq} .= " ";
					$p++; $OffSet++; $s++;
				}
			}
			if ($CigCom[$c] =~ "D")
			{
				for my $i (1..$CigLen[$c])
				{
					$Posns[$a][$p] = $OffSet;
					$HP->{AlleleRuler}[$a]->{AlleleSeq} .= "-";
					$p++; $OffSet++;
				}
			}
			elsif ($CigCom[$c] eq "I") # special case, interpolate posns
			{
				$OffSet-- if $OffSet > 0;
				for my $i (1..$CigLen[$c])
				{
					$Posns[$a][$p] = $OffSet + $i/100;
					$HP->{AlleleRuler}[$a]->{AlleleSeq} .= $Seq[$s];
					$p++; $s++;
				}
				$OffSet++;
			}
		}
		#print "$HP->{IntAlleleSeqs}[$a]\n";
		#print "$HP->{AlleleRuler}[$a]->{AlleleSeq}\n";
	}
	for my $a (0..2)
	{
		next unless ($HP->{Alleles} =~ /$a/ or $a == 0); # don't process if no data for this allele
		@{$HP->{AlleleRuler}[$a]->{Positions}} = @{$Posns[$a]}; # attach interval
	}
	return;
}

sub BuildAlleleRuler
{
    my $HP = shift; # pointer to ID rec ref
	my @Ref;

	my $ScaleStart = 10 * int($HP->{Begin}/10); # get first pos divisible by 10
	#print "ScaleStart: $ScaleStart\n";
		@Ref = split //, `cgatools decodecrr --reference $RefFile --range $HP->{Chr}:$HP->{Begin}-$HP->{End}`;
		#print join("",@Ref); exit;
		pop @Ref; # remove trailing return
	#print join("",@Ref),"\n"; #exit;
	# Set up numbering, using refpos data
	for my $a (0..2) # each allele
	{
		next unless ($HP->{Alleles} =~ /$a/ or $a == 0); # don't process if no data for this allele
		my $P = $HP->{AlleleRuler}[$a]; # ptr to ruler rec for this allele
		my @Legend; # numbering positions
		my $Scale = ""; # scale
		my $RefSeq = ""; # reference sequence
		my $AlleleSeq = ""; # allele sequence
		my $r = 0; # counter for ref base
		my $s = 0; # counter for allele base
		for my $m (0..int@{$P->{Positions}}-1)
		{
			#print "$P->{Positions}[$m]\n"; exit;
			if ($P->{Positions}[$m] == int $P->{Positions}[$m]) # integer ie ref
			{
				$RefSeq .= $Ref[$r++]; # add ref seq base, move counter to next base
				if ($P->{Positions}[$m] % 10 == 0) # integer = 0
				{
					$Scale .= "|"; # insert marker
					if ($m >= 2) # room to insert nr length 3
					{
						$Legend[$m] = substr($P->{Positions}[$m],-1,1); # 1st nr
						$Legend[$m-1] = substr($P->{Positions}[$m],-2,1); # 2nd nr
						$Legend[$m-2] = substr($P->{Positions}[$m],-3,1); # 3rd nr;
					}
					else # no room to insert nr length 3
					{
						$Legend[$m] = " "; # nr spacer
					}
				}
				else # integer ­ 0
				{
					$Scale .= "."; # line spacer
					$Legend[$m] = " "; # nr spacer
				}
			}
			else # non-integer ie insertion
			{
				$RefSeq .= " "; # add ref seq base
				$Legend[$m] = " "; # nr spacer
				$Scale .= "_"; # not ref, so empty spacer
			}
		}
		foreach my $Val (@Legend) {$HP->{AlleleRuler}[$a]->{Legend} .= $Val;}
		$HP->{AlleleRuler}[$a]->{Scale} = $Scale;
		$HP->{AlleleRuler}[$a]->{RefSeq} = $RefSeq;
	}

}

sub PrintAlleleRuler
{
    my $HP = shift;
	print $OUT1 "Aligned interval allele sequences\n\n";
	for my $a (0..2) # for each allele
	{
		next unless ($HP->{Alleles} =~ /$a/ or $a == 0); # don't process if no data for this allele
		print $OUT1 "Allele $a:\n";
		print $OUT1 $HP->{AlleleRuler}[$a]->{Legend},"\n"; # ptr to scale nrs
		print $OUT1 $HP->{AlleleRuler}[$a]->{Scale},"\n"; # ptr to scale line
		print $OUT1 $HP->{AlleleRuler}[$a]->{RefSeq},"\n"; # ptr to scale ruler
		print $OUT1 $HP->{AlleleRuler}[$a]->{AlleleSeq},"\n\n"; # ptr to scale ruler
	}
}

sub GetEvDNBs
{
    my $FH = shift;
    my $HP = shift; # ref to hash of evint ids
    my @Temp;
    foreach my $IntID (sort {$a <=> $b} keys %$HP) # sorting keys ie intvl ids
    {
		push @Temp, $IntID; # stored sorted intvl ids temporarily
    }
    my $P1 = $Temp[0]; # get first ID in array
    my $P2 = $Temp[$#Temp]; # get last ID in array, all IDs should be contiguous
    my $Found = 0;
    while (<$FH>)
    {
		my ($ID, undef) = split /\t/, $_, 2;
		if ($ID >= $P1 and $ID <= $P2) # searching for entries with these IDs
		{
			chomp;
			push @{$$HP{$ID}->{DNBs}}, $_; # add evDNB entry to array of DNBs
			my (undef, undef, undef, undef, undef, undef, $Allele, undef,
				undef, undef, undef, $Offset, undef) = split /\t/, $_, 13; # extract allele and pos
			$$HP{$ID}->{DNBCounts}[$Allele]++; # increment DNB count for allele
			$$HP{$ID}->{MinOffset} = $Offset if $$HP{$ID}->{MinOffset} > $Offset; # get min pos
			$$HP{$ID}->{MaxOffset} = $Offset + 35 if $$HP{$ID}->{MaxOffset} < $Offset + 35; # get max pos
			$Found = 1;
		}
		else
		{
			last if $Found;
		}
    }
    return $Found; # to signal that some were found [dnbs are in the hash arrays]
}

sub PrintEvDNBs
{
    my $FH = shift;
    my $HRef = shift;
    foreach my $ID (sort {$a <=> $b} keys %$HRef) # sorting keys ie intvl ids
    {
	print $FH "Evidence DNBs for Interval $ID\n";
	for my $n (0..int(@{$$HRef{$ID}->{DNBCounts}})-1)
	{
	    print $FH "$$HRef{$ID}->{DNBCounts}[$n] reads for allele $n\n"; # if $$HRef{$ID}->{DNBs}[$n];
	}
	print $FH "\n";
	print $FH $ColHeader;
	for my $n (0..int(@{$$HRef{$ID}->{DNBs}})-1)
	{
	    print $FH $$HRef{$ID}->{DNBs}[$n]."\n";
	}
    }
}

sub PrintEvCorrs
{
    my $FH = shift;
    my $ARef = shift;
    #print $FH "Correlated Evidence Entries for Intervals\n";
    foreach my $Corr (@$ARef) # sorting keys ie intvl ids
    {
		print $FH "$Corr\n"; # could just do a join?
    }
	print $FH "\n";
}

sub ExtractEvReads
{
    my $HRef = shift; # ref to hash
    my $Version = shift; # sw version nr, need to choose how to split intdnb entries
    my ($V1, $V2, $V3, $V4) = split  /\./, $Version;
    foreach my $ID (keys %$HRef) # starting with each intvl hash
    {
	my $Tot = @{$$HRef{$ID}->{DNBs}};
	for my $n (0...$Tot-1)
	{
	    my $P = $$HRef{$ID}->{Reads}[$n] = NewReadRecord ();
	    # intDNB records differ between versions, so testing for version and splitting record accordingly
	    if ($V2 == 7 and $V3 == 0) # ie 1.7.0.x
	    {
			($P->{ID}, undef, undef, undef, undef, $P->{Allele}, $P->{Arm}, $P->{Strand}, $P->{Offset}, undef, $P->{OffRef}, $P->{Alignment}, undef, undef, undef, undef, undef, undef, $P->{Seq}, $P->{Qual}, undef) = split ( /\t/, $$HRef{$ID}->{DNBs}[$n]);
	    }
	    else # later version
	    {
			($P->{ID}, undef, undef, undef, undef, undef, $P->{Allele}, $P->{Arm}, $P->{Strand}, $P->{Offset}, undef, $P->{OffRef}, $P->{Alignment}, undef, undef, undef, undef, undef, undef, $P->{Seq}, $P->{Qual}, undef) = split ( /\t/, $$HRef{$ID}->{DNBs}[$n]);
	    }

	    if ($P->{Arm} eq "L")
	    {
		$P->{Seq} = substr($P->{Seq}, 0, 35);
		$P->{Qual} = substr($P->{Qual}, 0, 35);
	    }
	    else
	    {
		$P->{Seq} = substr($P->{Seq}, 35, 35);
		$P->{Qual} = substr($P->{Qual}, 0, 35);
	    }

	    if ($P->{Strand} eq "-")
	    {
		$P->{Seq} = reverse $P->{Seq};
		$P->{Seq} =~ tr/AGCT/TCGA/;
		$P->{Qual} = reverse $P->{Qual};
	    }

	    @{$P->{CigCom}} = ($P->{Alignment} =~ m/(\D+)/g); # CIGAR commands
	    @{$P->{CigLen}} = ($P->{Alignment} =~ m/(\d+)/g); # CIGAR lengths used for offsets for base counting
	    $$HRef{$ID}->{MinOffset} = $P->{OffRef} if $P->{OffRef} < $$HRef{$ID}->{MinOffset}; # Calculating now because needed in next part
	}
    }
    return; #reads attached to hash records
}

sub BuildReadsRuler
{
    my $HP = shift; # pointer to ID rec ref

	# build array of ref posns, including insertion posns
	for my $a (0..2)
	{
		next unless ($HP->{Alleles} =~ /$a/ or $a == 0); # don't process if no data for this allele ie 2, as 0 and 1 must exist
		for my $n ($HP->{MinOffset}..$HP->{Begin}-1) {push @{$HP->{Ruler}[$a]->{Positions}}, $n;} # attach pre interval
		push @{$HP->{Ruler}[$a]->{Positions}}, @{$HP->{AlleleRuler}[$a]->{Positions}}; # attach interval from allele ruler
		for my $n ($HP->{End}..$HP->{MaxOffset}) {push @{$HP->{Ruler}[$a]->{Positions}}, $n;} # attach post interval
		#print "allele $a\n";
	}
	#for my $a (0..2)
	#{
	#	next unless ($HP->{Alleles} =~ /$a/ or $a == 0); # don't process if no data for this allele
	#	print "allele $a\n";
	#	print join("\t",@{$HP->{Ruler}[$a]->{Positions}}),"\n";
	#	print "\n\n";
	#}
	#exit;
	# build legend, scale, reference seq
	my @Ref; # get and store reference sequence
		@Ref = split //, `cgatools decodecrr --reference $RefFile --range $HP->{Chr}:$HP->{MinOffset}-$HP->{MaxOffset}`;
		pop @Ref; # remove trailing return
	#print join("",@Ref),"\n\n";


	# create
	for my $a (0..2)
	{
	my @Legend; # numbering positions
	my $Scale = ""; # scale
	my $RefSeq = ""; # reference sequence
	my $r = my $s = 0; # counters for ref base and allele base

		next unless ($HP->{Alleles} =~ /$a/ or $a == 0); # don't process if no data for this allele
		my $P = $HP->{Ruler}[$a]; # pointer to allele specific struct
		for my $m (0..int@{$P->{Positions}}-1) # loop through positions for this allele
		{
			if ($P->{Positions}[$m] == int $P->{Positions}[$m]) # integer ie ref
			{
				$RefSeq .= $Ref[$r++]; # add ref seq base, move counter to next base
				if ($P->{Positions}[$m] % 10 == 0) # integer = 0, insert pos in legend
				{
					$Scale .= "|"; # insert marker in scale
					if ($m >= 2) # room to insert nr length 3
					{
						$Legend[$m] = substr($P->{Positions}[$m],-1,1); # 1st nr in legend
						$Legend[$m-1] = substr($P->{Positions}[$m],-2,1); # 2nd in to legend
						$Legend[$m-2] = substr($P->{Positions}[$m],-3,1); # 3rd in to legend
					}
					else # no room to insert nr length 3
					{
						$Legend[$m] = " "; # nr spacer in legend
					}
				}
				else # integer ­ 0
				{
					$Scale .= "."; # line spacer in scale
					$Legend[$m] = " "; # nr spacer in legend
				}
			}
			else # non-integer ie insertion
			{
				$Legend[$m] = " "; # nr spacer
				$Scale .= "_"; # not ref, so empty spacer
				$RefSeq .= " "; # add space to ref seq
			}
		}
		foreach my $Val (@Legend) {$HP->{Ruler}[$a]->{Legend} .= $Val;}
		$HP->{Ruler}[$a]->{Scale} = $Scale;
		$HP->{Ruler}[$a]->{RefSeq} = $RefSeq;
		#print "$HP->{Ruler}[$a]->{Legend}\n";
		#print "$HP->{Ruler}[$a]->{Scale}\n";
		#print "$HP->{Ruler}[$a]->{RefSeq}\n";
	}
}


sub AlignEvReads
{
    my $HP = shift; # reads hash ref
    my $ID = shift;

    my $Tot = @{$HP->{Reads}};
	my %Posns; # hash to store positions, including insertions
	#my $Posns; # array of alleles (1 & 2) to store insertion posns and lengths
	# Check for insertions
	#if ($HP->{IntEntry} =~ /(\d{1,2}I)/)
	#{
	#	#print "$HP->{IntEntry}\n";
	#	$Posns = SetPositions ($HRef, $ID); # need to adjust ruler, storage for insertion
	#}
    # Loop through entries to align reads using CIGAR
    for my $r (0..$Tot-1) # loop through reads
    {
		my $ReadPos = 0; # counter for pos in read
		my $RefPos = 0; # counter for pos in "ref"
		my $Curr = 0; # counter for formatting the string
		my $Line = 1; # start with line 1 ,switch o line 2 after B statement
		my $P = $HP->{Reads}[$r]; # pointer to reads array ie intvl -> reads []
		#$Debug = 1 if $P->{Allele} > 0;
		#exit if $P->{Allele} > 1;
		my $PadLen = $P->{OffRef} - $HP->{MinOffset}; # length to pad at start to align read
		#if ($P->{Alignment} =~ /(\d{1,2}I)/) {$Debug = 0; print "$P->{Alignment}\t$P->{Seq}\n";} else {$Debug = 0;}
		print "Min: $HP->{MinOffset}\tOffRef: $P->{OffRef}\tPadLen: $PadLen\n" if $Debug;
		my (@FormattedRead1, @FormattedRead2);
		my (@FormattedScores1, @FormattedScores2);
		for my $n (1..$PadLen) # pad start of formatted read
		{
			$FormattedRead1[$RefPos] = " ";
			$FormattedScores1[$RefPos] = " ";
			print "$RefPos\t$FormattedRead1[$RefPos]\t$FormattedScores1[$RefPos]\n" if $Debug; #exit;
			$RefPos++;
		}
		print "Here: RefPos $RefPos\n" if $Debug;
		for my $c (0..$#{$P->{CigCom}}) # loop through CIGAR elements
		{
			if ($P->{CigCom}[$c] =~ /I/) # ins, count moving only read counter forward
			{
				print "$P->{CigCom}[$c]\t$P->{CigLen}[$c]\n" if $Debug;
				for my $i (1..$P->{CigLen}[$c])
				{
					if ($Line == 1)
					{
						$FormattedRead1[$RefPos] = substr ($P->{Seq},$ReadPos, 1);
						$FormattedScores1[$RefPos] = substr ($P->{Qual},$ReadPos, 1);
						#print "$RefPos\t$FormattedRead1[$RefPos]\t$FormattedScores1[$RefPos]\n" if $Debug; exit;
					}
					else
					{
						$FormattedRead2[$RefPos] = substr ($P->{Seq},$ReadPos, 1);
						$FormattedScores2[$RefPos] = substr ($P->{Qual},$ReadPos, 1);
						#print "$RefPos\t$FormattedRead2[$RefPos]\t$FormattedScores2[$RefPos]\n" if $Debug; exit;
					}
					$ReadPos++;
				}
			}
			elsif ($P->{CigCom}[$c] =~ /M/) # match, count moving both counters forward
			{
				print "$P->{CigCom}[$c]\t$P->{CigLen}[$c]\n" if $Debug;
				for my $i (1..$P->{CigLen}[$c])
				{
					if ($Line == 1)
					{
						$FormattedRead1[$RefPos] = substr ($P->{Seq},$ReadPos, 1);
						$FormattedScores1[$RefPos] = substr ($P->{Qual},$ReadPos, 1);
						print "$RefPos\t$FormattedRead1[$RefPos]\t$FormattedScores1[$RefPos]\n" if $Debug; #exit;
					}
					else
					{
						$FormattedRead2[$RefPos] = substr ($P->{Seq},$ReadPos, 1);
						$FormattedScores2[$RefPos] = substr ($P->{Qual},$ReadPos, 1);
						print "$RefPos\t$FormattedRead2[$RefPos]\t$FormattedScores2[$RefPos]\n" if $Debug; #exit;
					}
					$RefPos++; $ReadPos++;
				}
			}
			elsif ($P->{CigCom}[$c] eq "B") # overlap, move ref counter back, read counter not changed
			{
				$Line = 2; # moving to new line after B
				print "\nLine $Line\n$P->{CigCom}[$c]\t$P->{CigLen}[$c]\n" if $Debug;
				$PadLen = $RefPos - $P->{CigLen}[$c] - 1; # allow for backspace
				if ($P->{CigCom}[$c-1] eq "D") {$PadLen -= $P->{CigLen}[$c-1];} # corr for del in overlap, increase backspace by prev D
				print "PadLen $PadLen\n" if $Debug;
				$RefPos = 0;
				for my $n (0..$PadLen) # pad start of formatted read, resetting counter refpos, inlcuding backspace etc
				{
					$FormattedRead2[$RefPos] = " ";
					$FormattedScores2[$RefPos] = " ";
					print "$RefPos\t$FormattedRead2[$RefPos]\t$FormattedScores2[$RefPos]\n" if $Debug; #exit;
					$RefPos++;
				}
			}
			elsif ($P->{CigCom}[$c] eq "N") # gap, move ref counter forward, read counter not changed
			{
				print "$P->{CigCom}[$c]\t$P->{CigLen}[$c]\n" if $Debug;
				for my $i (1..$P->{CigLen}[$c])
				{
					if ($Line == 1)
					{
						$FormattedRead1[$RefPos] = ".";
						$FormattedScores1[$RefPos] = " ";
						print "$RefPos\t$FormattedRead1[$RefPos]\t$FormattedScores1[$RefPos]\n" if $Debug; #exit;
					}
					else
					{
						$FormattedRead2[$RefPos] = ".";
						$FormattedScores2[$RefPos] = " ";
						print "$RefPos\t$FormattedRead2[$RefPos]\t$FormattedScores2[$RefPos]\n" if $Debug; #exit;
					}
					$RefPos++;
				}
			}
			elsif ($P->{CigCom}[$c] eq "P") # padding, move ref counter forward, read counter not changed
			{
				print "Line $Line\n$P->{CigCom}[$c]\t$P->{CigLen}[$c]\n" if $Debug;
				for my $i (1..$P->{CigLen}[$c])
				{
					if ($Line == 1)
					{
						$FormattedRead1[$RefPos] = " ";
						$FormattedScores1[$RefPos] = " ";
						print "$RefPos\t$FormattedRead1[$RefPos]\t$FormattedScores1[$RefPos]\n" if $Debug; #exit;
					}
					else
					{
						$FormattedRead2[$RefPos] = " ";
						$FormattedScores2[$RefPos] = " ";
						print "$RefPos\t$FormattedRead2[$RefPos]\t$FormattedScores2[$RefPos]\n" if $Debug; #exit;
					}
					$RefPos++;
				}
			}
			elsif ($P->{CigCom}[$c] eq "D") # deletion, move ref counter forward, read counter not changed
			{
				print "$P->{CigCom}[$c]\t$P->{CigLen}[$c]\n" if $Debug;
				for my $i (1..$P->{CigLen}[$c])
				{
					if ($Line == 1)
					{
						$FormattedRead1[$RefPos] = "-";
						$FormattedScores1[$RefPos] = " ";
						print "$RefPos\t$FormattedRead1[$RefPos]\t$FormattedScores1[$RefPos]\n" if $Debug; #exit;
					}
					else
					{
						$FormattedRead2[$RefPos] = "-";
						$FormattedScores2[$RefPos] = " ";
						print "$RefPos\t$FormattedRead2[$RefPos]\t$FormattedScores2[$RefPos]\n" if $Debug; #exit;
					}
					$RefPos++;
				}
			}
			else
			{
				print "CIGAR format command $P->{CigCom}[$c] is missing from formatting instructions\n";
				print "Found in CIGAR string $P->{Alignment} of read $r of evidence interval $ID\n";
			}
		}
		$HP->{MaxOffset} = $RefPos if $HP->{MaxOffset} < $RefPos;

		my $FormattedRead = join ("",@FormattedRead1)."\n".join ("",@FormattedRead2);
		my $FormattedScores = join ("",@FormattedScores1)."\n".join ("",@FormattedScores2);
		#$Debug = 1;
		if ($Debug)
		{
			print "$FormattedRead\n";
			print "$FormattedScores\n";
		} #exit;
		#foreach my $Pos (sort {$a <=> $b} keys %FormattedScores1) {print "$FormattedScores1{$Pos}";} print "\n";
		#foreach my $Pos (sort {$a <=> $b} keys %FormattedScores2) {print "$FormattedScores2{$Pos}";} print "\n";
		#exit;

		if ($P->{Allele} == 0) { push @{$HP->{Allele0}},$FormattedRead; push @{$HP->{AlleleQ0}},$FormattedScores;}
		elsif ($P->{Allele} == 1) { push @{$HP->{Allele1}},$FormattedRead; push @{$HP->{AlleleQ1}},$FormattedScores;}
		elsif ($P->{Allele} == 2) { push @{$HP->{Allele2}},$FormattedRead; push @{$HP->{AlleleQ2}},$FormattedScores;}
    }

    @{$HP->{Allele0}} = sort { $b cmp $a } @{$HP->{Allele0}};
    @{$HP->{Allele1}} = sort { $b cmp $a } @{$HP->{Allele1}};
    @{$HP->{Allele2}} = sort { $b cmp $a } @{$HP->{Allele2}} if @{$HP->{Allele2}};
    # may need to change this, to sort both on the seq array
    @{$HP->{AlleleQ0}} = sort { $b cmp $a } @{$HP->{AlleleQ0}};
    @{$HP->{AlleleQ1}} = sort { $b cmp $a } @{$HP->{AlleleQ1}};
    @{$HP->{AlleleQ2}} = sort { $b cmp $a } @{$HP->{AlleleQ2}} if @{$HP->{AlleleQ2}};

    return; # aligned reads and scores in records
}

sub GetRuler
{
    my $HP = shift; # pointer to ID rec ref
    my $Chr = shift; # pointer to chr
	#print "Begin: $HP->{Begin}\tEnd: $HP->{End}\n";
	print "Min: $HP->{MinOffset}\tMax: $HP->{MaxOffset}\n";
    #my $Stop = $HP->{MaxOffset};
    my $Diff = $HP->{MaxOffset} - $HP->{MinOffset};
	# Set up ruler records for each allele
	for my $a (0..2) {$HP->{Ruler}[$a] = NewRulerRecord ();} # recs to store ruler

	my $ScaleStart = 10 * int($HP->{MinOffset}/10);
	print "ScaleStart: $ScaleStart\n";
	my @Ref;
    @Ref = split '', `cgatools decodecrr --reference $RefFile --range $Chr:$HP->{MinOffset}-$HP->{MaxOffset}`;
	# Set up numbering, using refpos data
	for my $a (0..2) #
	{
		next unless ($HP->{Alleles} =~ /$a/ or $a == 0); # don't process if no data for this allele
		print "n: $a\talleles: $HP->{Alleles}\n";
		my $P1 = $HP->{Ruler}[$a]; # ptr to Positions array
		my @Legend;
		my $Legend = "";
		my $Scale = "";
		my $RefSeq = "";
		#my $P1 = $HP->{Ruler}[$a]->{Positions}; # ptr to Positions array
		my $s = 0;
		for my $m (0..int(@{$P1->{Positions}})-1)
		{
			if ($P1->{Positions}[$m] == int $P1->{Positions}[$m]) # integer ie ref
			{
				$RefSeq .= $Ref[$s++]; # add ref seq base
				if ($P1->{Positions}[$m] % 10 == 0) # integer = 0
				{
					$Scale .= "|"; # insert marker
					if ($m >= 2) # room to insert nr length 3
					{
						$Legend[$m] = substr($P1->{Positions}[$m],-1,1); # 1st nr
						$Legend[$m-1] = substr($P1->{Positions}[$m],-2,1); # 2nd nr
						$Legend[$m-2] = substr($P1->{Positions}[$m],-3,1); # 3rd nr;
						print "$Legend[$m]\t",$Legend[$m-1],"\t",$Legend[$m-2],"\n"; exit;
					}
					else # no room to insert nr length 3
					{
						$Legend[$m] = " "; # nr spacer
					}
				}
				else # integer ­ 0
				{
					$Scale .= "."; # line spacer
					$Legend[$m] = " "; # nr spacer
				}
			}
			else # non-integer ie insertion
			{
				$RefSeq .= " "; # add ref seq base
				$Legend[$m] = " "; # nr spacer
				$Scale .= "_"; # not ref, so empty spacer
			}
		}
		foreach my $Val (@Legend) {$HP->{Ruler}[$a]->{Legend} .= $Val;}
		print "$HP->{Ruler}[$a]->{Legend}\n";
		print "$Scale\n";
		print "$RefSeq\n\n";
		$HP->{Ruler}[$a]->{Scale} = $Scale;
		$HP->{Ruler}[$a]->{RefSeq} = $RefSeq;
		print $HP->{Ruler}[$a]->{Legend},"\n"; # ptr to scale nrs array
		print $HP->{Ruler}[$a]->{Scale},"\n"; # ptr to scale nrs array
		print $HP->{Ruler}[$a]->{RefSeq},"\n"; # ptr to scale nrs array
	}
	# exit;
}

sub GetEvReadBaseCounts
{
    my $HRef = shift; # reads array ref
    my $ID = shift;
    my $HP = $$HRef{$ID}; # pointer to intvl
    my $Tot = int (@{$$HRef{$ID}->{Reads}});
    # Loop through entries to extract bases by position and count them
    for my $n (0..$Tot-1)
    {
	my $P = $$HRef{$ID}->{Reads}[$n]; # pointer to read array (intvl rec -> read [])
	my $Tot2 = int (@{$P->{CigCom}}); # nr CIGAR commands
	my $ReadPos = 0; # counter for pos in read
	my $RefPos = $P->{OffRef}; # counter for pos in ref
	my $Count = "Count".$P->{Allele};
	my $BaseStr = "Base".$P->{Allele};
	for my $m (0..$Tot2-1) # loop through CIGAR commands
	{
	    if ($P->{CigCom}[$m] eq "M") # match: count, move both counters forward
	    {
		for my $p (1..$P->{CigLen}[$m]) #
		{
		    $HP->{$Count}{$RefPos}++; # increment count for approp allele
		    ${$HP->{$BaseStr}}{$RefPos} .= substr ($P->{Seq},$ReadPos, 1);
		    $RefPos++; # must do inside loop
		    $ReadPos++; # must do inside loop
		}
	    }
	    elsif ($P->{CigCom}[$m] eq "B") # overlap, move ref counter back, read counter not changed
	    {
			$RefPos -= $P->{CigLen}[$m];
	    }
	    elsif ($P->{CigCom}[$m] eq "N" or $P->{CigCom}[$m] eq "P") # gap, move ref counter forward, read counter not changed
	    {
			$RefPos += $P->{CigLen}[$m];
	    }
	    elsif ($P->{CigCom}[$m] eq "I") # insertion, move read counter forward, ref counter not changed
	    {
			$ReadPos += $P->{CigLen}[$m];
	    }
	    elsif ($P->{CigCom}[$m] eq "D") # deletion, move ref counter forward, read counter not changed
	    {
			$RefPos += $P->{CigLen}[$m];
	    }
	    else
	    {
			print "CIGAR format command $P->{CigCom}[$m] is missing from formatting instructions\n";
			print "formatting $P->{Alignment} was being applied";
	    }
	}
    }
    return;
}

sub GetEvReadBases
{
    my $HRef = shift; # intvl hash ref
    for my $Allele (0..2)
    {
	my $BaseStr = "Base".$Allele;
	my $BaseCount = "BasesCounts".$Allele;
	my $P = $HRef->{$BaseStr}; # pointer to ref to hash for storing bases - hash ref eg int hash -> base
	my $P2 = $HRef->{$BaseCount}; # pointer to ref to hash for storing base count strings - hash ref eg int hash -> base
	foreach my $Pos (sort {$a <=> $b} keys %$P) # sorting posns ascending
	{
	    my %Hash; $Hash{"A"} = $Hash{"C"} = $Hash{"G"} = $Hash{"T"} = $Hash{"N"} = 0; # store hash
	    # Count different bases for each position
	    for $n (0..length($$P{$Pos})-1 ) # looping through the bases found at this pos
	    {
		$Hash{substr ($$P{$Pos},$n,1)}++; # add count for each base
	    }
	    # Generate string of base and count for each base present
	    $$P2{$Pos} = "";
	    foreach my $Base (sort {if ($a eq 'N' and $Hash{$a} == $Hash{$b}) { return 1; }
				    elsif ($b eq 'N' and $Hash{$a} == $Hash{$b}) { return -1; }
				    else { return $Hash{$b} <=> $Hash{$a}; }} keys %Hash) # sorting keys posns, largest first, but N always last when equal
	    {
		if ($Hash{$Base})
		{
		    $$P2{$Pos} .= "$Base:$Hash{$Base} ";
		}
	    }
	    $HRef->{MaxOffset} = $Pos if $Pos > $HRef->{MaxOffset}; # waiting to cal till now, to use full range of posns
	}
    }
}

sub PrintEvAlignmentHeader
{
    my $FH = shift;
    print $FH "\nAlignment of reads covering Evidence Interval(s)\n";
}

sub PrintEvAlignments
{
    my $FH = shift; # file handler
    my $ID = shift; # intvl ID
    my $HP = shift; # pointer to ID rec ref
	my $Allele2Exists = int @{$HP->{Allele2}};

	# Print header
    print $FH "\nInterval $ID (interval $HP->{Begin} to $HP->{End}), aligned region from $HP->{MinOffset} to $HP->{MaxOffset})\n";
	# Print allele 0 base alignments
    print $FH "Allele 0 reads\n";
    #print $FH "$Ruler";
	print $FH $HP->{Ruler}[0]->{Legend},"\n";
	print $FH $HP->{Ruler}[0]->{Scale},"\n";
	print $FH $HP->{Ruler}[0]->{RefSeq},"\n";
	print $FH (join ("\n",@{$HP->{Allele0}}),"\n\n");
	# Print allele 1 base alignments
    print $FH "Allele 1 reads\n";
	print $FH $HP->{Ruler}[1]->{Legend},"\n";
	print $FH $HP->{Ruler}[1]->{Scale},"\n";
	print $FH $HP->{Ruler}[1]->{RefSeq},"\n";
    print $FH (join ("\n",@{$HP->{Allele1}}),"\n\n");
	# Print allele 2 base alignments, if exists
    if ($Allele2Exists)
    {
		print $FH "Allele 2 reads\n";
		print $FH $HP->{Ruler}[2]->{Legend},"\n";
		print $FH $HP->{Ruler}[2]->{Scale},"\n";
		print $FH $HP->{Ruler}[2]->{RefSeq},"\n";
		print $FH (join ("\n",@{$HP->{Allele2}}),"\n\n");
    }

	# Print allele 0 qscore alignments
    print $FH "Allele 0 qscores\n";
	print $FH $HP->{Ruler}[0]->{Legend},"\n";
	print $FH $HP->{Ruler}[0]->{Scale},"\n";
	print $FH $HP->{Ruler}[0]->{RefSeq},"\n";
    print $FH (join (@{$HP->{AlleleQ0}}),"\n");
	# Print allele 1 qscore alignments
    print $FH "Allele 1 qscores\n";
	print $FH $HP->{Ruler}[1]->{Legend},"\n";
	print $FH $HP->{Ruler}[1]->{Scale},"\n";
	print $FH $HP->{Ruler}[1]->{RefSeq},"\n";
    print $FH (join (@{$HP->{AlleleQ1}}),"\n\n");
	# Print allele 2 qscore alignments, if exists
    if ($Allele2Exists)
	{
		print $FH "Allele 2 qscores\n";
		print $FH $HP->{Ruler}[2]->{Legend},"\n";
		print $FH $HP->{Ruler}[2]->{Scale},"\n";
		print $FH $HP->{Ruler}[2]->{RefSeq},"\n";
		print $FH (join (@{$HP->{AlleleQ2}}),"\n\n");
	}
}

sub PrintTabbedEvAlignments
{
    my $FH = shift; # file handler
    my $ID = shift; # intvl ID
    my $HP = shift; # pointer to ID rec ref
	my $Allele2Exists = int @{$HP->{Allele2}};

    my ($Count0 , $Count1, $Count2, $Tabbed) = "";

    foreach my $Pos ($HP->{MinOffset}..$HP->{MaxOffset}) # sorting keys posns
    {
		$Count0 .= ${$HP->{Count0}}{$Pos}."\t";
		$Count1 .= ${$HP->{Count1}}{$Pos}."\t";
		$Count2 .= ${$HP->{Count2}}{$Pos}."\t" if int keys %{$HP->{Count2}};
    }

	for my $a (0..2)
	{
		$HP->{Ruler}[$a]->{Scale} = join ("\t", split //, $HP->{Ruler}[$a]->{Scale}); # tabbing ruler
		$HP->{Ruler}[$a]->{Scale} =~ s /$\t//; # correct for extra starting tab
		$HP->{Ruler}[$a]->{RefSeq} = join ("\t", split //, $HP->{Ruler}[$a]->{RefSeq}); # tabbing ruler
		$HP->{Ruler}[$a]->{RefSeq} =~ s /$\t//; # correct for extra starting tab
	}

    # Allele 0 base count
    $Tabbed = join ("\n", @{$HP->{Allele0}});
    $Tabbed = join ("\t", split //, $Tabbed); # print matched seq
    $Tabbed =~ s /\n\t( ){0,1}/\n/g; # correct for extra tab after return

    print $FH "\nInterval $ID (interval $HP->{Begin} to $HP->{End}), aligned region from $HP->{MinOffset} to $HP->{MaxOffset})\n";

    # Print allele 0
    print $FH "\nAllele 0 reads\n";
	print $FH join("\t",@{$HP->{Ruler}[0]->{Positions}}),"\n";
	print $FH $HP->{Ruler}[0]->{Scale},"\n";
	print $FH $HP->{Ruler}[0]->{RefSeq},"\n";
    print $FH "$Count0\n";
    print $FH "$Tabbed\n";

    # Print allele 1
    $Tabbed = join ("\n", @{$HP->{Allele1}});
    $Tabbed = join ("\t", split //, $Tabbed); # print matched seq
    $Tabbed =~ s /\n\t( ){0,1}/\n/g; # correct for extra tab after return

    print $FH "\nAllele 1 reads\n";
	print $FH join("\t",@{$HP->{Ruler}[1]->{Positions}}),"\n";
	print $FH $HP->{Ruler}[1]->{Scale},"\n";
	print $FH $HP->{Ruler}[1]->{RefSeq},"\n";
    print $FH "$Count1\n";
    print $FH "$Tabbed\n";

    # Print allele 2
    if ($Allele2Exists)
    {
		$Tabbed = join ("\n", @{$HP->{Allele2}});
		$Tabbed = join ("\t", split //, $Tabbed); # print matched seq
		$Tabbed =~ s /\n\t( ){0,1}/\n/g; # correct for extra tab after return

		print $FH "\nAllele 2 reads\n";
		print $FH join("\t",@{$HP->{Ruler}[1]->{Positions}}),"\n";
		print $FH $HP->{Ruler}[2]->{Scale},"\n";
		print $FH $HP->{Ruler}[2]->{RefSeq},"\n";
		print $FH "$Count2\n\n";
		print $FH "$Tabbed\n";
	}

	# Allele 0 qscore count
    $Tabbed = join ("\n", @{$HP->{AlleleQ0}});
    $Tabbed = join ("\t", split //, $Tabbed); # print matched seq
    $Tabbed =~ s /\n\t( ){0,1}/\n/g; # correct for extra tab after return

    # Print allele 0
    print $FH "\nAllele 0 qscores\n";
	print $FH join("\t",@{$HP->{Ruler}[0]->{Positions}}),"\n";
	print $FH $HP->{Ruler}[0]->{Scale},"\n";
	print $FH $HP->{Ruler}[0]->{RefSeq},"\n";
    print $FH "$Count0\n\n";
    print $FH "$Tabbed\n";

	# Allele 1 qscore count
    $Tabbed = join ("\n", @{$HP->{AlleleQ1}});
    $Tabbed = join ("\t", split //, $Tabbed); # print matched seq
    $Tabbed =~ s /\n\t( ){0,1}/\n/g; # correct for extra tab after return

    # Print allele 1
    print $FH "\nAllele 1 qscores\n";
	print $FH join("\t",@{$HP->{Ruler}[1]->{Positions}}),"\n";
	print $FH $HP->{Ruler}[1]->{Scale},"\n";
	print $FH $HP->{Ruler}[1]->{RefSeq},"\n";
    print $FH "$Count1\n\n";
    print $FH "$Tabbed\n";

	if ($Allele2Exists)
	{
		# Allele 2 qscore count
		$Tabbed = join ("\n", @{$HP->{AlleleQ2}});
		$Tabbed = join ("\t", split //, $Tabbed); # print matched seq
		$Tabbed =~ s /\n\t( ){0,1}/\n/g; # correct for extra tab after return

		# Print allele 0
		print $FH "\nAllele 2 qscores\n";
		print $FH join("\t",@{$HP->{Ruler}[1]->{Positions}}),"\n";
		print $FH $HP->{Ruler}[2]->{Scale},"\n";
		print $FH $HP->{Ruler}[2]->{RefSeq},"\n";
		print $FH "$Count2\n\n";
		print $FH "$Tabbed\n";
	}

}

sub PrintEvReadBaseCountHeader
{
    my $FH = shift;
    print $FH "\nCounts of bases per position across Evidence Interval(s)\n";
}

sub PrintEvReadBaseCounts
{
    my $FH = shift;
    my $ID = shift;
    my $HP = shift; # reads array ref
    print $FH "\nInterval $ID\n";
    print $FH "Pos\tCount 0\tCount 1\tCount 2\tBases 0\tBases 1\tBases 2\n";
    for my $Pos ($HP->{MinOffset}..$HP->{MaxOffset})
    {
		print $FH $Pos."\t${$HP->{Count0}}{$Pos}\t${$HP->{Count1}}{$Pos}\t${$HP->{Count2}}{$Pos}\t";
		print $FH ${$HP->{BasesCounts0}}{$Pos}."\t".${$HP->{BasesCounts1}}{$Pos}."\t".${$HP->{BasesCounts2}}{$Pos}."\n";
    }
}

sub PrintVarHeaders
{
    print $OUT1 "Var Records in file $VarFile\n";
    my $H = shift;
    print $OUT1 $H;
}

sub PrintHeaders1
{
    my $File = shift;
    my $Term = shift;
    my $ColH = shift;
    print $OUT1 "\n$Term Record(s) from file $File\n";
    print $OUT1 $ColH;
}

sub PrintData1
{
    my $ARef = shift;
    print $OUT1 join ("\n", @$ARef)."\n";
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
	return my @Array if $Count++ > 50; # too many lines for a header, something is wrong, return empty array
    }
}

sub GetVersion
{
    my $ARef = shift;
	my $Key;
	my $Val;
    foreach my $Line (@$ARef)
    {
		($Key, $Val) = split /\t/, $Line;
		return $Val if ($Key =~ /#SOFTWARE_VERSION/);
    }
    print "Version $Key not recognised\n";
	print "Data generated by this Software Version may not be analysed correctly\n";
	return 0;
}

sub GetColHeader
{
    my $FH = shift;
    my $Count = 0;
    while (<$FH>)
    {
	return $_ if (/^>/); # loop until a line starts with >
	return "" if $Count++ > 3; # too many lines for a header, something is wrong, return empty array
    }
}

sub FindCoverageRefScoreFile
{
    my $Dir = shift;
    my $Chr = shift;
    opendir my $DH, $Dir or die "Cannot open dir $Dir";
	#print "$Dir\n";
    my @Files = readdir $DH;
    close $DH;
    foreach my $File (@Files)
    {
		#print "$File\n";
		if ($File =~ /^coverageRefScore-$Chr-/) {return $File;}
    }
    return 0;
}

sub Max2Ints
{
    my $One = shift;
    my $Two = shift;
    if ($One > $Two) {return $One;}
    else 	     {return $Two;}
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
    elsif ($File =~ /.tsv$/)
    {
		open ($FH, "cat $File |") or die ("$!: can't open file $File");
    }
    else
    {
		die ("$!: do not recognise file type $File");
    }
    return $FH;
}


sub NewIntervalRec
{
    my $Record =
    {
	IntEntry => "", # full evint entry
	Chr => "", # chr
	Alleles => "", # alleles, 2 sep by ;
	Begin => 10**10, # first pos in evint
	End => -1, # last pos in evint
	IntAlleleSeqs => [], # allele sequences
	IntAlleleCIGARs => [], # allele CIGARs
	DNBs => [], # reads from evDNB file
	DNBCounts => [], #
	Reads => [], # array to store struct of read data extracted from DNBs: both mate-pairs and qual scores etc
	#ScoreStrings => [],
	MinOffset => 10**10, # first pos of left most read
	MaxOffset => -1, # last pos of right most read
	AlleleRuler => [], # struct to store a ruler for each allele, to display evint alleles
	Ruler => [], # struct to store a ruler for each allele, to display reads
	Allele0 => [],
	Allele1 => [],
	Allele2 => [],
	AlleleQ0 => [],
	AlleleQ1 => [],
	AlleleQ2 => [],
	Count0 => {},
	Count1 => {},
	Count2 => {},
	Bases0 => {},
	Bases1 => {},
	Bases2 => {},
	BasesCounts0 => {},
	BasesCounts1 => {},
	BasesCounts2 => {},
    };
    return $Record;
}

sub NewReadRecord
{
    my $Record =
    {
        ID =>,
        Allele =>,
		Arm => "",
        Strand => "",
        Offset =>,
        OffRef =>,
        Alignment => "",
        Seq => "",
		Qual => "",
		FSeq => "",
		FQual => "",
		CigCom => [],
		CigLen => [],

    };
    return $Record;
}

sub NewRulerRecord
{
    my $Record =
    {
        CIGAR => 0,
        Positions => [],
        Seq => [],
		Legend => "",
        Scale => "",
        RefSeq => "",
        AlleleSeq => "",

    };
    return $Record;
}

sub NewCIGARRecord
{
    my $Record =
    {
        CigCom => [],
        CigLen => [],
    };
    return $Record;
}


############################

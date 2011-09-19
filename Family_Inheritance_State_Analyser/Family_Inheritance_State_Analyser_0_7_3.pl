#!/usr/bin/perl
use strict;
#use feature "say";
use File::Basename;
$| = 1;

# Family_Inheritqnce_State_Analyser
# From parental genotype data, detemine which genotypes can be assigned unambiguously (usually aa/ab)
# [Not analysing haploids at the moment]
# For these aa/ab parental genotypes, look at child inheritance patterns
# From aa/ab, children can be aa or ab ie can be the same or different to each other
# This pattern should be stable in blocks, which will reflect recombination events
# This block information can be used to look for variants which are the same or different in the children

# Looking at parental genotyes:
# aa/bb and aa/aa are uninformative
# following are MIEs
# aa/aa giving ab or bb in a child is MIE unexpected allele
# aa/bb for aa or bb, in a child as MIE impossible combination
# aa/ab for bb, in a child as MIE impossible combination
# aa/ab potentially informative for state inheritance

# input is a testvariants results file

# Flow
# 1. Take output of listvariants/testvariants
# 2. For all loci, check whether the childrens' genotypes can arise from the parents' genotypes
# 3. Remove those loci which are discordant ie MIEs
# 4. Mark loci where one parent is aa (or a) and other is ab (or b) ie potentially informative
# 5. For these loci, are the children are identical (S) or not (D) for the informative allele
# 4. Mark loci where one both parents are ab ie potentially informative
# 6. For these loci are the children identical (S) for the informative allele
# 7. Save informative aa/ab and ab/ab loci
# 8. Concatenate blocks of loci with the same state
# 9. Remove singleton loci and concatenate state blocks again
# 10. Mark blocks of short runs (user defined) as having an uncertain state (U)
# 11. Mark gaps between states (G)

# To Do: extend to start and end of chromosomes (ie begin and end of called chrs)

# Format:
# perl prog file dir nr nr nr [nr]
# ie
# perl Family_Inheritance_State_Analyser \
# --Input input_file \
# --Output_Dir output_dir \
# --Mat_Field mat_field_nr \
# --Pat_Field pat_field_nr \
# --Child_Fields child_field_nr [child_field_nr] \
# --Short_Block nr
# eg
# perl /perl/Family_Inheritance_State_Analyser_0_7_1.pl \
# --Input /data/Family_Quartet_testvariants.tsv \
# --Output_Dir /data/ \
# --Mat_Field 8  --Pat_Field 9  --Child_Fields 10 11  --Short_Block 5 \

# Rick Tearle 2010-11

# testvariants fields:
# variantId chromosome begin end varType reference alleleSeq xRef GS000000XX1-ASM GS000000XX2-ASM [GS000000XXN-ASM] .... sum0 sum1 [sumN]
my $Time -= time; # start time
my $ScriptName = basename ($0); # get script name
my $Length2Ignore = 1000; # not used at the moment
my $MIEBlockWindow = 1_000_000; # for binning MIEs
my $GapNot2Cross = 1_000_000;
my $Debug = 0;
my %Colours = SetColours ();

# Parsing and storing input parameters
# Only childfields can be repeated
print "\nProcessing input parameters\n";
my %ExpectedParams =  GetExpectedParams ();
my %EnteredParams = GetEnteredParams ();
TestParameters (%EnteredParams); # no return val used, dies in subroutine if no correct

# Setting up prog paras from input paras
my $FileIn = $EnteredParams{input};
unless (-f $FileIn) {die "Testvariants input file $FileIn not found\n";} # requires existing file
my $DirectoryOut = $EnteredParams{output_dir}; #
$DirectoryOut =~ s/\/$//; # remove trailing slash if present
unless (-d $DirectoryOut) {die "Output directory $DirectoryOut not found\n";} # requires existing file
my $MatField = int $EnteredParams{mat_field};
unless ($MatField) {die "Need nr of field containing maternal genome data\n";}
my $PatField = int $EnteredParams{pat_field};
unless ($PatField) {die "Need nr of field containing paternal genome data\n";}
my @ChildFields = split / /, $EnteredParams{child_fields};
my $NrChildren = int @ChildFields;
if ($NrChildren == 0) {die "Need nr of field containing at least one child's genome data\n";}
elsif ($NrChildren == 1) {print "Only 1 field containing a child's genome data provided\nWill only calculate MIEs\n";}
my $ShortBlock = $EnteredParams{short_block} || $ExpectedParams{short_block};
if ($ShortBlock != int $ShortBlock) {die "Length of block to resolve must be an integer, not $ShortBlock\n";}
elsif ($ShortBlock == 0) 			{die "Length of block to resolve must be at least 1, currently set at $ShortBlock\n";}
elsif ($ShortBlock == 1) 			{print "Blocks of Length 1 are automatically resolved, will not resolve further\n";}
#print "$FileIn\n$DirectoryOut\n$MatField\n$PatField\n",join(" ",@ChildFields),"\n$ShortBlock\n";

my ($FileOut, undef, $Ext) = fileparse ($FileIn); #, qr/\.[^.]*/);
$FileOut =~ s/(?:\.[^.]+)+$//;
#print "$FileOut\n"; exit;
my $n = 1;
my $Suff = "";
while (-f $DirectoryOut."/".$FileOut.$Suff."_State.tsv" and -f $DirectoryOut."/".$FileOut.$Suff."_MIE.tsv" )
{$Suff = "_$n"; $n++;} # loop t	ill we have a new unique filename

print "Assigning Fields:\n";
print "\tMaternal Field\t$MatField\n\tPaternal Field\t$PatField\n\tChild Fields\t",join("\t",@ChildFields),"\n";

print "\nOpening input file $FileIn\n";
my $IN = OpenFile ($FileIn); # open the file with correct file format
my $ColHeader = <$IN>; # get col header
chomp $ColHeader;
my @ColHeader = split /\t/, $ColHeader;

print "Creating subdirs in dir $DirectoryOut\n";
my $DirOutBed = $DirectoryOut."/bed";
unless (-d $DirOutBed) {mkdir $DirOutBed or die "Cannot create output directory $DirOutBed";}
my $DirOutStates = $DirectoryOut."/States";
unless (-d $DirOutStates) {mkdir $DirOutStates or die "Cannot create output directory $DirOutStates";}
my $DirOutStatesGT1 = $DirectoryOut."/States_GT_1";
unless (-d $DirOutStatesGT1) {mkdir $DirOutStatesGT1 or die "Cannot create output directory $DirOutStatesGT1";}
my $DirOutStatesGT5;
if ($ShortBlock > 1)
{
	$DirOutStatesGT5 = $DirectoryOut."/States_GT_$ShortBlock";
	unless (-d $DirOutStatesGT5) {mkdir $DirOutStatesGT5 or die "Cannot create output directory $DirOutStatesGT5";}
}
print "Opening output files:\n";
my $NrStateFields;
my $StateFileOut = $FileOut.$Suff."_State.tsv";
my $StateColHeader = "";
open my $OUT1, ">", $DirectoryOut."/".$StateFileOut;
print "\t",$StateFileOut,"\n"; #exit;
print $OUT1 join("\t",@ColHeader),"\t";
#print $OUT1 join("\t",@ColHeader[0..3]),"\t",join("\t",@ColHeader[8..int(@ColHeader)-1]),"\t";

for my $n (1..$NrChildren-1) # Generating col header for state file
{
    for my $m ($n+1..$NrChildren)
    {
	$StateColHeader .= "Ch$n"."_Ch$m"."_Mat\t"."Ch$n"."_Ch$m"."_Pat\t";
	$NrStateFields += 2;
    }
}
print $OUT1 "$StateColHeader\n"; # col header to file

my $FileOut2 = $FileOut.$Suff."_MIE.tsv";
my $MIEColHeader = "";
open my $OUT2, ">", $DirectoryOut."/".$FileOut2;
print "\t",$FileOut2,"\n"; #exit;
print $OUT2 "$ColHeader\t";
for my $n (1..$NrChildren) {$MIEColHeader .= "MIE_Ch$n\t";} # col header to file
print $OUT2 "$MIEColHeader\n"; # col header to file

my $FileOut3 = $FileOut.$Suff."_MIE.bed";
my $MIEBedColHeader = "track name=MIEs description=MIEs itemRgb=on";
open my $OUT3, ">", $DirOutBed."/".$FileOut3;
print "\t",$FileOut3,"\n"; #exit;
print $OUT3 "$MIEBedColHeader\n";

my $FileOut4 = $FileOut.$Suff."_MIE_Blocks.bed";
my $MIEBlockBedColHeader = "track name=MIE_Blocks description=MIE_Blocks itemRgb=on";
open my $OUT4, ">", $DirOutBed."/".$FileOut4;
print "\t",$FileOut4,"\n"; #exit;
print $OUT4 "$MIEBlockBedColHeader\n";

# Loading chr nrs - note that chrX, chrY and chrM are not loaded
my @ChrNames = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
	   'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22');
my %Vars; # hash to store var records in an array for each chr
my %MIE_Blocks; # hash to store MIE block records in hash for each chr

foreach my $Chr (@ChrNames)
{
    $Vars{$Chr} = [];  # set up array for each chr, as hash
    $MIE_Blocks{$Chr} = {};  # set up hash for each chr, as hash
}

######## Process Test Variants ########
# Takes first and second nominated columns in Args as maternal, paternal
# Takes any other nominated columns in Args as children
# new columns for MIE file
# col 1 is for MIEs for child1 (1 = MIE)
# ..
# col n is for MIEs for childn (1 = MIE)

# new columins for state file
# col 1 states whether child1 and child2 are same (S) or diff (D) for an informative maternally derived allele
# col 2 states whether child1 and child2 are same (S) or diff (D) for indormative paternally derived allele
# ..
# col n!-1 states whether childn-1 and childn are same (S) or diff (D) for an informative maternally derived allele
# col n! states whether childn-1 and childn are same (S) or diff (D) for an informative paternally derived allele

my %States; # hash to store state block records in an array for each chr
foreach my $Chr (@ChrNames)
{
    next if $Chr eq "chrX" or $Chr eq "chrY" or $Chr eq "chrM"; # not processing these yet
    $States{$Chr} = [];  # set up ptr to array for each chr, as hash
}

print "\nProcessing MIEs and States....\n";
my $LastFieldIn = (() = $ColHeader =~ /\t/g); # get count of fields in testvariants output
my $FirstNewField = $LastFieldIn + 1;
my @Variants;
my $MIECount = 0; # MIE counter
my $VariantCount = 0; # variant locus counter
my $PotInfCount = 0; # potentially state informative counter ie right parental genotypes
my $InfCount = 0; # state informative counter
my @AnalyseStateSubs = ("",\&AnalyseState1, \&AnalyseState2); # array storing refs to 2 analysestate subs - not using yet?

while (<$IN>)
{
    chomp;
    my @Fields = split /\t/, $_;
#last if $Fields[1] eq "chr3";
    next if $Fields[1] eq "chrX" or $Fields[1] eq "chrY" or $Fields[1] eq "chrM";
    # not processing X yet, Y and M show no recombination and should not be processed for mapping, but should be for MIEs, not processing yet
    $VariantCount++; # increment variant locus counter - only those that are processed

    # Find MIEs #
    # type 1: aa aa-> b
    # type 2: aa bb -> aa
    # type 3: aa ab -> bb

	my ($MIE_Flag, $Informative, $Fields) = CheckForMIEs (\@Fields); # pass testvariant fields, return MIE, Informative, status, updated fields

	if ($MIE_Flag) # is an MIE, save to file but do not use for state mapping
	{
		$MIECount++; # increment MIE counter
		$MIE_Blocks{$Fields[1]}{int($Fields[2]/$MIEBlockWindow)}++; # increment counter for MIE Block
		SaveRec ($OUT2, $Fields); # save rec to MIE file
		SaveRecAsMIEBed ($OUT3, $Fields); # save bed rec to file
	}
	elsif ($Informative) # not an MIE, but state informative so process
	{
		$PotInfCount++; # increment potentially state informative  counter
		my $Use = 0;
		#($Use, $Fields) = &{$AnalyseStateSubs[$Informative]};

		if ($Informative == 1) # parents aa ab from analysis above, check kids
		{
			($Use, $Fields) = AnalyseState1 ($Fields); # pass fields, return whether to use ie state informative, and updated fields
		}
		elsif ($Informative == 2) # parents ab ab from analysis above, check kids
		{
			($Use, $Fields) = AnalyseState2 ($Fields); # pass fields, return whether to use ie state informative, and updated fields
		}

		if ($Use) # record is informative for state
		{
			$InfCount++; # increment state informative counter
			SaveRec ($OUT1, $Fields); # save rec to state file
			delete @$Fields[4..7]; # remove unneeded fields, to reduce size of kept record array
			push @{$Vars{$Fields[1]}}, $Fields; # adding ptr to this record to array for this chr
		}
	}

}
close $OUT1; close $OUT2; close $OUT3;

print "\nNr MIE Containing Loci\t$MIECount\n";
print "Nr Potentially Informative Loci\t$PotInfCount\n";
print "Nr Actual Informative Loci\t$InfCount\n";
print "Nr Variant Loci\t$VariantCount\n"; #exit;

# Generate counts for MIE bins
print "Generating stats for MIEs bins....\n";
# Get mean and SD
my ($N, $Mean, $SD, $LTail, $RTail) = GetMIEBinStats (\%MIE_Blocks);
print "Bin Stats:\n\tN\tMean\tSD\n";
print "\t$N\t$Mean\t$SD\n";
print "Saving MIEs bins as bed file....\n";
SaveMIEBins (\%MIE_Blocks, $OUT4, $LTail, $RTail);
close $OUT4;

######### Concatenate States into Runs ########
print "\nGenerating and Concatenating State Files for Child Combinations\n\t";
my $StateBlockHeader = "Chr\tBegin\tEnd\tState\tLength\tRecords";
my @StateColHeader = split /\t/,$StateColHeader;
print join ("\n\t",@StateColHeader),"\n";
for my $n (0..$NrStateFields-1) # processing each child combo
{
    # Generate Blocks with the same state
	#print "\t$StateColHeader[$n]\n"; exit;
    my $FileOut1 = $FileOut.$Suff."_".$StateColHeader[$n]."_States.tsv";
    open $OUT1, ">", "$DirOutStates/$FileOut1";
    print $OUT1 $StateBlockHeader,"\n"; #exit;
    my $FileOut2 = $FileOut.$Suff."_".$StateColHeader[$n]."_States.bed";
    open $OUT2, ">", "$DirOutBed/$FileOut2";
    my $BedHeader = "track name=$StateColHeader[$n] description=States itemRgb=on";
    print $OUT2 "$BedHeader\n";
    foreach my $Chr (@ChrNames)
    {
		next unless int(@{$Vars{$Chr}}); # don't process if no entries
		#print "here 1a\t$Chr\t$States{$Chr}\t$Vars{$Chr}\n";
		$States{$Chr} = ConcatenateVariants ($Vars{$Chr}, $FirstNewField+$n); # concatenating child states
		#print "here 1b\t$Chr\n";
		SaveStateBlocks ($OUT1, $States{$Chr});
		SaveStateBlocksAsBed ($OUT2, $States{$Chr});
    }

    # Resolving Blocks of 1
    $FileOut1 = $FileOut.$Suff."_".$StateColHeader[$n]."_States_GT_1.tsv";
    open $OUT1, ">", "$DirOutStatesGT1/$FileOut1";
    print $OUT1 "$StateBlockHeader\n";
    $FileOut2 = $FileOut.$Suff."_".$StateColHeader[$n]."_States_GT_1.bed";
    open $OUT2, ">", "$DirOutBed/$FileOut2";
    $BedHeader = "track name=$StateColHeader[$n]"."_GT_1 description=States itemRgb=on";
    print $OUT2 "$BedHeader\n";
    foreach my $Chr (@ChrNames)
    {
		next unless int(@{$States{$Chr}}); # don't process if no entries
		my $Count = 1; # set to true for testing
		my $Loops = 0;
		while ($Count) # keep processing till no more singletons
		{
			($Count, $States{$Chr}) = ResolveSingletons ($States{$Chr});
			last if $Loops++ > 20; # catch all for infinite loops, should not occur
		}
		SaveStateBlocks ($OUT1, $States{$Chr});
		SaveStateBlocksAsBed ($OUT2, $States{$Chr});
    }

    # Converting blocks below length cutoff to U
	next unless $ShortBlock > 1; # already doing singletons, only proceed if > 1
    #print "here 3\n";
    $FileOut1 = $FileOut.$Suff."_".$StateColHeader[$n]."_States_GT_".$ShortBlock.".tsv";
    open $OUT1, ">", "$DirOutStatesGT5/$FileOut1";
    $FileOut2 = $FileOut.$Suff."_".$StateColHeader[$n]."_States_GT_".$ShortBlock.".bed";
    print $OUT1 "$StateBlockHeader\n";
    open $OUT2, ">", "$DirOutBed/$FileOut2";
    $BedHeader = "track name=$StateColHeader[$n]"."_GT_".$ShortBlock." description=States\ itemRgb=on";
    print $OUT2 "$BedHeader\n";
    foreach my $Chr (@ChrNames)
    {
		my $Count = 1;
		my $Loops = 0;
		next unless int(@{$States{$Chr}}); # don't process if no entries
		while ($Count) # keep processing till no more singletons
		{
			#print "HERE\n";
			($Count, $States{$Chr}) = ResolveShortBlocks ($States{$Chr}, $ShortBlock);
			last if $Loops++ > 20; # catch all for infinite loops, should not occur
		}
		#print "\n";
		$States{$Chr} = FillInGaps ($States{$Chr});
		SaveStateBlocks ($OUT1, $States{$Chr});
		SaveStateBlocksAsBed ($OUT2, $States{$Chr});
    }
}


$Time += time;
print "\ntime $Time\n";

###########################################################################
#                                   SUBS                                  #
###########################################################################

sub GetExpectedParams
{
	my %Hash =
	(
	input => "",
	output_dir => "",
	mat_field => -1,
	pat_field => -1,
	child_fields => "",
	short_block => 5, # optional, default = 5
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
		if ($Key eq "child_fields") # special case, may have repeated params
		{
			$Hash{$Key} .= $Val; # make a hash entry out of key and val, possible multiple entries
		}
		else # all others expect one entry
		{
			$Hash{$Key} = $Val; # make a hash entry out of key and val
		}
		#print "$Key\t$EnteredParams{$Key}\n" if $Debug;
	}
		#print int(keys %Hash),"\n" if $Debug;
		#foreach my $Arg (keys %Hash) {print "Arg: $Arg\t",$ExpectedParams{$Arg},"\n";}
		#print "Arg string:\t",join (" ",@ARGV),"\n" if $Debug;
		#exit if $Debug;
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

sub SetColours
{
	my %Colours;
	$Colours{Red} = "255,0,0"; # red
	$Colours{Blue} = "0,0,255"; # blue
	$Colours{Green} = "0,255,0"; # green
	$Colours{Purple} = "255,0,255"; # purple
	$Colours{PaleBlue} = "150,150,255"; # pale blue
	$Colours{Grey} = "150,150,150"; # grey
	$Colours{PaleGrey} = "50,50,50"; # pale grey

	$Colours{SameClr} = $Colours{Red}; # S = red
	$Colours{DiffClr} = $Colours{Blue}; # D = blue
	$Colours{UncertClr} = $Colours{Grey}; # U = grey
	$Colours{GapClr} = $Colours{Purple}; # G = purple

	return %Colours;
}

sub CheckForMIEs
{
	my $Fields = shift;
	my $MIE_Flag = my $Informative = 0;

    if ($$Fields[$MatField] eq $$Fields[$PatField]) # parents have same genoytpe #
    {
		if ($$Fields[$MatField] eq "00") # parents both 00, kids must be 00
		{
			for my $n (0..$NrChildren-1)
			{
				if ($$Fields[$ChildFields[$n]] =~ /1/) # child n carries at least 1 allele 1
				{
					$$Fields[$FirstNewField+$n] = 1; # MIE in child n
					$MIE_Flag = 1;
				}
				else
				{
					$$Fields[$FirstNewField+$n] = 0; # no MIE in child n
				}
			}
		}
		elsif ($$Fields[$MatField] eq "11") # parents both 11, kids must be 11
		{
			for my $n (0..$NrChildren-1)
			{
				if ($$Fields[$ChildFields[$n]] =~ /0/) # child n carries at least 1 allele 0
				{
					$$Fields[$FirstNewField+$n] = 1; # MIE in child n
					$MIE_Flag = 1;
				}
				else
				{
					$$Fields[$FirstNewField+$n] = 0; # no MIE in child n
				}
			}
		}
		elsif ($$Fields[$MatField] eq "01") # parents both 01, kids can be 00 01 or 11 so not looking at MIEs, looking for informative
		{
			$Informative = 2; # parents are ab/ab so maybe informative for hapmapping, will check children below
		}
    }
    elsif (($$Fields[$MatField] eq "00" and $$Fields[$PatField] eq "11") or
	   ($$Fields[$MatField] eq "11" and $$Fields[$PatField] eq "00")) # one parent 00, other 11 so kids must be 01 #
    {
		for my $n (0..$NrChildren-1)
		{
			if ($$Fields[$ChildFields[$n]] !~ /N/ and $$Fields[$ChildFields[$n]] ne "01") # child n is not heterozygous and not no-called
			{
				$$Fields[$FirstNewField+$n] = 2; # MIE in child n
				$MIE_Flag = 1;
			}
			else
			{
				$$Fields[$FirstNewField+$n] = 0; # no MIE in child n
			}
		}
	}
	elsif (($$Fields[$MatField] eq "00" and $$Fields[$PatField] eq "01") or
		   ($$Fields[$MatField] eq "01" and $$Fields[$PatField] eq "00")) # parents are 00 01, kids cannot be 11
	{
			$Informative = 1; # parents are aa/ab ie informative for hapmapping
			for my $n (0..$NrChildren-1)
			{
				if ($$Fields[$ChildFields[$n]] eq "11") # child n is homozygous for 1
				{
					$$Fields[$FirstNewField+$n] = 3; # MIE in child n
					$MIE_Flag = 1;
				}
				else
				{
					$$Fields[$FirstNewField+$n] = 0; # no MIE in child n
				}
			}
	}
	elsif (($$Fields[$MatField] eq "11" and $$Fields[$PatField] eq "01") or
		   ($$Fields[$MatField] eq "01" and $$Fields[$PatField] eq "11")) # parents are 11 01
	{
		$Informative = 1; # parents are aa/ab ie informative for hapmapping
		for my $n (0..$NrChildren-1)
		{
			if ($$Fields[$ChildFields[$n]] eq "00") # child n is homozygous for 0
			{
				$$Fields[$FirstNewField+$n] = 3; # MIE in child n
				$MIE_Flag = 1;
			}
			else
			{
				$$Fields[$FirstNewField+$n] = 0; # no MIE in child n
			}
		}
    }
	return ($MIE_Flag, $Informative, $Fields);
}

sub SaveRec
{
	my $FH = shift;
	my $Fields = shift;
	print $FH join("\t",@$Fields),"\n";
}

sub SaveRecAsMIEBed
{
    my $FH = shift;
    my $Array = shift;
    my $RGB;
    #print join("\t",@$Array),"\n"; exit;
	my $Name;
	if ($$Array[12] == 1 or $$Array[13] == 1) {$RGB = $Colours{Red}; $Name = 1;}
	elsif ($$Array[12] == 2 or $$Array[13] == 2) {$RGB = $Colours{Blue}; $Name = 2;}
	elsif ($$Array[12] == 3 or $$Array[13] == 3) {$RGB = $Colours{Green}; $Name = 3;}
	else                       {$RGB = $Colours{UncertClr}; $Name = 4;}

	print $FH $$Array[1],"\t",$$Array[2],"\t",$$Array[3],"\t",$$Array[4],"\t",$Name,"\t";
	print $FH "+","\t",$$Array[2],"\t",$$Array[3],"\t",$RGB,"\n";
}

sub AnalyseState1 # parents 00 01 or 00 10 or 01 00 or 01 11 or 11 01 or 11 10
{
	my $Fields = shift; # ptr to modified testvariant fields
	my $Use = 0; # whether to use for state analysis
	my $Combo = 0; # counter for which field to enter data for, increments by 2 to allow 2 fields per combination, one maternal, one paternal
	#print "Nr Fields: ",int(@$Fields),"\nNr Children: ",$NrChildren,"\n";
	#$$Fields[int(@$Fields)-1+($NrChildren*2)] = "";
	#print "Nr Fields: ",int(@$Fields),"\n"; exit;
	if ($$Fields[$MatField] eq "01" and $$Fields[$PatField] eq "00") # allele 1 is unique for mother
	{
		my $Combo = 0; # counter for which field to enter data for, increments by 2 to allow 2 fields per combination, one maternal, one paternal
		for my $Ch1 (0..$NrChildren-2) # loop for first child in pair
		{
			for my $Ch2 ($Ch1+1..$NrChildren-1) # loop for second child in pair
			{
				if (($$Fields[$ChildFields[$Ch1]] =~ /1/ and $$Fields[$ChildFields[$Ch2]] =~ /1/) or
					($$Fields[$ChildFields[$Ch1]] eq "00" and $$Fields[$ChildFields[$Ch2]] eq "00")) # two kids same genotype, autosomes only for now
				{
					$$Fields[$FirstNewField+$Combo] = "S"; # kids have same genoytope from mother
					$Use = 1;
				}
				elsif ($$Fields[$ChildFields[$Ch1]] !~ /N/ and $$Fields[$ChildFields[$Ch2]] !~ /N/) # only want fully called
				{
					$$Fields[$FirstNewField+$Combo] = "D"; # kids have diff genoytope from mother
					$Use = 1;
				}
				else
				{
					$$Fields[$FirstNewField+$Combo] = "N"; # not informative genotype from mother
				}
				$$Fields[$FirstNewField+$Combo+1] = "N"; # not informative genotype from father
				$Combo += 2; # go to next set of fields for next child
			}
		}
	}
	elsif ($$Fields[$MatField] eq "01" and $$Fields[$PatField] eq "11") # allele 0 is unique for mother
	{
		my $Combo = 0;
		for my $Ch1 (0..$NrChildren-2)
		{
			for my $Ch2 ($Ch1+1..$NrChildren-1)
			{
				if (($$Fields[$ChildFields[$Ch1]] =~ /0/ and $$Fields[$ChildFields[$Ch2]] =~ /0/) or
					($$Fields[$ChildFields[$Ch1]] eq "11" and $$Fields[$ChildFields[$Ch2]] eq "11")) # two kids same genotype, autosomes only for now
				{
					$$Fields[$FirstNewField+$Combo] = "S"; # kids have same genoytope from mother
					$Use = 1;
				}
				elsif ($$Fields[$ChildFields[$Ch1]] !~ /N/ and $$Fields[$ChildFields[$Ch2]] !~ /N/) # only want fully called
				{
					$$Fields[$FirstNewField+$Combo] = "D"; # kids have diff genoytope from mother
					$Use = 1;
				}
				else
				{
					$$Fields[$FirstNewField+$Combo] = "N"; # # not informative genotype from mother
				}
				$$Fields[$FirstNewField+$Combo+1] = "N"; # not informative genotype from father
				$Combo += 2;
			}
		}
	}
	elsif ($$Fields[$MatField] eq "00" and $$Fields[$PatField] eq "01") # allele 1 is unique for father
	{
		my $Combo = 0;
		for my $Ch1 (0..$NrChildren-2)
		{
			for my $Ch2 ($Ch1+1..$NrChildren-1)
			{
				if (($$Fields[$ChildFields[$Ch1]] =~ /1/ and $$Fields[$ChildFields[$Ch2]] =~ /1/) or
					($$Fields[$ChildFields[$Ch1]] eq "00" and $$Fields[$ChildFields[$Ch2]] eq "00")) # two kids same genotype, autosomes only for now
				{
					$$Fields[$FirstNewField+$Combo+1] = "S"; # kids have same genoytope from father
					$Use = 1;
				}
				elsif ($$Fields[$ChildFields[$Ch1]] !~ /N/ and $$Fields[$ChildFields[$Ch2]] !~ /N/) # only want fully called
				{
					$$Fields[$FirstNewField+$Combo+1] = "D"; # kids have diff genoytope from father
					$Use = 1;
				}
				else
				{
					$$Fields[$FirstNewField+$Combo+1] = "N"; # not informative genotype from father
				}
				$$Fields[$FirstNewField+$Combo] = "N"; # not informative genotype from mother
				$Combo += 2;
			}
		}
	}
	elsif ($$Fields[$MatField] eq "11" and $$Fields[$PatField] eq "01") # allele 0 is unique for father
	{
		my $Combo = 0;
		for my $Ch1 (0..$NrChildren-2)
		{
			for my $Ch2 ($Ch1+1..$NrChildren-1)
			{
				if (($$Fields[$ChildFields[$Ch1]] =~ /0/ and $$Fields[$ChildFields[$Ch2]] =~ /0/) or
					($$Fields[$ChildFields[$Ch1]] eq "11" and $$Fields[$ChildFields[$Ch2]] eq "11")) # two kids same genotype, autosomes only for now
				{
					$$Fields[$FirstNewField+$Combo+1] = "S"; # kids have same genoytope from father
					$Use = 1;
				}
				elsif ($$Fields[$ChildFields[$Ch1]] !~ /N/ and $$Fields[$ChildFields[$Ch2]] !~ /N/) # only want fully called
				{
					$$Fields[$FirstNewField+$Combo+1] = "D"; # kids have diff genoytope from father
					$Use = 1;
				}
				else
				{
					$$Fields[$FirstNewField+$Combo+1] = "N"; # not informative genotype from father
				}
				$$Fields[$FirstNewField+$Combo] = "N"; # not informative genotype from mother
				$Combo += 2;
			}
		}
	}
	return ($Use, $Fields)
}

sub AnalyseState2 # parents 01 01
{
	my $Fields = shift; # ptr to testvariant fields
	my $Use = 0; # whether to use for state analysis
	my $Combo = 0; # counter for which field to enter data for, increments by 2 to allow 2 fields per combination, one maternal, one paternal
	#$$Fields[int(@$Fields)-1+($NrChildren*2)] = "";

	for my $Ch1 (0..$NrChildren-2) # loop for first child in pair
	{
		for my $Ch2 ($Ch1+1..$NrChildren-1) # loop for second child in pair
		{
			if (($$Fields[$ChildFields[$Ch1]] eq "00" and $$Fields[$ChildFields[$Ch2]] eq "00") or
				($$Fields[$ChildFields[$Ch1]] eq "11" and $$Fields[$ChildFields[$Ch2]] eq "11"))# two kids same genotype, autosomes only for now
			{
				$$Fields[$FirstNewField+$Combo] = "S"; # kids have same genoytope from mother, and
				$$Fields[$FirstNewField+$Combo+1] = "S"; # kids have same genoytope from father
				$Use = 1;
			}
			elsif (($$Fields[$ChildFields[$Ch1]] eq "00" and $$Fields[$ChildFields[$Ch2]] eq "11") or
				($$Fields[$ChildFields[$Ch1]] eq "11" and $$Fields[$ChildFields[$Ch2]] eq "00"))# two kids same genotype, autosomes only for now
			{
				$$Fields[$FirstNewField+$Combo] = "D"; # kids have different genoytope from mother, and
				$$Fields[$FirstNewField+$Combo+1] = "D"; # kids have different genoytope from father
				$Use = 1;
			}
			else
			{
				$$Fields[$FirstNewField+$Combo] = "N"; # not informative
				$$Fields[$FirstNewField+$Combo+1] = "N"; # not informative
			}
			$Combo += 2; # go to next set of fields for next child
		}
	}
	return ($Use, $Fields)
}

sub GetMIEBinStats
{
	# chr array and blocks hash are globals
	my $MIEB_TotSq = my $MIEB_Tot = my $MIEB_N = 0;
	foreach my $Chr (@ChrNames)
	{
		if (defined $MIE_Blocks{$Chr})
		{
			foreach my $Bin (keys %{$MIE_Blocks{$Chr}})
			{
				$MIEB_N++;
				$MIEB_Tot += $MIE_Blocks{$Chr}{$Bin};
				$MIEB_TotSq += $MIE_Blocks{$Chr}{$Bin}**2;
			}
		}
	}
	my $MIEB_Mean = eval{$MIEB_Tot/$MIEB_N};
	my $MIEB_SD = eval{(($MIEB_TotSq - $MIEB_Tot*$MIEB_Mean)/($MIEB_N-1))**0.5};
	my $LTail = $MIEB_Mean - $MIEB_SD;
	my $RTail = $MIEB_Mean + $MIEB_SD;
	return ($MIEB_N, $MIEB_Mean, $MIEB_SD, $LTail, $RTail);
}

sub SaveMIEBins
{
	my $P = shift; # ptr to hash of MEI blocks
	my $FH = shift; # output file
	my $LTail = shift; # threshold for left tail of distribution
	my $RTail = shift; # threshold for right tail of distribution
	my $RGB; # colour
	foreach my $Chr (@ChrNames) # each chr in turn, order defined when loading array
	{
		if (defined $MIE_Blocks{$Chr}) # catches eg Y in females, chr subsets of data
		{
			foreach my $Bin (sort {$a <=> $b} keys %{$$P{$Chr}})
			{
				if    ($$P{$Chr}{$Bin} < $LTail) {$RGB = $Colours{PaleBlue};} # low MIE bins
				elsif ($$P{$Chr}{$Bin} > $RTail) {$RGB = $Colours{Blue};} # high MIE bins
				else  {$RGB = $Colours{UncertClr};} # remainder
				print $FH $Chr,"\t",$Bin*$MIEBlockWindow,"\t",($Bin+1)*$MIEBlockWindow,"\t\t",$MIE_Blocks{$Chr}{$Bin},"\t";
				print $FH "+","\t",$Bin*$MIEBlockWindow,"\t",($Bin+1)*$MIEBlockWindow,"\t",$RGB,"\n";
			}
		}
	}
}

sub ConcatenateVariants
{
    my $ArrayIn = shift; # ptr to array
    my $StateFieldNr = shift; # field to process
    #print int(@$ArrayIn),"\n";
    my @ArrayOut; # array to store records out
    my $Nr = -1;
    foreach my $Entry (@$ArrayIn)
    {
		if ($$Entry[$StateFieldNr] ne "N")  # S or D entry present for chilren for parent genome represented by field $StateFieldNr
		{
			if ($Nr == -1)
			{
				$ArrayOut[++$Nr] = NewStateRecord (); # increment counter and attached record
				LoadStateRecord ($ArrayOut[$Nr], $Entry, $StateFieldNr);
			}
			elsif ($$Entry[$StateFieldNr] eq $ArrayOut[$Nr]->{State}) # state is same as current record
			{
				#print join("\t",@$Entry),"\n";
				#print "$ArrayOut[$Nr]->{State}\t$ArrayOut[$Nr]->{Chr}\t$ArrayOut[$Nr]->{Begin}\t$ArrayOut[$Nr]->{End}\t$ArrayOut[$Nr]->{Records}\n";
				#print "$$Entry[2]\t$ArrayOut[$Nr]->{End}\n"; exit;
				if ($$Entry[2] - $ArrayOut[$Nr]->{End} < $GapNot2Cross) # distance between blocks is less than 1Mb
				{
					$ArrayOut[$Nr]->{End} = $$Entry[3]; # update of end of state range
					$ArrayOut[$Nr]->{Records}++; # record added to count
				}
				else
				{
					$ArrayOut[++$Nr] = NewStateRecord (); # increment counter and attached record
					LoadStateRecord ($ArrayOut[$Nr], $Entry, $StateFieldNr);
				}
			}
			else # state is different to current state record
			{
				$ArrayOut[++$Nr] = NewStateRecord (); # increment counter and attached record
				LoadStateRecord ($ArrayOut[$Nr], $Entry, $StateFieldNr);
			}
		}
    }
    return \@ArrayOut; # return ptr to array
}

sub LoadStateRecord
{
	my $Out = shift;
	my $In = shift;
	my $StateFieldNr = shift;

				$Out->{State} = $$In[$StateFieldNr]; # get state for new record
				$Out->{Chr} = $$In[1]; # get chr
				$Out->{Begin} = $$In[2]; # get begin of state range
				$Out->{End} = $$In[3]; # get current end of state range
				$Out->{Records}++; # record added to new count
}

sub SaveStateBlocks
{
    my $FH = shift;
    my $Array = shift;

    #print $FH "Chr\tBegin\tEnd\tState\tLength\tRecords\n";
    for my $Entry (@$Array)
    {
		my $P = $Entry;
		print $FH $P->{Chr},"\t",$P->{Begin},"\t",$P->{End},"\t",$P->{State},"\t",($P->{End} - $P->{Begin}),"\t",$P->{Records},"\n";
    }
}

sub SaveStateBlocksAsBed
{
    my $FH = shift;
    my $Array = shift;

    for my $Entry (@$Array)
    {
		my $P = $Entry;
		my $RGB;
		if ($P->{State} eq "S")    {$RGB = $Colours{SameClr};}
		elsif ($P->{State} eq "D") {$RGB = $Colours{DiffClr};}
		elsif ($P->{State} eq "U") {$RGB = $Colours{UncertClr};}
		elsif ($P->{State} eq "G") {$RGB = $Colours{GapClr};}
		else 					   {$RGB = $Colours{PaleGrey};}

		print $FH $P->{Chr},"\t",$P->{Begin},"\t",$P->{End},"\t",$P->{State},"\t",$P->{Records},"\t";
		print $FH "+","\t",$P->{Begin},"\t",$P->{End},"\t",$RGB,"\n";
    }
}

sub ResolveSingletons
{
    my $ArrayIn = shift;
    my @ArrayOut = ();
    my $Out = -1;
    my $In = -1;
	my $Count = 0; # count nr singletons that are not Us
    # Either ignore singletons that lie between 2 blocks of reasonable nr variants,
    # or mark if start of short blocks of uncertainty

    # Load first record_in
    $ArrayOut[++$Out] = $$ArrayIn[++$In];

    # Process remainder of records_in
    for my $In (1..int(@$ArrayIn)-2)
    {
		if (($ArrayOut[$Out]->{State} eq $$ArrayIn[$In]->{State}) and
			($$ArrayIn[$In]->{Begin} - $ArrayOut[$Out]->{End} < $GapNot2Cross) ) # same state, so extend old block to include new block
		{
			$ArrayOut[$Out]->{End} = $$ArrayIn[$In]->{End};
			$ArrayOut[$Out]->{Records} += $$ArrayIn[$In]->{Records};
			$ArrayOut[$Out]->{StateErrors} += $$ArrayIn[$In]->{StateErrors}; # increment count of state errors for prev block
		}
		else # diff state
		{
			if ($$ArrayIn[$In]->{Records} == 1) # new rec is a singleton
			{
				if ($$ArrayIn[$In]->{State} ne "U") # is not a U
				{
					$Count++;
					if ($ArrayOut[$Out]->{State} eq $$ArrayIn[$In+1]->{State}) # prev block and next block the same
					{
						if ($$ArrayIn[$In+1]->{Records} > $ShortBlock) # next block > cutoff, add singleton to prev block
						{
							$ArrayOut[$Out]->{End} = $$ArrayIn[$In]->{End}; # increment length of prev block
							$ArrayOut[$Out]->{Records} += $$ArrayIn[$In]->{Records}; # increment count of recs in prev block
							$ArrayOut[$Out]->{StateErrors}++; # increment count of state errors for prev block
						}
						else # next block ² cutoff, so save singleton and make it a U
						{
							$ArrayOut[++$Out] = $$ArrayIn[$In];
							$ArrayOut[$Out]->{State} = "U"; # state is uncertain
							$ArrayOut[$Out]->{StateErrors}++; # store nr state errors for this block
							# save singleton to a bad state file? not done at the moment
						}
					}
				}
				else # singleton is a U, so just save it
				{
					$ArrayOut[++$Out] = $$ArrayIn[$In];
				}
			}
			else # not a singleton, save it
			{
				$ArrayOut[++$Out] = $$ArrayIn[$In];
				$ArrayOut[$Out]->{Length} = $ArrayOut[$Out]->{End} - $ArrayOut[$Out]->{Begin}; # clean up current rec
			}
		}
    }

	# process last rec
	$In = int(@$ArrayIn) - 1;
	my $P = $$ArrayIn[$In];
	#print "Last: $P->{State}\t$P->{Begin}\t$P->{End}\t$P->{Records}\t$P->{StateErrors}\n"; #exit;
	if ($$ArrayIn[$In]->{Records} == 1) # new rec is a singleton
	{
		if ($$ArrayIn[$In]->{State} ne "U") # is not a U
		{
			# add to prev
			#print "Old Out: $P->{State}\t$P->{Begin}\t$P->{End}\t$P->{Records}\t$P->{StateErrors}\n"; #exit;
			$ArrayOut[$Out]->{End} = $$ArrayIn[$In]->{End}; # increment length of prev block
			$ArrayOut[$Out]->{Records} += $$ArrayIn[$In]->{Records}; # increment count of recs in prev block
			$ArrayOut[$Out]->{StateErrors}++; # increment count of state errors for prev block
			my $P = $ArrayOut[$Out];
			#print "New Out: $P->{State}\t$P->{Begin}\t$P->{End}\t$P->{Records}\t$P->{StateErrors}\n"; #exit;
		}
		else
		{
			# save rec
			$ArrayOut[++$Out] = $$ArrayIn[$In];
		}
	}
	else
	{
		# save rec
		$ArrayOut[++$Out] = $$ArrayIn[$In];
	}

    return ($Count, \@ArrayOut);
}

sub ResolveShortBlocks
{
    my $ArrayIn = shift;
    my $ShortRun = shift; # length of runs to resolve
    my @ArrayOut = ();
    my $In = 0;
	my $Count = 0; # count nr singletons that are not Us

    # Process records_in situ, finding short runs and making uncertain
	# Process first rec
	if ($$ArrayIn[$In]->{Records} <= $ShortRun) # if small block, make U
	{
		#print "$$ArrayIn[$In]->{Records}\t$$ArrayIn[$In]->{Begin}\t$$ArrayIn[$In]->{End}\t$$ArrayIn[$In]->{State}\n";
		$$ArrayIn[$In]->{State} = "U"; # state is uncertain, restated in case first block was not U, but was short
		$$ArrayIn[$In]->{End} = $$ArrayIn[$In+1]->{Begin};
		#exit;
	}

	# Loop through all but last rec
    for my $In (1..int(@$ArrayIn)-2) # first record is header
    {
		if ($$ArrayIn[$In]->{Records} <= $ShortRun) # small block, look at it
		{
			$$ArrayIn[$In]->{State} = "U"; # state is uncertain, restated in case first block was not U, but was short
			$$ArrayIn[$In]->{Begin} = $$ArrayIn[$In-1]->{End};
			$$ArrayIn[$In]->{End} = $$ArrayIn[$In+1]->{Begin};
		}
    }

	# Process last rec
	$In = int(@$ArrayIn) - 1;
	if ($$ArrayIn[$In]->{Records} <= $ShortRun) # if small block, make U
	{
		$$ArrayIn[$In]->{State} = "U"; # state is uncertain, restated in case first block was not U, but was short
	}

	# Join contiguous U recs
    my $Out = $In = 0;
    $ArrayOut[$Out] = $$ArrayIn[$In];
    # Process remainder of records_in
    $In = 0;
    for $In (1..int(@$ArrayIn)-1)
    {
		if ($ArrayOut[$Out]->{State} eq "U" and $$ArrayIn[$In]->{State} eq "U") # both blocks are U, join
		{
			$ArrayOut[$Out]->{End} = $$ArrayIn[$In]->{End};
			$ArrayOut[$Out]->{Records} += $$ArrayIn[$In]->{Records};
			$ArrayOut[$Out]->{StateErrors} += $$ArrayIn[$In]->{StateErrors}; # increment count of state errors for prev block
			$Count++;
		}
		else # save
		{
			$ArrayOut[++$Out] = $$ArrayIn[$In] if $$ArrayIn[$In]->{End} > 0; # save, but only if a real rec, catch all for overrun
		}
    }
    return ($Count, \@ArrayOut);

}

sub FillInGaps
{
    my $ArrayIn = shift;
    my @ArrayOut = ();

    my $Out = -1;
    my $In = -1;
    my $NrIns = int(@$ArrayIn);
	# Only need to put gaps between S and D blocks
    while ($In < $NrIns - 1)
    {
		$ArrayOut[++$Out] = $$ArrayIn[++$In]; # loading current input record into output record
		if ($$ArrayIn[$In]->{State} =~ /S|D/ and $$ArrayIn[$In+1]->{State} =~ /S|D/) #  if one is S and other is D
		{
			my $P = NewStateRecord ();
			$P->{Chr} = $$ArrayIn[$In]->{Chr}; # cur chr
			$P->{State} = "G"; # gap, diff to other uncertain bits in having no recs
			$P->{Records} = 0; # filler, no records
			$P->{Begin} = $$ArrayIn[$In]->{End}; # begin is end of out rec
			$P->{End} = $$ArrayIn[$In+1]->{Begin}; # end is begin of in rec
			$P->{Length} = $P->{End} - $P->{Begin}; # recalc length
			$ArrayOut[++$Out] = $P; # save new gap record
		}
		else # save current in rec to out
		{
			#$ArrayOut[++$Out] = $$ArrayIn[++$In]; # loading current input record into output record
		}
    }
	#$ArrayOut[++$Out] = $$ArrayIn[++$In]; # save last rec
    return \@ArrayOut;
}

#sub FillInGaps
#{
#    my $ArrayIn = shift;
#    my @ArrayOut = ();
#    my $Out = 0;
#    my $In = 0;
#
#    # Process records_in to records_out
#    # Load first record_in - script looks back 1 rec so need it loaded
#    $ArrayOut[$Out] = NewStateRecord (); # put in to first rec out
#    LoadNewRecord ($$ArrayIn[$In], $ArrayOut[$Out]);
#    if ($$ArrayIn[$In]->{State} eq "U") # if first rec is U
#    {
#	$$ArrayIn[$In]->{End} = $$ArrayIn[$In+1]->{Begin}; # extend end to cover gap to next rec
#	$$ArrayIn[$In]->{Length} = $$ArrayIn[$In]->{End} - $$ArrayIn[$In]->{Begin}; # recalc length
#    }
#
#    # Remainder of recs
#    my $NrIns = int(@$ArrayIn);
#    while ($In < $NrIns)
#    {
#	if ($$ArrayIn[$In]->{State} eq "U") # if rec is U
#	{
#	    $$ArrayIn[$In]->{End} = $$ArrayIn[$In+1]->{Begin}; # extend end to cover gap to next rec
#	    $$ArrayIn[$In]->{Length} = $$ArrayIn[$In]->{End} - $$ArrayIn[$In]->{Begin}; # recalc length
#
#	    if ($$ArrayIn[$In-1]->{State} eq "U") # if prev rec is U, concat recs
#	    {
#		$ArrayOut[$Out]->{End} = $$ArrayIn[$In]->{End}; # extend prev rec to end of this rec
#		$ArrayOut[$Out]->{Records} += $$ArrayIn[$In]->{Records}; # add nr records
#		$ArrayOut[$Out]->{Length} = $ArrayOut[$Out]->{End} - $ArrayOut[$Out]->{Begin}; # recalc length
#		# nb not saving this rec, just adding to previous
#	    }
#	    else # prev rec not a U, save this rec
#	    {
#		$ArrayOut[++$Out] = NewStateRecord ();
#		LoadNewRecord ($$ArrayIn[$In], $ArrayOut[$Out]);
#	    }
#	}
#	else # rec not a U
#	{
#	    $ArrayOut[++$Out] = NewStateRecord ();
#	    LoadNewRecord ($$ArrayIn[$In], $ArrayOut[$Out]);
#	    if ($$ArrayIn[$In+1]->{State} eq "U") # if next rec is U, back extend next rec
#	    {
#		$$ArrayIn[$In+1]->{Begin} = $$ArrayIn[$In]->{End};
#		$$ArrayIn[$In+1]->{Length} = $$ArrayIn[$In]->{End} - $$ArrayIn[$In+1]->{Begin}; # recalc length
#	    }
#	    else # next rec not a U, need to fill gap with new rec
#	    {
#		$ArrayOut[++$Out] = NewStateRecord ();
#		$ArrayOut[$Out]->{Chr} = $$ArrayIn[$In]->{Chr}; # cur chr
#		$ArrayOut[$Out]->{State} = "U"; # uncertain, diff to others in no recs
#		$ArrayOut[$Out]->{Records} = 0; # filler, no records
#		$ArrayOut[$Out]->{Begin} = $$ArrayIn[$In]->{End}; # begin is end of curr rec
#		$ArrayOut[$Out]->{End} = $$ArrayIn[$In+1]->{Begin}; # end is begin of next rec
#		$ArrayOut[$Out]->{Length} = $ArrayOut[$Out]->{End} - $ArrayOut[$Out]->{Begin}; # recalc length
#		$ArrayOut[$Out]->{Records} = 0; # filler, no records
#	    }
#	}
#	$In++; # move counter to next rec_in
#    }
#    return \@ArrayOut;
#}

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

sub LoadNewRecord
{
    my $In = shift;
    my $Out = shift;
    $Out->{Chr} = $In->{Chr};
    $Out->{State} = $In->{State};
    $Out->{Begin} = $In->{Begin};
    $Out->{End} = $In->{End};
    $Out->{Records} = $In->{Records};
}

sub NewStateRecord
{
    my $Record =
    {
	Chr => "",
        Begin => -1,
        End => -1,
		State => "",
        Records => 0,
		MIEs => 0,
		StateErrors => 0,
        Length => -1,
    };
    return $Record;
}

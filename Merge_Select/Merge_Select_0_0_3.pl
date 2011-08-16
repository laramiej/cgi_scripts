#!/usr/bin/perl
use strict;
use File::Basename;
$| = 1;
#use feature "state";
my $Debug = 0;

# Version 0_0_3

# Merge 2 files based on matching fields
# Select matched, unmatched or all records of file A to output
# Choose whether to output fields present in both files (duplicate fields)
# Choose fields of both files to output

# perl Merge_Select \
# --Input_FileA input_file_1 \
# --Input_FileB input_file_2 \
# --Output_File filename \
# --Select	all|matched|unmatched
# --Match fieldname_in_fileA:fieldname_in_FileB [repeated]
# --Remove_Dup_Fields [param present or absent]
# --Output_Fields A.*,B.* [or eg A.fieldname1,A.fieldname2,B.fieldname1 etc]

# eg
# perl /Users/rtearle/Documents/Programming/Perl/Scripts/Dev/Merge_Select.pl \
# --Input_FileA /Users/rtearle/Documents/Scratch/CEPH_Family_3Gen_Dominant_1/CEPH_Family-MGM_F_V_C1-Het_Het_Het-Coding-Protein_Affecting-Edited.tsv \
# --Input_FileB  /Users/rtearle/Documents/TBF/CG_Public_Genomes/54_Genomes/Public_Genomes_Unrelated_54-testvariants.tsv_Freq_Short.tsv.bz2 \ \
# --Output_File /Users/rtearle/Documents/TBF/ \
# --Select	all
# --Match chromosome:chromosome
# --Match begin:begin
# --Match end:end
# --Match varType:varType \
# --Match alleleSeq:alleleSeq \
# --Remove_Dup_Fields
# --Output_Fields A.*,B.*


# Parsing and storing input parameters
# Only input_file fields can be repeated
# input paramaters are case insensitive

print "\nProcessing input parameters\n";
my $NrParams;
my %ExpectedParams =  GetExpectedParams (); # list of expected parms
my %EnteredParams = GetEnteredParams (); # list of entered params
TestParameters (%EnteredParams); # compares expectd and entered params

# Assigning input params
# Input Files
my $FileA = $EnteredParams{Input_FileA}; # ptr to list of input files
unless (-f $FileA) {die "Input file $FileA not found\n";} # requires existing files
my $FileB = $EnteredParams{Input_FileB}; # ptr to list of input files
unless (-f $FileB) {die "Input file $FileB not found\n";} # requires existing files
# Ouput File
my $FileOut = $EnteredParams{Output_File}; # output file
if (-f $FileOut) # ouput file exists, create a new one based on the name
{
	print "Output file $FileOut exists, modifying to unique file name ";
	$FileOut =~ /^(.+?)\./; my $Stub = $1; # set stub to name without extensions
	$FileOut =~ /(\..+)?$/; my $Ext = $1; # set ext to extensions
	my $n = 1; # n will increment to find a unique name
	my $Suff = ""; # suff tracks n
	while (-f $Stub.$Suff.$Ext) {$Suff = "-$n"; $n++;} # loop till we have a new unique filename
	$FileOut = $Stub.$Suff.$Ext; # file out now has same name, same extensions, but also -n at the end of the name, making it unique
	print "$FileOut\n";
}
# Match fields
my %MatchFields;
foreach my $Val (@{$EnteredParams{Match}})
{
	$Val =~ /(.+):(.+)/; # split entry based on :
	unless ($1 and $2) {die "Field names to match must be of the form name:name; entry $Val not valid\n";}
	$MatchFields{$1} = $2; # make a hash with first part, second part from match entry
}
my $NrMatchFields = int %MatchFields;
# Select
my $Select = $EnteredParams{Select}; # can output all|merged|unmerged
unless ($Select =~ /^(all|matched|unmatched)$/) {die "You must choose to output all, matched or unmatched records from file A";}
# Remove Duplicate Fields ie if file B has fields already in file A
my $RemoveDupFields = $EnteredParams{Remove_Dup_Fields}; # remove dup fields in B?
#if ($RemoveDupFields and $RemoveDupFields !~ /^(A|B)$/) {die "If you must choose to remove duplicate fields you must select either A or B to remove; entry $RemoveDupFields not valid\n";}
my $OuputFieldString = $EnteredParams{Output_Fields}; # can output all|merged|unmerged

print "Input Files:\n\t$FileA\n\t$FileB\n";
print "Ouput File: $FileOut\n";
print "Match Fields: \n"; foreach my $Field (keys %MatchFields) {print "\t$Field $MatchFields{$Field}\n";};
print "Records to output: $Select\n";
print "Remove duplicate fields: ";{$RemoveDupFields == 1 ? print "true" : print "false"}; print "\n";

# Loading chr nrs, setting up var hash
my @ChrNames = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
				'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
				'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'); # using this array forces the output order of chrs into the correct order
my %FileARecs; # hash for storing file A recs
foreach my $Chr (@ChrNames) {$FileARecs{$Chr} = {};} # print "$Chr\t";}  # set up hash of hashes, one for each chr

# Open file, remove header, get field names and any start char, leave file open ready to process data, returns FH
my $FieldsA = []; my $FieldsB = []; # arrays of fields in file A and B
my ($StartCharA, $StartCharB); # start char for column line, if any, to store them
my $IN_A = GetFieldNames ($FileA, $FieldsA, $StartCharA);
my $IN_B = GetFieldNames ($FileB, $FieldsB, $StartCharB);

# Find field nr for chr, begin, end fields
my $ChrFieldNrA = FindFieldNrByName ($FieldsA, "chr|chrom|chromosome"); # field nr containing chr
if ($ChrFieldNrA == -1) {die "Cannot find a chr field called 'chr', 'chrom' or 'chromosome' in file $FileA\n";}
my $ChrFieldNrB = FindFieldNrByName ($FieldsB, "chr|chrom|chromosome"); # field nr containing chr
if ($ChrFieldNrB == -1) {die "Cannot find a chr field called 'chr', 'chrom' or 'chromosome' in file $FileA\n";}
my $BeginFieldNrA = FindFieldNrByName ($FieldsA, "begin|offset"); # field nr containing chr
if ($ChrFieldNrA == -1) {die "Cannot find a begin field called 'begin', or 'offset' in file $FileA\n";}
my $EndFieldNrA = FindFieldNrByName ($FieldsA, "end"); # field nr containing chr
if ($ChrFieldNrA == -1) {print " Warning: cannot find an end field called 'end' in file $FileA\n";}

# Find field numbers for match fields
my ($MatchFieldsNrsA, $MatchFieldsNrsB) = GetMatchFieldNumbers (\%MatchFields, $FieldsA, $FieldsB);
$NrMatchFields = int@$MatchFieldsNrsA-1;
#foreach my $n (0..$NrMatchFields-1) {print $$MatchFieldsNrsA[$n],"\t",$$MatchFieldsNrsA[$n],"\n";}

# Find field numbers of duplicate fields if remove dups flag is set
my $DupFieldNrsB;
if ($RemoveDupFields) {$DupFieldNrsB = GetDupFieldNumbers ($FieldsA, $FieldsB);} #foreach my $n (@$DupFieldNrsB) {print "$n\n";}; exit;

# Get field nrs for customer selected output fields
my $OutputFields = ExtractOutputFields ($FieldsA, $FieldsB, $OuputFieldString, $RemoveDupFields, $DupFieldNrsB);

# Load File A in to hash, create key from match fields, store whole rec as val struct
print "\nLoading file $FileA\n";

while (<$IN_A>) # loop through recs
{
	chomp $_; # remove trailing return
	my @Fields = split "\t", $_; # get fields as array
	my $Key = ""; # key for hash
	foreach my $FieldNr (@$MatchFieldsNrsA) {$Key .= $Fields[$FieldNr]."\t";} # add match field to key
	$Key =~ s/\t$//; # remove trailing tab
	$FileARecs{$Fields[$ChrFieldNrA]}->{$Key} = NewRec (); # val is a struct
	$FileARecs{$Fields[$ChrFieldNrA]}->{$Key}->{rec} = $_; # store rec
	$FileARecs{$Fields[$ChrFieldNrA]}->{$Key}->{begin} = $Fields[$BeginFieldNrA]; # store begin
	$FileARecs{$Fields[$ChrFieldNrA]}->{$Key}->{end} = $Fields[$EndFieldNrA]; # store end
	print "FileA: $Key\n$FileARecs{$Fields[$ChrFieldNrA]}->{$Key}->{rec}\n" if $Debug; #exit;
	last if $Debug;
}

# Process File B
print "\nComparing to file $FileB\n";
while (<$IN_B>) # loop through recs
{
	chomp (my $Rec = $_); # assign rec, remove trailing return (use rec througout loop)
	my @Fields = split "\t", $Rec; # get fields
	my $Key = ""; # key for comparing to keys of file A
	foreach my $FieldNr (@$MatchFieldsNrsB) {$Key .= $Fields[$FieldNr]."\t";} # add field to key
	$Key =~ s/\t$//; # remove trailing tab
	#print "$Key\n$Rec\n" unless $Debug++; #exit;

	# Match to keys in file A
	if ($FileARecs{$Fields[$ChrFieldNrB]}->{$Key}) # key for this rec present in file A
	{
		next if $FileARecs{$Fields[$ChrFieldNrB]}->{$Key}->{match} == 1; # skip if rec in FileA is already matched by rec in File B
		$FileARecs{$Fields[$ChrFieldNrB]}->{$Key}->{match} = 1; # rec in FileA is matched by rec in FileB
		next if $Select =~ /unmatched/; # no point modifying matches if we are not keeping them

		# add file B fields to file A rec $OutputFields
		my $PA = $FileARecs{$Fields[$ChrFieldNrB]}->{$Key}; # ptr to hash entry for this rec,
		$PA->{rec} = BuildOuputString($PA->{rec}, $Rec, $OutputFields);
		print "\n\nIn FileB: key: $Key\nFile A: $PA->{rec}\nFileB: $Rec\n if $Debug"; #exit;
	}
}

# Select and print desired recs
open my $OUT, ">", $FileOut; # open output file
my $TabsForFieldsB = "\t" x (int(@$FieldsB)-1); # to add tabs for fields B to unmatched recs
my %Class = {"M" => 0, "N" => 0, "A" => 0, "MP" => 0, "NAP" => 0, "NUP" => 0};

# Output header
print $OUT BuildOuputString(join("\t",@$FieldsA), join("\t",@$FieldsB), $OutputFields)."\n";

foreach my $Chr (@ChrNames)
{
	#print "$FileARecs{$Chr}\n";
	foreach my $Key (sort {SortStringsasArrays ($FileARecs{$Chr}, $a, $b)} keys %{$FileARecs{$Chr}}) # using sub to sort on being, end fields)
	{
		$Class{A}++;
		#print "$Chr\n";
		print "In compare: Key: $Key\nRec: $FileARecs{$Chr}->{$Key}->{rec}\nmatch: $FileARecs{$Chr}->{$Key}->{match}\n\n" if $Debug;
		if ($FileARecs{$Chr}->{$Key}->{match}) # rec is matched, so has already been processed
		{
			$Class{M}++;
			if ($Select eq "all" or $Select eq "matched") # want all or matched
			{
				print $OUT $FileARecs{$Chr}->{$Key}->{rec},"\n";
				$Class{MP}++;
				#$Class{PA}++;
			}
		}
		else # rec is not matched
		{
			$Class{N}++;
			if ($Select eq "all") # want all
			{
				$FileARecs{$Chr}->{$Key}->{rec} = BuildOuputString($FileARecs{$Chr}->{$Key}->{rec}, "", $OutputFields); # get selected fields for rec from file A, no fiels for B as not matched
				print $OUT $FileARecs{$Chr}->{$Key}->{rec},$TabsForFieldsB,"\n"; # adding empty fields for file B, as no match
				$Class{UAP}++;
			}
			elsif ($Select eq "unmatched") # want unmatched
			{
				$FileARecs{$Chr}->{$Key}->{rec} = BuildOuputString($FileARecs{$Chr}->{$Key}->{rec}, "", $OutputFields); # get selected fields for rec from file A, no fiels for B as not matched
				print $OUT $FileARecs{$Chr}->{$Key}->{rec},"\n"; # just fields of A in output
				#$Class{PA}++;
				$Class{UUP}++;

			}
		}
	}
}

foreach my $Class (keys %Class) {print "$Class\t$Class{$Class}\n";}

###########################################################################
#                                   SUBS                                  #
###########################################################################

sub GetExpectedParams
{
	my %Hash = # hash to store expected params
	(
		"Input_FileA" => -1,
		"Input_FileB" => -1,
		"Output_File" => -1,
		"Select" => -1,
		"Match" => [],
		"Remove_Dup_Fields" => 0,
		"Output_Fields" => "A.*,B.*",
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
		chomp ($ARGVs[$n]);
		my ($Key, $Val) = split / /, $ARGVs[$n], 2; # put first element into key, any other elements into val
		if ($Key eq "Match") # multiple entries expected, setting up array
		{
			push @{$Hash{$Key}}, $Val; # add input to input hash
		}
		elsif ($Key eq "Remove_Dup_Fields") # multiple entries expected, setting up array
		{
			$Hash{$Key} = 1; # set it to true
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
		if (!defined $EnteredParams{$ExpArg}) #  if there is no entry for this expected param
		{print "Parameter --$ExpArg appears to be missing from arguments\n"; $ArgErrors++;}
		elsif ($EnteredParams{$ExpArg} == -1) # if no entry has been entered from the input list
		{print "No value for parameter --$ExpArg in arguments"; $ArgErrors++;}
	}
	if (int keys %EnteredParams > $NrParams) {die "There are too many parameter in the arguments";} # die if too many params
	if ($ArgErrors) {die "$ArgErrors detected in input parameters";} # die if a param not entered
}

sub GetFieldNames
{
	my $File = shift;
	my $Fields = shift;
	my $StartChar = shift;

	my $ColHeader;
	my $IN = OpenFile ($File); # open file with correct file format
	my $Header = GetHeader ($IN); # get header, not keeping

	if ($Header) # is a CG file
	{
		my $ColHeader = <$IN>; # get col header, is next line
	}
	else # not a CG file
	{
		close $IN;
		$IN = OpenFile ($File); # close and reopen file, ie start again
		$ColHeader = <$IN>; # get col header, should be first line
	}

	if ($ColHeader =~ /^(>|<|#)/) # remove a starting char if is >, < or #
	{
		$StartChar = $1; # start char saved
		$ColHeader =~ s/^.//; # and removed from col header
	}
	#print "$ColHeader";

	chomp $ColHeader; # remove trailing return
	@$Fields = split "\t", $ColHeader; # split field names on tab and save

	return $IN; # returning file handle
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

sub GetMatchFieldNumbers
{
	my $MatchFields = shift;
	my $FieldsA = shift;
	my $FieldsB = shift;
	my $MatchFieldsNrsA = [];
	my $MatchFieldsNrsB = [];

	foreach my $Field (keys %$MatchFields)
	{
		my $Found = 0;
		for my $n (0..int(@$FieldsA)-1)
		{
			if ($$FieldsA[$n] eq $Field)
			{
				push @{$MatchFieldsNrsA}, $n;
				$Found = 1; last; # keeping first match
			}
		}
		unless ($Found) {die "$Field is a not a field in $FileA\n";} #

		$Found = 0;
		my $MatchNrB = -1;
		for my $n (0..int(@$FieldsB)-1)
		{
			if ($$FieldsB[$n] eq $MatchFields{$Field})
			{
				push @{$MatchFieldsNrsB}, $n;
				$Found = 1; last; # keeping first match
			}
		}
		unless ($Found) {die "$Field is a not a field in $FileB\n";} #
	}
	return ($MatchFieldsNrsA, $MatchFieldsNrsB); # returning ptr to hash containing fieldnrs to match
}

sub GetDupFieldNumbers
{
	my $FieldsA = shift;
	my $FieldsB = shift;
	my $DupFieldNrsB = []; # array to store fields to delete

	for my $n (0..int(@$FieldsA)-1)
	{
		for my $m (0..int(@$FieldsB)-1)
		{
			if ($$FieldsA[$n] eq $$FieldsB[$m])
			{
				if ($RemoveDupFields eq "A") {push @$DupFieldNrsB, $n;} # only keeping first match, not as clean but very clear
				else 					{push @$DupFieldNrsB, $m;}
			}
		}
	}
	return $DupFieldNrsB; # returning ptr to array containing fieldnrs to not output
}

sub ExtractOutputFields
{
	my %AllFields;
	$AllFields{A} = shift;
	$AllFields{B} = shift;
	my $OuputFieldString = shift;
	my $RemoveDupFields = shift;
	my $DupFieldNrsB = shift;

	my %Keep; #$Keep{A} = []; $Keep{B} = []; # hash of fields to keep
	my %KeepNrs;
	$KeepNrs{A} = []; # array to store nrs of fields to keep for file A
	$KeepNrs{B} = []; # and B

	# Dissect and check output field string
	my @Bits1 = split /,/, $OuputFieldString; # break string on ','
	foreach my $Bit (@Bits1) # processing each part
	{
		my ($File, $Field) = split /\./, $Bit; # split part on '.'
		#print "Bit: $Bit\tFile: $File\tField: $Field\n";
		unless ($File =~ /^(A|B)$/) {die "$Bit is not a field reference in the correct form, needs to be eg 'A.*' or 'A.chromosome'\n";}
		#push @{$Keep{$File}}, $Field; # add to
		$Keep{$File} .= $Field."\t";
	}

	# Extract fields to keep
	foreach my $File ("A", "B") # for each file
	{
		my $P = $AllFields{$File};  # ptr to array with full list of fields for this file
		$Keep{$File} =~ s/\t$//; # remove trailing tab
		#print "File: $File\tFields: $Keep{$File}\n";
		if ($Keep{$File} =~ /\*/)
		{
			#print "$File all\n";
			for my $n (0..int(@{$P})-1) # all field nrs for this file
			{
				push @{$KeepNrs{$File}}, $n; # keep all
			}
		}
		else
		{
			my @Fields = split /\t/, $Keep{$File};
			foreach my $Field (@Fields)
			{
				my $Found = 0;
				#print "$File $Field...\t";
				for my $n (0..int(@{$P})-1) # all field nrs for this file
				{
					#print "$$P[$n]\t"; sleep 1;
					if ($Field eq $$P[$n]) # in keep list
					{
						push @{$KeepNrs{$File}}, $n; # so keep
						#print "yes\t$n\t";
						$Found = 1;
						last;
					}
				}
				die "Field $Field was not found in the list of fields for file $File\n" unless $Found;
			}
		}
		# print "File: $File\t"; foreach my $n (@{$KeepNrs{$File}}) {print "$n ";}; print "\n";

		# Remove dups if flag set
		# print "$File\t$RemoveDupFields\n";
		if ($File eq "B" and $RemoveDupFields)
		{
			#print int(@$DupFieldNrsB),"\n"; exit;
			#for (my $n = int(@{$KeepNrs{$File}})-1; $n > 0; $n--) # all field nrs to be kept for this file backwards
			#print int(@{$KeepNrs{$File}})-1,"\n";
			for my $n (0..int(@{$KeepNrs{$File}})-1) # all field nrs to be kept for this file backwards
			{
				for my $m (0..int(@$DupFieldNrsB)-1) # field nrs in dup list
				{
					#print "${$KeepNrs{$File}}[$n]\t$$DupFieldNrsB[$m]\n"; sleep 1;
					if (${$KeepNrs{$File}}[$n] == $$DupFieldNrsB[$m])
					{
						#print "$n\t$m\t${$KeepNrs{$File}}[$n]\t$$DupFieldNrsB[$m]\n";
						$KeepNrs{$File}[$n] = -1;
					} # remove if same
				}
			}
			@{$KeepNrs{$File}} = grep ($_ != -1, @{$KeepNrs{$File}});

		}
		@{$KeepNrs{$File}} = sort {$a <=> $b} @{$KeepNrs{$File}}; # put them in to original order
		#print "File: $File\t"; foreach my $n (@{$KeepNrs{$File}}) {print "$n ";}; print "\n";
	}

	return (\%KeepNrs); # returning ptr to hash of arrays containing fieldnrs to
}

sub BuildOuputString
{
	my %Recs;
	$Recs{A} = shift; # string A
	$Recs{B} = shift; # string B
	my $FieldNrs = shift; # ptr to hash of arrays of field nrs to use
	#my $FullOutputString = shift; # ptr to hash, whether to use full output strings

	my %Fields; # hash to store fields for A, B as arrays
	my $String = ""; # string to output

#print "$A\n";
#print "$A\n";
#print "$Fields\n";
#print $FieldNrs,"\n",%$FieldNrs,"\t",$$FieldNrs{A},"\n";

	for my $File ("A", "B") # two files
	{
		$Fields{$File} = [];
		chomp $Recs{$File}; @{$Fields{$File}} = split "\t", $Recs{$File}; # put string in to array
		for my $n (0..int(@{$$FieldNrs{$File}})-1) # elements of field nr array
		{
			my $m = $$FieldNrs{$File}->[$n]; # get contents of element n
			$String .= ${$Fields{$File}}[$m]."\t"; # add contents of field m of string A|B into output string
		}
	}
	$String =~ s/\t$//;
	#print $String,"\n"; exit if $Count++ > 10;
	return $String;
}

sub FindFieldNrByName
{
	my $Fields = shift;
	my $SearchString = shift;

	for my $n (0..int(@$Fields)-1)
	{
		#print "$n\t$$Fields[$n]\n"; sleep 1;
		if ($$Fields[$n] =~ /$SearchString/)
		{
			return $n; # field nr containing chr
		}
	}
	return -1; # if no field found
}

sub SortStringsasArrays # sorts based on begin and end of two recs
{
	my $Hash = shift; # ptr to hash of recs
	my $Key1 = shift; # first string
	my $Key2 = shift; # second string


	# array[1] is begin, array[2] is end, returning order based on these fields
	if ($$Hash{$Key1}->{begin} < $$Hash{$Key2}->{begin}) # begin of 1 < begin of 2
	{
		return -1;
	}
	elsif ($$Hash{$Key1}->{begin} == $$Hash{$Key2}->{begin}) # begin of 1 == begin of 2
	{
		if ($$Hash{$Key1}->{end} < $$Hash{$Key2}->{end}) # end of 1 < end of 2
		{
			return -1;
		}
		elsif ($$Hash{$Key1}->{end} == $$Hash{$Key2}->{end}) # end of 1 == end of 2
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

sub NewRec
{
	my $Struct =
	{
		rec => "",
		match => 0,
		begin => -1,
		end => -1,
	};
	return $Struct;

}
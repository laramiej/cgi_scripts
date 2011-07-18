#!/usr/bin/perl
use strict;
use File::Basename;
$| = 1;

# Evidence Interval to BED file
# takes an evidence interval file, extracts data and formats for BED

# Format:
# perl prog file dir string
# ie
# perl EvInt2BED.pl evint_file output_dir name
# eg
# perl /Users/ricktearle/Documents/Programming/Perl/Scripts/Dev/EvInt2BED_0_1.pl \
# /Users/ricktearle/Documents/Training_Data/NA19240/GS00028-DNA_C01/ASM/EVIDENCE/evidenceIntervals-chr21-GS19240-180-36-ASM.tsv.bz2 \
# /Users/ricktearle/Documents/TBF  19240
# blocks coloured in greyscale based on interval score

# Setting variables from Args - expect snp_file output_dir
unless (@ARGV >= 2) {print "Need full name of evidence interval file,  directory for output and name to describe BED file\n"; exit};
my $FileIn = $ARGV[0];
unless (-f $FileIn) {print "Cannot find file $FileIn\n"; exit};
my $DirOut = $ARGV[1];
unless (-d $DirOut) {mkdir $DirOut or die "$!: Cannot find/create output directory $DirOut\n";}
my $Name;
if ($ARGV[2]) {$Name =  $ARGV[2];}
else          {$Name = "";}

my $FileOut = $FileIn;
$FileOut =~ s/^.+\///; # remove path
$FileOut =~ s/\..*//; # remove extensions
print "$FileOut\n";

$FileOut .= ".bed";

open (my $OUT, ">", "$DirOut/$FileOut") or die ("$!: can't open file $DirOut/$FileOut.tsv\n");


#### Evidence Int file ####
print "Opening evInt File\n";
my $IN = OpenFile ($FileIn);
my $Header = GetHeader ($IN);
unless (int @$Header) {die "File $FileIn does not appear to have a header";}
my $ColHeader = GetColHeader ($IN);
unless ($ColHeader) {die "File $FileIn does not appear to have a column header";}
    
# track name=pairedReads description="Clone Paired Reads" useScore=1
print $OUT "track\tname=$Name\tdescription=EvidenceIntervals\tuseScore=1\titemRgb=0\n";

while (<$IN>)
{
    
    #chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRGB, blockCount, blockSizes, blockStarts
    # not using itemRGB or greater fields
    # IntervalId Chromosome OffsetInChromosome Length Ploidy AlleleIndexes Score Allele0 Allele1 Allele2 Allele1Alignment Allele2Alignment
    my ($ID, $Chr, $P1, $P2, undef, $Alleles, $Score, $Allele0, $Allele1, undef) = split /\t/, $_, 10; # using offset and length
    $P2 += $P1;

    print $OUT "$Chr\t$P1\t$P2\t$ID\t$Score\t+\t$P1\t$P2\n"
    
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
	open ($FH, "zcat $File |") or die ("$!: can't open file $File");
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


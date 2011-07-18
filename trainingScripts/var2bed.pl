#!/usr/bin/perl
use strict;

# A Very Simple CG Variants to BED Converter

# A smarter one could color code by variant type or by score

my $usage=  "usage:   perl var2bed.pl [option] var_file_name chromosome start_value end_value\n".
	    "options: -r   output reference calls in addition to variants\n".
	    "         -s   add allele sequences to call name (does not apply with -r)\n";

# process -r optional arg which says we should output ref calls too
my $output_ref= 0;
my $output_call_types= 1;
my $add_allele_seqs= 0;

if ($ARGV[0] eq "-r") {
    $output_ref= 1;
    $output_call_types= 0;
    shift @ARGV;
} elsif ($ARGV[0] eq "-s") {
    $add_allele_seqs= 1;
    shift @ARGV;
} elsif ($ARGV[0]=~/^\-/) {
    die $usage;
}

# check that we got 4 args
die $usage if (!defined($ARGV[0]) || !defined($ARGV[1]) || !defined($ARGV[2]));

# arg #1 ($ARGV[0]) is the file name

# arg #2 is the chromosome
die "arg 2 '$ARGV[1]' must be a chromosome name\n$usage"
  unless ($ARGV[1]=~/^chr\d+$/ || $ARGV[1]=~/^chr[MXY]$/ || $ARGV[1]=~/^all$/);

# check that args 3-4 are all integers
die "arg 3 '$ARGV[2]' must be an integer\n$usage" unless ($ARGV[2]=~/^\d+$/);
# this can be defaulted to all
die "arg 4 '$ARGV[3]' must be an integer or not used\n$usage" unless ($ARGV[3]=~/^\d+$/ || !defined($ARGV[3]));

#if(!defined($ARGV[3])){
#	print "Setting stop to 999999999 to output all variants\n";
#}
# set values from args
my $file=    $ARGV[0];
my $chr=     $ARGV[1];
my $start=   $ARGV[2];
my $end=     $ARGV[3] || 999999999; 

# open file
if ($file=~/\.gz$/) { # if it ends with .gz
    open (FILE,"zcat $file |")  || die "can't open 'zcat $file'";
} elsif ($file=~/\.bz2$/) { # if it ends with .bz2
    open (FILE,"bzcat $file |") || die "can't open 'bzcat $file'";
} else {
    open (FILE,"<$file") || die "can't open '$file'";
}

# read var file until we get to a line which starts with >
while (<FILE>) {
    last if (/^\>/);
}

# print header line for BED file
print "track name=var2bed itemRgb=On\n";

# this subroutine getline allows us to get and ungets lines from FILE
my $unget= 0;
sub getline() {
    if ($unget) { 
		$unget= 0;
		return 1;
    } elsif(eof FILE) {
		return 0;
    } else {
		$_= <FILE>;
		return 1;
    }
}

# Here we parse the body of the var file, which looks like:
#
#    0       1      2      3      4      5         6      7          8         9           10     11
# >locus   ploidy allele chrom  begin   end     varType reference  alleleSeq  totalScore hapLink xRef
# 21030571   2      1    chr21  9719767 9719773   ref   AATTCT     AATTCT     112

# read loci
my $count = 0;
while (getline()) {	# read the 1st (possibly only) line of the locus
    chomp;
    $count++;
    
    my @column= split(/\t/); 	# split on tab into an array
    my $locus_id= $column[0];	# assign to variable names for clarity
    my $allele=   $column[2];
    my $type=     $column[6];
    my $seq=      $column[8];
    my $score=    $column[9];
    my $locus_ploidy= $column[1];		
    my $locus_chrom=  $column[3];		
    my $locus_start=  $column[4];
    my $locus_end=    $column[5];
        if($count % 100000 == 0){
        	#print "$count: $_\n";
        	#print "$chr : $locus_chrom\n";
        }
    # handle allele="all" loci, which are are always ref, no-call, or PAR-called-in-X
    if ($allele eq  "all") {
	if ($type eq "ref" && $output_ref) {
	    # We could print as a ref entry, or not
	    print "$locus_chrom\t$locus_start\t$locus_end\tref\t150\n";
	} else {
	    # otherwise leave no entry for "all no-call" region
	    # print "$locus_chrom\t$locus_start\t$locus_end\tno-call\t0\n";
	}
	next;
    }
    
    # Not an allele="all" locus, therefore this is the first line, however it
    # still could be the only line in a ploidy=1 region (like chrM or most of chrY)
    
    my $no_call= ($type=~/^no-call/); # true if matches no-call-rc, no-call-ri or no-call
    my @strs= ("","","");  # collect each of the allele types
    my @seqs= ("","","");  # collect each of the allele sequences
    my $min_score= $score;
    $strs[$allele]=$type;
    $seqs[$allele]=$seq;
    
    # Read any additional lines for this locus
    while (getline()) {
	@column= split(/\t/); 		# split on tab into an array
	if ($column[0] ne $locus_id) {  # then this line is from the next locus
	    $unget= 1;
	    last;
	} else { # this is another row for the same locus	    
	    # in this code we aassume that the first line of a locus starts at the
	    # leftmost position for the whole locus, and that chrom does not change
	    # within a locus.  So we don't try to update either of these from
	    # lines after the first one.  We do have to update other items:
	    $allele=    $column[2];
	    $type=      $column[6];	
	    $seq=       $column[8];
	    $score=     $column[9];
	    my $end=    $column[5];	
	    $locus_end= $end   if ($end>$locus_end);   # take rightmost position as end
	    $min_score= $score if ($score<$min_score); # find lowest score
	    $no_call||= ($type=~/^no-call/); # match no-call-rc, no-call-ri or no-call
	    my $divider= ($strs[$allele] ? "+":""); # if this is a second row for this allele
	    $strs[$allele].= $divider.$type;
	    $seqs[$allele].= $seq;
	}
    }

    # Make a string allele1/allele2
    # Canonical form puts the ref on the 2nd allele for hets
    my $str; 
    my @color;
    my $bed_score= $output_call_types ? $min_score : 0;
    
    if (!$output_call_types) {
	$str= "";
    } elsif ($locus_ploidy==2) {
	if ($strs[1] eq "ref") {
	    $strs[1]= $strs[2]; $strs[2]="ref";
	    my $temp= $seqs[1];
	    $seqs[1]= $seqs[2]; $seqs[2]= $temp;
	}
	$str= "$strs[1]/$strs[2]";
	

    } else { # ploidy==1, no slash
	$str= $strs[1];
    }
    
    if($strs[1] eq "snp" || $strs[2] eq "snp"){
		$color[0] = 255;
	}else{
		$color[0] = 0;
	}
	if($strs[1] eq "sub" || $strs[2] eq "sub"){
		$color[1] = 255;
	}else{
		$color[1] = 0;
	}
	if($strs[1] eq "del" || $strs[2] eq "del"){
		$color[2] = 255;
	}else{
		$color[2] = 0;
	}
    if ($add_allele_seqs) {
	$str.=":";
	$seqs[1]="-" if ($seqs[1] eq ""); # these deal with deletions
	$seqs[2]="-" if ($seqs[2] eq "");
	if ($locus_ploidy==2) {
	    $str.= "$seqs[1]/$seqs[2]";
	} else {
	    $str.= "$seqs[1]";	    
	}
    }
    
    # print out as a locus unless no-called
    unless ($no_call) {
    	if($chr eq $locus_chrom || $chr eq "all"){
        print "$locus_chrom\t$locus_start\t$locus_end\t$str\t$bed_score\t+\t$locus_start\t$locus_end\t"; # print as a ref entry
        #concatenate to create colors
        #SNP = red
        #DEL = Blue
        #SUB = green
        #Other is a combination of the two
        print "$color[0],$color[1],$color[2]\n";
    	}
    } else {
	# otherwise leave no entry for "all no-call" region
        # print "$locus_chrom\t$locus_start\t$locus_end\tno-call\t0\n";	
    }
}

close FILE;
exit 0;





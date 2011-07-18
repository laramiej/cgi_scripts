#!/usr/bin/perl
use strict;

my $usage= "usage:  perl range.pl [option] tsv_file_name column_number start_value end_value\n".
           "option: -h	output CG header (header must exist in file if this is given)\n".
           "note:   column_number starts at 1\n";

# process -h optional arg which says we should output CG header
my $output_header= 0;
if ($ARGV[0] eq "-h") {
    $output_header= 1;
    shift @ARGV;
}

# check that we got 4 args
die $usage if (!defined($ARGV[0]) || !defined($ARGV[1]) || !defined($ARGV[2]) || !defined($ARGV[3]) || defined($ARGV[4]));

# arg 1 ($ARGV[0]) is the file name

# check that args 2-4 are all integers, and that arg 2 is >0
die "arg 2 '$ARGV[1]' must be an integer\n$usage" unless ($ARGV[1]=~/^\d+$/);
die "arg 3 '$ARGV[2]' must be an integer\n$usage" unless ($ARGV[2]=~/^\d+$/);
die "arg 4 '$ARGV[3]' must be an integer\n$usage" unless ($ARGV[3]=~/^\d+$/);
die "arg 2 '$ARGV[1]' must not be zero\n$usage"   unless ($ARGV[1]>0);

# set values from args
my $file=    $ARGV[0];
my $col_num= $ARGV[1]-1;  # decrement by one, as split output will be 0 indexed
my $start=   $ARGV[2];
my $end=     $ARGV[3];

# open file
if ($file=~/\.gz$/) { # if it ends with .gz
    open (FILE,"zcat $file |")  || die "can't open 'zcat $file'";
} elsif ($file=~/\.bz2$/) { # if it ends with .bz2
    open (FILE,"bzcat $file |") || die "can't open 'bzcat $file'";
} else {
    open (FILE,"<$file") || die "can't open '$file'";
}

# read until a line which starts with >
while (<FILE>) {
    print $_ if ($output_header);
    last if (/^\>/);
}

# parse and skip until we find the start
my $prev= "";	# this will hold the previous line
my $val;	# temp var for the position

while (<FILE>) {			# loop through lines in file
    chomp;	  			# remove training newline
    my @columns= split(/\t/); 		# split on tab into an array
    $val= $columns[$col_num];		# value in the specified column
    next unless ($val=~/^\d+$/);	# skip lines where val is not a number
    $prev= $_;				# save previous line
    last if ($val>=$start);		# stop loop at specified position
}

# output the last line
# also output previous line if we are already ahead of specified position
print "$prev\n" if ($prev && $val>$start);	
#print $_    if ($_);


while (<FILE>) { 			# now print lines until we find the end
	chomp;
	my @columns= split(/\t/); 	
	$val= $columns[$col_num];
	next unless ($val=~/^\d+$/);
	last if ($val>$end);		# this time stop after specified position
	print $_,"\n";
}

# close and quit
close FILE;
exit 0;





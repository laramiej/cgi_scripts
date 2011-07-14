Author: Jason Laramie jlaramie@completegenomics.com

findOverlappingJunctions.pl using the cgatools function junctionDiff will find overlapping structural variation junctions.
This is accomplished by reading in the regular junction file prior to running junctionDiff.  Any junctions that are outputted from 
junctionDiff are by definition different therefore removing these from the orginal SV junction list will find SV overlapping junctions.

Installation:
Three reference tracks need to be downloaded prior to running the script.  These are the repeatMask track, RefSeq track and the reference genome.  These can all be found on ftp.completegenomics.com.

Once downloaded, at the top of the script you will need to set full paths to these resources.  If you are only planning on 
using one genome build then the other genome build doesn't need to be set.
if($options{b} eq "hg18"){
	$path2repeatMask = "/path/to/referenceBuild/rmsk36.tsv.gz";
	$path2genes = "/path/to/referenceBuild/gene36.tsv.gz";
	$path2reference = "/path/to/referenceBuild/hg18.crr";
}elsif($options{b} eq "hg19"){
	$path2repeatMask = "/path/to/referenceBuild/rmsk37.tsv.gz";
	$path2genes = "/path/to/referenceBuild/gene37.tsv.gz";
	$path2reference = "/path/to/referenceBuild/build37.crr";

The path to cgatools (unless it is part of your environment (e.g. $PATH)) needs to be set at the top of the script.
my $path2cgatools = "/path/to/cgatoolsBinary/cgatools";

Additionally, a space delimited file (samples.txt example also downloaded) will need to be created that points at the genomes of interest.  The columns of this file need to be names "DNAid	normal	tumor".  The "tumor" column is optional if only step 2 is run.  If Step 1 is run (see below), the matched tumor genome location needs be given.

This script has the potential to run in two steps.  The first step will run junctionDiff in the way it was intended.  That is to compare, for instance, a tumor SV junctions to a normal SV junctions.  The junctions that are different between these two (e.g. in the tumor not in the normal) are then used in step 2 to find overlapping junctions.  Therefore this analysis answers the question, what junctions are not seen in the matched normal but are common across tumor sample?  


To run:
perl /home/jlaramie/scripts/findOverlappingJunctions.pl

perl findOverlappingJunctions.pl -s -d -h -f -t
	-s	A tab-delim sample manifest
	-d	full path of the output directory
	-b	genome build [hg18 or hg19]
	-h	flag to use high confidence calls
	-k	flag to skip step one which is to compare tumor to matched normal
	-t	(optional) for debugging
Died at /home/jlaramie/scripts/findOverlappingJunctions.pl line 245.

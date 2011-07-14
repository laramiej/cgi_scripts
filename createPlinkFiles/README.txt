Author: Jason Laramie 

createPlinkFiles.pl will create a set of plink files from a testvariants output containing only SNPs.  Good for cryptic relatedness or ethnicity (when combined with HapMap data) calculations or stats in PLINK.  PLINK can be found here:
	http://pngu.mgh.harvard.edu/~purcell/plink/

This assumes unrelated people.  Change the output'd tped file if this is not the case.

To run this:
	perl createPlink.pl -i -d
		-i	A test variant file
		-d	output directory
		-m	flag to include missing data (these will be encoded as NN even is it is a half-call)

Possibly a more complete way of converting this is to use our VCF converter located here:
	http://community.completegenomics.com/tools/m/cgtools/197.aspx
and use the VCF formatted file with plinks located here:
	http://atgu.mgh.harvard.edu/plinkseq/
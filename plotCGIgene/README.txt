Author: Jason M Laramie, jlaramie@completegenomics.com

plotCGIgene_cgatools1-3.pl will create a pretty plot of a single gene overlapping a lot of CGI data (see attached png (LIPI gene from NA19240); the real plot is slightly different).  Perl script needs to be run with plotCGIgene.R

The plots that are created are PDFs and require some UCSF tracks (In the future I will change this to refseq to match our annotation).  This script will only work with mastervar files generated from cgatools 1.3. I will update this when I get a chance to the new mastervars made from 1.4.

Compadibility: 
	Tested with Linux/Unix and Mac.  Not tested on windows and I am pretty sure it will not work due to the perl script needing to be able to execute the R script from the command line.  Finally, the current implementation is only compatible with cgatools version 1.3 and assembly version 1.11.  

Installation
	The following R libraries need to be installed before executing:
		gplots
	No R modules are required to be installed that would be different then the base installation.
	External data that is required to be downloaded before running the program:
		KnownGene (http://hgdownload.cse.ucsc.edu/downloads.html#human) - use whatever build you are on hg18 or hg19 and click the annotation link
		kgXref (http://hgdownload.cse.ucsc.edu/downloads.html#human)
		Known2locusId (http://hgdownload.cse.ucsc.edu/downloads.html#human)
		refLink (http://hgdownload.cse.ucsc.edu/downloads.html#human)
	A mastervar file will need to be created using cgatools 1.3 (http://cgatools.sourceforge.net/docs/1.3.0/) and placed into the ASM directory of its full genome directory.
	At the top of the plotCGIgene_cgatools1-3.pl script the path to specific files will need to added to these variables:
	my $pathToKnownGene = "/home/jlaramie/refs/$options{r}/knownGene.txt";
		$pathTokgXref
		$pathToKnown2locusId
		$pathTorefLink
		$pathToPlotter
		Finally, the file prefix for the way the mastervar is named needs to be here. If you change this to "var" it will clash with the CGI var file:
		$masterVarPrefix = "^mastervar"

How to run:
	perl /Users/jlaramie/Downloads/plotCGIgene/plotCGIgene_cgatools1-3.pl -i -g -d -n

	The following flags are available
		-i	/path/to/genome/data/ASM directory
		-g	input genelist
		-d	/path/to/output/directory
		-r	reference genome version hg18 or hg19
		-x	extra bp to add to start/stop of transcript (optional)
		-n	integer number basepairs. If any transcripts are bigger than this number the entire plot will not be created.(optional)

Written by Jason M Laramie jlaramie@completegenomics.com

Future development:
	Update code to cgatools 1.4 and assembly version 1.12
	Allow the user to choose which coverage metric is plotted
	Change transcript input to use Refseq to match our annotation and remove all of the excess UCSC tracks.
	
	



import sys
import bz2
import csv
import re

import pysam

references = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM')
referencelengths = (247249719L, 242951149L, 199501827L, 191273063L, 180857866L, 170899992L, 158821424L, 146274826L, 140273252L, 135374737L, 134452384L, 132349534L, 114142980L, 106368585L, 100338915L, 88827254L, 78774742L, 76117153L, 63811651L, 62435964L, 46944323L, 49691432L, 154913754L, 57772954L, 16571L)
cigar_pattern = re.compile('([0-9]*)([DMN])')

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class CigarParserError(Error):
    pass

def cigar_to_tuples(cigar_string):
    cigar_tuples = []
    for n, o in cigar_pattern.findall(cigar_string):
        if o == "M":
            cigar_tuples.append((0, int(n)))
        elif o == "N":
            cigar_tuples.append((3, int(n)))
        else:
            raise CigarParserError
    return cigar_tuples

#test cigar
#print cigar_to_tuples('10M5N10M1N13M')
#sys.exit()

f = bz2.BZ2File("ejdb_test.tsv.bz2", "r")

#f = bz2.BZ2File("evidenceJunctionDnbsBeta-GS000004357-ASM.tsv.bz2", "r")

csv_reader = csv.reader(f, delimiter='\t')

#dr_bam = pysam.Samfile("5137.bam", "wb", referencenames=references, referencelengths=referencelengths)
dr_bam = pysam.Samfile("5137.sam", "wh", referencenames=references, referencelengths=referencelengths)

for line in csv_reader:
    if line:
        if not line[0].startswith('#'):
            if line[0] == '5137':

                dr = pysam.AlignedRead()

                lseq = line[18][:35]
                rseq = line[18][35:]
                lscore = line[19][:35]
                rscore = line[19][35:]

                # Set qname 
                dr.qname = "".join([line[1], "-", line[2], "-", line[3], ":", line[4]])

                # Set sequence
                dr.seq = lseq

                # Set flag
                dr.flag = 97

                # Set rname
                dr.rname = dr_bam.gettid(line[7])

                # Set position
                dr.pos = int(line[8])

                # Set MAPQ
                dr.mapq = ord(line[10])

                # Set Cigar
                dr.cigar = cigar_to_tuples(line[9])

                dr.mrnm = dr_bam.gettid(line[13])
                dr.mpos= int(line[14])
                #dr.isize=int(line[17])
                dr.isize=0
                dr.qual=lscore
                dr.tags = ( ("R2", rseq), ("Q2", rscore), ("XS",1))
                dr_bam.write(dr)
                print repr(dr)

dr_bam.close()

"""
a.qname = "read_28833_29006_6945"
a.seq="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
a.flag = 99
a.rname = 0
a.pos = 32
a.mapq = 20
#a.cigar = ( (0,10), (2,1), (0,25) )
a.cigar = [(0, 10), (3, 5), (0, 10), (3, 1), (0, 13)]
a.mrnm = 0
a.mpos=199
a.isize=167
a.qual="<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"
a.tags = ( ("NM", 1),
           ("RG", "L1") )
"""

#!/usr/bin/env python
import sys
import bz2
import csv
import re
import os

import argparse

import pysam

from Bio.Seq import Seq


references = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM')
referencelengths = (247249719L, 242951149L, 199501827L, 191273063L, 180857866L, 170899992L, 158821424L, 146274826L, 140273252L, 135374737L, 134452384L, 132349534L, 114142980L, 106368585L, 100338915L, 88827254L, 78774742L, 76117153L, 63811651L, 62435964L, 46944323L, 49691432L, 154913754L, 57772954L, 16571L)

cigar_pattern = re.compile('([0-9]*)([MN])')
cg_cigar_pattern = re.compile('([0-9]*)([MNB])')

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class CigarParserError(Error):
    pass

def cigar_to_tuples(cigar_string):
    ct = []
    for n, o in cigar_pattern.findall(cigar_string):
        if o == "M":
            ct.append((0, int(n)))
        elif o == "N":
            ct.append((3, int(n)))
        else:
            raise CigarParserError, o
    return ct

def cg_cigar_to_tuples(cigar_string):
    ct = []
    for n, o in cg_cigar_pattern.findall(cigar_string):
        if o == "M":
            ct.append((0, int(n)))
        elif o == "N":
            ct.append((3, int(n)))
        elif o == "B":
            ct.append((9, int(n)))
        else:
            raise CigarParserError, o
    return ct

def cg_cigar_tuples_pack(ct):
    #print ct
    if ct[3][1] == 0:
        ct_n = [(0, (ct[0][1] - ct[1][1]) + ct[2][1] + ct[4][1]), ct[5], ct[6]] 
    else:
        ct_n = [(0, (ct[0][1] - ct[1][1]) + ct[2][1]), ct[3], ct[4], ct[5], ct[6]] 
    #print "CT_N: ", ct_n
    return ct_n

def proc_cg_seq(s, cs):
    # [:5] + [5+2:]
#    if cs[1][1] == 1:
#        print "s[ : %s ] + s[ %s + %s : ]" % (cs[0][1], cs[0][1], cs[1][1])
#        print rid
    new_s = s[ : cs[0][1] ] + s[ cs[0][1] + cs[1][1] : ]
    if rid == "19993257":
        print "s[ : %s ] + s[ %s + %s : ] len = %s" % (cs[0][1], cs[0][1], cs[1][1], len(new_s))
    return new_s
def proc_cg_qual(q, cs):
    return q[ : cs[0][1] ] + q[ cs[0][1] + cs[1][1]: ].strip()

class Junction: 
    def __init__(self, id):
        self.j_id = id
        self.a_pos = None
        self.a_dnb = None
        self.a_strand = None
        self.b_pos = None
        self.b_dnb = None
        self.b_strand = None

def convert_ej_to_sam(ej_file_handle, output_sam, junction_id = None):

    ej_csv_reader = csv.reader(ej_file_handle, delimiter='\t')

    for line in ej_csv_reader:
        if line:
            if not line[0].startswith('#'):
                if line[0] == junction_id:

                    dr_1 = pysam.AlignedRead()
                    dr_2 = pysam.AlignedRead()

                    #j = Junction(juncion_id)
                    #j.a_pos = 12020253
                    #j.a_dnb = 'R'
                    #j.a_strand = "-"
                    #j.b_pos = 26296211
                    #j.b_dnb = 'L'
                    #j.b_strand = "+"

                    s_1 = Seq(line[18][:35])
                    s_2 = Seq(line[18][35:])
                    #print s_1, s_2
                    cs_1 = line[9]
                    cs_2 = line[15]
                    qs_1 = line[19][:35]
                    qs_2 = line[19][35:]
                    strand_1 = line[6]
                    strand_2 = line[12]
                    dnb_1 = line[5]
                    dnb_2 = line[11]

                    # Set qname, same for both

                    global rid
                    rid = line[4]

                    dr_1.qname = "".join([line[1], "-", line[2], "-", line[3], ":", line[4]])
                    dr_2.qname = "".join([line[1], "-", line[2], "-", line[3], ":", line[4]])

                    if (dnb_1 == 'L' and strand_1 == "+") and (dnb_2 == 'R' and strand_2 == "+"):
                        continue
                    elif (dnb_1 == 'L' and strand_1 == "+") and (dnb_2 == 'R' and strand_2 == "-"):
                        continue
                    elif (dnb_1 == 'R' and strand_1 == "-") and (dnb_2 == 'L' and strand_2 == "+"):
                        # create the sequence first
                        # a is the RHS, so reverse strand
                        # take the second half and revcompl
                        #process the string according to the CGI cigar for DR 1

                        s_2_rev = s_2.reverse_complement()
                        cs_2_t = cg_cigar_to_tuples(cs_2)
                        # fix the sequence according to the 2b overlap
                        s_2_rev_cor = proc_cg_seq(s_2_rev, cs_2_t)

                        if rid == "19993257":
                            print len(s_2_rev_cor)

                        dr_1.seq = s_2_rev_cor.tostring()
                        dr_1.flag = 145
                        qs_2_cor = proc_cg_qual(qs_2, cs_2_t)
                        #print len(qs_2_cor)
                        dr_1.qual= qs_2_cor[::-1]

                        cs_1_t = cg_cigar_to_tuples(cs_1)
                        # Same for DR 2
                        s_1_cor = proc_cg_seq(s_1, cs_1_t)
                        if rid == "19993257":
                            print len(s_1_cor)
                        dr_2.seq = s_1_cor.tostring()
                        qs_1_cor = proc_cg_qual(qs_1, cs_1_t)
                        #print len(qs_1_cor)
                        dr_2.qual = qs_1_cor

                        dr_2.flag = 97

                        # Set rname
                        dr_1.rname = output_sam.gettid(line[13])
                        dr_2.rname = output_sam.gettid(line[13])

                        # Set position
                        dr_1.pos = int(line[8])
                        dr_2.pos = int(line[14])

                        # Set MAPQ, ASCII-33 ecoded

                        #print ord(line[10]), line[10]
                        #print ord(line[16]), line[16]

                        dr_1.mapq = ord(line[10])-33
                        dr_2.mapq = ord(line[16])-33
                        # Set Cigar
                        dr_1.cigar = cg_cigar_tuples_pack(cs_1_t)
                        dr_2.cigar = cg_cigar_tuples_pack(cs_2_t)

                        dr_1.mrnm = output_sam.gettid(line[13])
                        dr_2.mrnm = output_sam.gettid(line[13])

                        dr_1.mpos= int(line[14])
                        dr_2.mpos= int(line[8])

                        #dr.isize=int(line[17])
                        dr_1.isize=0
                        dr_2.isize=0

                        #dr.tags = ( ("R2", rseq), ("Q2", rscore), ("XS",1))
                        dr_1.tags = (("XS",1),)
                        dr_2.tags = (("XS",1),)
                    output_sam.write(dr_1)
                    output_sam.write(dr_2)

    output_sam.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert evidenceJunctionsDnb-GSXXX-ASM.tsv to SAM format.')

    #parser.add_argument('-i', dest='ej_filename', help='evidenceJunctionsDnb-GSXXX-ASM.tsv file ')
    #parser.add_argument('--input', dest='ej_filename', action='store_const', const=sum, default=max, help='')

    parser.add_argument('-j', dest='ej_junction_id', help='Junction ID to export')
    parser.add_argument('ej_filename', help='videnceJunctionsDnb-GSXXX-ASM.tsv file')
    parser.add_argument('ej_sam_filename', help='name of sam output file')

    args = parser.parse_args()

    ej_file_handle = bz2.BZ2File(args.ej_filename, "r")

    output_sam = pysam.Samfile(args.ej_sam_filename, "wh", referencenames=references, referencelengths=referencelengths)

    convert_ej_to_sam(ej_file_handle, output_sam, junction_id = args.ej_junction_id)

    #pysam.sort(args.ej_sam_filename, 'sorted_')


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

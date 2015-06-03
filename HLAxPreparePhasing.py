#!/usr/bin/env python
"""

Read input of format:
var.f.pos	ref	alt	length	f.pos   m.pos   c.pos
60		T	C	1	60      NA      NA
2305		G	A	1       2305    NA      NA
5181		AAAAAAA	AAAAAA	7    	5181    1541    1716
6817		T	C	1       6817    3170    3343

Read the matching VCF and Fasta for the data above.

For every entry in the primary format, find the variant in either the VCF file or the fasta.

Output final results for phasing.

Author: Rune Moellegaard Friborg <runef@birc.au.dk>
"""


import getopt
import sys


##############################################################
####################### Configuration ########################

VERSION="0.01"
UPDATED="2015-06-03"

##############################################################
####################### Main #################################

def main(args):
    pass

##############################################################
######################### Help ###############################

def usage():
    print("""HLAxPreparePhasing version """+VERSION+""" by Rune M. Friborg (updated """+UPDATED+""")
Usage:
  HLAxPreparePhasing --vcf=<file> --fa=<fasta file> --input=<input file>
""")

class ArgContainer():
    def __init__(self):
        self.vcffile    = ""
        self.fafile     = ""
        self.inputfile  = ""

if __name__ == '__main__':

    try:
        opts, dirs = getopt.getopt(sys.argv[1:], "", ["help", "vcf=", "fa=", "input="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    args = ArgContainer()
    for o, a in opts:
        if a and a[0] == "=":
            a = a[1:]

        if o == "--vcf":
            args.vcffile = a
        elif o == "--fa":
            args.fafile = a
        elif o == "--input":
            args.inputfile = a
        elif o == "--help":
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"

    if args.vcffile and args.fafile and args.inputfile:
        try:
            main(args)
        except KeyboardInterrupt:
            sys.exit(1)
    else:
        usage()
        sys.exit()

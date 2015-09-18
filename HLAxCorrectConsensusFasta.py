#!/usr/bin/env python
"""

Read fasta file and output a new fasta, where variants have been inserted based on the info in two columns

Author: Rune Moellegaard Friborg <runef@birc.au.dk>
"""

import getopt
import sys
import gzip
from Bio import SeqIO
import vcf
import itertools
import tempfile
import subprocess
import os


##############################################################
####################### Configuration ########################

VERSION="0.01"
UPDATED="2015-09-18"
PID=str(os.getpid())

##############################################################
####################### Classes #################################

class Fasta():
    """
    Reads a fasta files with a single sequence and stores it in self.seq
    """
    def __init__(self, filename):
        sys.stdout.write("Reading '"+filename+"'...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")    
        else:
            handle = open(filename, "r")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        if not len(records) == 1:
            sys.stderr.write("Error!: "+str(len(records))+" sequences read. Fasta file '"+ filename +"' must have exactly one sequence!\n")
            sys.exit(1)
        
        self.val = records[0].seq
        self.id = records[0].id
        sys.stdout.write("done\n")



##############################################################
####################### Main #################################

def main(args):
    # Read Fasta file
    faSequence = Fasta(args.fa)

    # Read variants and output new fasta
    
    
##############################################################
######################### Help ###############################

def usage():
    print("""HLAxCorrectConsensusFasta version """+VERSION+""" by Rune M. Friborg (updated """+UPDATED+""")
Usage:
  HLAxCorrectConsensusFasta \
    --fa-input=<fasta file> --var-pos=<file:column> --var-ref-len=<file:column>
    --var-new-seq=<file:column> --fa-output=<fasta file>
""")


class ArgContainer():
    def __init__(self):
        self.fa_input     = ""
        self.fa_output    = ""
        self.var_pos_file = ""
        self.var_pos_col  = None
        self.var_ref_len_file = ""
        self.var_ref_len_col  = None
        self.var_new_seq_file = ""
        self.var_new_seq_col  = ""

    def ok(self):
        err = 0
        if not self.fa_input:
            sys.stderr.write("Missing argument: --fa-input\n")
            err = 1
        if not self.var_pos_file:
            sys.stderr.write("Missing file part argument: --var-pos\n")
            err = 1
        if not self.var_pos_col:
            sys.stderr.write("Missing column part argument: --var-pos\n")
            err = 1
        if not self.var_ref_len_file:
            sys.stderr.write("Missing file part argument: --var-ref-len\n")
            err = 1
        if not self.var_ref_len_col:
            sys.stderr.write("Missing column part argument: --var-ref-len\n")
            err = 1
        if not self.var_new_seq_file:
            sys.stderr.write("Missing file part argument: --var-new-seq\n")
            err = 1
        if not self.var_new_seq_col:
            sys.stderr.write("Missing column part argument: --var-new-seq\n")
            err = 1
        if not self.fa_output:
            sys.stderr.write("Missing argument: --fa-output\n")
            err = 1
        if err:
            sys.stderr.write("\n")

        return not err

        

if __name__ == '__main__':

    try:
        opts, dirs = getopt.getopt(sys.argv[1:], "", ["help", "fa-input=", "var-pos=", "var-ref-len=", "var-new-seq=", "fa-output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    args = ArgContainer()
    for o, a in opts:
        if a and a[0] == "=":
            a = a[1:]

     
        if o == "--fa-input":
            args.fa_input = a
        elif o == "--fa-output":
            args.fa_output = a
        elif o == "--var-pos":
            try:
                args.var_pos_file, args.var_pos_col = a.split(":")
            except ValueError:
                pass
        elif o == "--var-ref-len":
            try:
                args.var_ref_len_file, args.var_ref_len_col = a.split(":")
            except ValueError:
                pass
        elif o == "--var-new-seq":
            try:
                args.var_new_seq_file, args.var_new_seq_col = a.split(":")
            except ValueError:
                pass
        elif o == "--help":
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"

    if args.ok():
        try:
            main(args)
        except KeyboardInterrupt:
            sys.exit(1)
    else:
        usage()
        sys.exit()

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

VERSION="0.02"
UPDATED="2015-09-22"
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
        
        self.val = str(records[0].seq)
        self.id = records[0].id
        sys.stdout.write("done\n")


class Custom():
    """
    Reads a column in a file and stores the list
    """
    def __init__(self, filename, colnum, sep=' ', has_header=True):
        self.filename = filename
        self.colnum = colnum
        self.data = []

        sys.stdout.write("Parsing '"+filename+"':"+str(colnum)+"...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")

        if has_header:
            handle.next()

        # Parse lines
        i = 0
        for line in handle:
            try:
                i += 1

                if line == "\n":
                    continue

                l = line.strip().split(' ')
                self.data.append(l[colnum-1])
                
            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1)
        sys.stdout.write("done\n")
        handle.close()

            


##############################################################
####################### Main #################################

def main(args):
    # Read Fasta file
    faSequence = Fasta(args.fa_input)

    # Read variants and output new fasta
    posList = Custom(args.var_pos_file, args.var_pos_col).data
    refLenList = Custom(args.var_ref_len_file, args.var_ref_len_col).data
    newSeqList = Custom(args.var_new_seq_file, args.var_new_seq_col).data

    if (len(posList) == len(refLenList) and len(posList) == len(newSeqList)):
        sys.stdout.write("Writing '"+args.fa_output+"...")
        ohandle_fasta = open(args.fa_output, "w")
        ohandle_fasta.write(">"+str(faSequence.id)+"\n")

        region_start = 0
        region_end = 0
        prev_entry = (0, 0, "")
        for entry in zip(posList, refLenList, newSeqList):
            pos = int(entry[0])
            refLen = int(entry[1])
            newSeq = entry[2]

            if pos < region_start:
                sys.stderr.write("WARNING! " + str(entry) + " conflicts with previous variant " + str(prev_entry) +"\n")
                continue

            region_end = pos

            # Write sequence between variants
            ohandle_fasta.write(faSequence.val[region_start:region_end])

            # Write new seq
            ohandle_fasta.write(newSeq)
            
            # Set region start and skip ref seq
            region_start = region_end+refLen

            prev_entry = entry

        # Write final region
        ohandle_fasta.write(faSequence.val[region_start:])        
        ohandle_fasta.close()

        sys.stdout.write("done\n")

    else:
        sys.stderr.write("Error! Mismatch in length of columns")
        sys.stderr.write("pos: "+str(len(posList))+"\n")
        sys.stderr.write("seq len: "+str(len(refLenList))+"\n")
        sys.stderr.write("new seq: "+str(len(newSeqList))+"\n")

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
                args.var_pos_col = int(args.var_pos_col)
            except ValueError:
                pass
        elif o == "--var-ref-len":
            try:
                args.var_ref_len_file, args.var_ref_len_col = a.split(":")
                args.var_ref_len_col = int(args.var_ref_len_col)
            except ValueError:
                pass
        elif o == "--var-new-seq":
            try:
                args.var_new_seq_file, args.var_new_seq_col = a.split(":")
                args.var_new_seq_col = int(args.var_new_seq_col)
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
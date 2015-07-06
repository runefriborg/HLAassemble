#!/usr/bin/env python
"""


Read fasta input from parent and child.

Read the phased information and filter out any phased entries with positions crossing

Cut the fasta input up, by using the phased information, which provides a position for the child and the parent.

For every section, run global alignment and output the results

Author: Rune Moellegaard Friborg <runef@birc.au.dk>
"""

VCF_OUTPUT_TEMPLATE=""

import getopt
import sys
import gzip
from Bio import SeqIO
import vcf
import itertools
import nwalign as nw


##############################################################
####################### Configuration ########################

VERSION="0.01"
UPDATED="2015-07-05"

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
        sys.stdout.write("done\n")


class VCF():
    """
    Reads a VCF file and stores all records by position in self.dict
    """
    def __init__(self, filename):
        sys.stdout.write("Reading '"+filename+"'...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")
        vcf_reader = vcf.Reader(handle)
        vcfDict = {}
        for record in vcf_reader:
            vcfDict[record.POS] = record
        handle.close()

        self.val = vcfDict
        sys.stdout.write("done\n")


class PhasedPositions():
    """
    Read the input of format:
    #CHILDFATHERMOTHER
    3343   6817    3170

    And provide filtering methods
    """
    def __init__(self, filename, parent):

        self.content = []
        
        sys.stdout.write("Parsing '"+filename+"'...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")

        # Skip header
        handle.next()
        
        # Parse lines
        i = 0
        for line in handle:
            try:
                i += 1

                if line == "\n":
                    continue

                l = line.strip().split('\t')
                c_pos, f_pos, m_pos = l
                
                # Typecast
                c_pos = int(c_pos)

                if parent == 'f':
                    if f_pos != 'None':
                        f_pos = int(f_pos)
                        self.content.append((c_pos, f_pos))
                else:
                    if m_pos != 'None':
                        m_pos = int(m_pos)
                        self.content.append((c_pos, m_pos))
                
            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1)            

        sys.stdout.write("done\n")
        handle.close()
    
    def clean(self):
        """
        Removes crossed phased positions 
        """
        
        #L = []
        
        
        #for c_pos, m_pos, f_pos in self.content:
        pass

        
##############################################################
####################### Main #################################

def main(args):

    # Read Fasta files (child and parent)
    faSequence = []
    for filename in [args.c_fa, args.parent_fa]:
        faSequence.append(Fasta(filename))

    phasePos = PhasedPositions(args.phase_pos, args.parent)

##############################################################
######################### Help ###############################

def usage():
    print("""HLAxGlobalAlignParentAndChild version """+VERSION+""" by Rune M. Friborg (updated """+UPDATED+""")
Usage:
  HLAxAssembleFastaFromPhasing \
    --parent-vcf=<file> --parent-fa=<fasta file from previous HLAx stage> --parent=m|f \
    --c-vcf=<file> --c-fa=<fasta file from previous HLAx stage> \
    --phase-pos=<output from previous HLAx stage> --c-output-prefix=<file prefix>
""")


class ArgContainer():
    def __init__(self):
        self.parent_vcf    = ""
        self.parent_fa     = ""
        self.parent        = None
        self.c_vcf    = ""
        self.c_fa     = ""
        self.phase_pos  = ""
        self.c_output_prefix  = ""

    def ok(self):
        err = 0
        if not self.parent_vcf:
            sys.stderr.write("Missing argument: --parent-vcf\n")
            err = 1
        if not self.parent_fa:
            sys.stderr.write("Missing argument: --parent-fa\n")
            err = 1
        if not self.parent:
            sys.stderr.write("Missing argument: --parent\n")
            err = 1
        if not self.c_vcf:
            sys.stderr.write("Missing argument: --c-vcf\n")
            err = 1
        if not self.c_fa:
            sys.stderr.write("Missing argument: --c-fa\n")
            err = 1
        if not self.phase_pos:
            sys.stderr.write("Missing argument: --c-phase-pos\n")
            err = 1
        if not self.c_output_prefix:
            sys.stderr.write("Missing argument: --c-output-prefix\n")
            err = 1

        if err:
            sys.stderr.write("\n")

        return not err

        

if __name__ == '__main__':

    try:
        opts, dirs = getopt.getopt(sys.argv[1:], "", ["help", "parent-vcf=", "parent-fa=", "parent=", "c-vcf=", "c-fa=", "phase-pos=", "c-output-prefix="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    args = ArgContainer()
    for o, a in opts:
        if a and a[0] == "=":
            a = a[1:]

     
        if o == "--parent-vcf":
            args.parent_vcf = a
        elif o == "--parent-fa":
            args.parent_fa = a
        elif o == "--parent":
            args.parent = a
        elif o == "--c-vcf":
            args.c_vcf = a
        elif o == "--c-fa":
            args.c_fa = a
        elif o == "--phase-pos":
            args.phase_pos = a
        elif o == "--c-output-prefix":
            args.c_output_prefix = a
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

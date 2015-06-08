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
import gzip
from Bio import SeqIO
import vcf

##############################################################
####################### Configuration ########################

VERSION="0.02"
UPDATED="2015-06-08"
VCF_OUTPUT_HEADER="CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFATHER\tMOTHER\tCHILD"
VCF_OUTPUT_TEMPLATE="..."

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


class JMJProcess():
    """
    Descripe the process.
    """
    def __init__(self, vcf, seq):
        self.child_vcf = vcf.val
        self.child_seq = seq.val
        self.custom_variants = []

    def preprocess(self, filename):
        """
        For every variant in the mother and father identify the variants in the child
        
        Use this method for the father and mother
        """
        new_variants= {}

        sys.stdout.write("Preprocessing '"+filename+"'...")
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
                pos, ref, alt, length, f_pos, m_pos, c_pos = l
            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1)            

            #record.alleles.append('K')
            #print(record.alleles)
            
            #sys.exit(0)
            # Skip if NA
            if (c_pos == 'NA'):
                continue

            # Fetch VCF record
            #record = vcf[int(pos)]

            # Type cast to int
            pos = int(pos)
            c_pos = int(c_pos)

            # Extract variants
            v = [ref, alt]

            if self.child_vcf.has_key(c_pos):
                # Child has a variant on the same position in it's VCF                
                child_record = self.child_vcf[c_pos]

                cv = [child_record.REF]
                cv.extend([x.sequence for x in child_record.ALT])

                if new_variants.has_key(c_pos):
                    # Variant conflict, Ignore variant!
                    new_variants[c_pos] = (None, None)
                else:
                    new_variants[c_pos] = (pos, cv)
            else:
                # Search fasta

                # Correct for one-indexing TODO!!! WARNING. Must be changed when vcf files have been corrected
                # fasta_c_pos = c_pos-1
                fasta_c_pos = c_pos

                # Get variants
                cv = []
                for item in v:
                    cv.append(str(self.child_seq[fasta_c_pos:(fasta_c_pos+len(item))]))

                if new_variants.has_key(c_pos):
                    # Variant conflict, Ignore variant!
                    new_variants[c_pos] = (None, None)
                else:
                    new_variants[c_pos] = (pos, cv)
        
        # Add newly found custom variants
        self.custom_variants.append((filename, new_variants))
        
        sys.stdout.write("done\n")
        handle.close()
        
    def process(self, filename, ofilename, parent_vcf_obj, parent_seq_obj):

        parent_vcf = [x.val for x in parent_vcf_obj]
        parent_seq = [x.val for x in parent_seq_obj]

        sys.stdout.write("Processing '"+filename+"'...")
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")
        else:
            handle = open(filename, "r")

        # Open output file
        if ofilename[-2:] == 'gz':
            ohandle = gzip.open(ofilename, "w")
        else:
            ohandle = open(ofilename, "w")

        # Skip header
        handle.next()

        
        # Parse lines
        i = 0
        data = {}
        for line in handle:
            try:
                i += 1

                if line == "\n":
                    continue

                l = line.strip().split('\t')
                pos, ref, alt, length, f_pos, m_pos, c_pos = l
                
                data[int(pos)] = (ref, alt, length, f_pos, m_pos)
            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1) 

        # Find max position
        max_pos = 0
        for seq in parent_seq + [self.child_seq]:
            if len(seq) > max_pos:
                max_pos = len(seq)

        # Do a full sweep of all positions and construct the resulting record
        for i in xrange(0, max_pos):
            ohandle.write(str(i) + "\n")
            
            pass
            

            #other_pos = getChildVariantPosBefore(pos)
            #while (other_pos != None):
            #    print(other_pos)
            #    other_pos = getChildVariantPosBefore(pos)        
        
            #record = childVCFRecords[int(pos)]

            #record.ALT.append(vcf.model._Substitution("A"))
            #print(record.ALT)

            #vcf_writer.write_record(record)

        sys.stdout.write("done\n")
        ohandle.close()
        handle.close()
    

        
        
##############################################################
####################### Main #################################

def main(args):

    # Read vcf file
    vcfRecords = []
    for filename in [args.f_vcf, args.m_vcf, args.c_vcf]:
        vcfRecords.append(VCF(filename))

    # Read Fasta files
    faSequence = []
    for filename in [args.f_fa, args.m_fa, args.c_fa]:
        faSequence.append(Fasta(filename))

    # The main object for storing and computing the phasing values.
    jmj = JMJProcess(vcfRecords[2], faSequence[2])

    # Preprocessing the data for father and mother
    for filename in [args.f_input, args.m_input]:
        jmj.preprocess(filename)
    
    #Processing the child
    jmj.process(args.c_input, args.c_output, vcfRecords[0:1], faSequence[0:1])
    



##############################################################
######################### Help ###############################

def usage():
    print("""HLAxPreparePhasing version """+VERSION+""" by Rune M. Friborg (updated """+UPDATED+""")
Usage:
  HLAxPreparePhasing --f-vcf=<file> --f-fa=<fasta file> --f-input=<input file> \
    --m-vcf=<file> --m-fa=<fasta file> --m-input=<input file> \
    --c-vcf=<file> --c-fa=<fasta file> --c-input=<input file> --c-output=<output file>
""")


class ArgContainer():
    def __init__(self):
        self.f_vcf    = ""
        self.f_fa     = ""
        self.f_input  = ""
        self.m_vcf    = ""
        self.m_fa     = ""
        self.m_input  = ""
        self.c_vcf    = ""
        self.c_fa     = ""
        self.c_input  = ""
        self.c_output  = ""

    def ok(self):
        err = 0
        if not self.f_vcf:
            sys.stderr.write("Missing argument: --f-vcf\n")
            err = 1
        if not self.f_fa:
            sys.stderr.write("Missing argument: --f-fa\n")
            err = 1
        if not self.f_input:
            sys.stderr.write("Missing argument: --f-input\n")
            err = 1
        if not self.m_vcf:
            sys.stderr.write("Missing argument: --m-vcf\n")
            err = 1
        if not self.m_fa:
            sys.stderr.write("Missing argument: --m-fa\n")
            err = 1
        if not self.m_input:
            sys.stderr.write("Missing argument: --m-input\n")
            err = 1
        if not self.c_vcf:
            sys.stderr.write("Missing argument: --c-vcf\n")
            err = 1
        if not self.c_fa:
            sys.stderr.write("Missing argument: --c-fa\n")
            err = 1
        if not self.c_input:
            sys.stderr.write("Missing argument: --c-input\n")
            err = 1
        if not self.c_output:
            sys.stderr.write("Missing argument: --c-output\n")
            err = 1

        if err:
            sys.stderr.write("\n")

        return not err

        

if __name__ == '__main__':

    try:
        opts, dirs = getopt.getopt(sys.argv[1:], "", ["help", "f-vcf=", "f-fa=", "f-input=", "m-vcf=", "m-fa=", "m-input=", "c-vcf=", "c-fa=", "c-input=", "c-output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        sys.exit(2)
    
    args = ArgContainer()
    for o, a in opts:
        if a and a[0] == "=":
            a = a[1:]

     
        if o == "--f-vcf":
            args.f_vcf = a
        elif o == "--f-fa":
            args.f_fa = a
        elif o == "--f-input":
            args.f_input = a
        elif o == "--m-vcf":
            args.m_vcf = a
        elif o == "--m-fa":
            args.m_fa = a
        elif o == "--m-input":
            args.m_input = a
        elif o == "--c-vcf":
            args.c_vcf = a
        elif o == "--c-fa":
            args.c_fa = a
        elif o == "--c-input":
            args.c_input = a
        elif o == "--c-output":
            args.c_output = a
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

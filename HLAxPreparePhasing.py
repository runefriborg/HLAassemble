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

VERSION="0.01"
UPDATED="2015-06-03"
VCF_OUTPUT_TEMPLATE="template.vcf"

##############################################################
####################### Main #################################

def main(args):

    # Read Fasta files
    faSequence = []
    for filename in [args.f_fa, args.m_fa, args.c_fa]:
        if filename[-2:] == 'gz':
            handle = gzip.open(filename, "r")    
        else:
            handle = open(filename, "r")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        if not len(records) == 1:
            sys.stderr.write("Error!: "+str(len(records))+" sequences read. Fasta file '"+ filename +"' must have exactly one sequence!\n")
            sys.exit(1)
        
        faSequence.append(records[0].seq)


    # Read vcf file
    vcfRecords = []
    for filename in [args.f_vcf, args.m_vcf, args.c_vcf]:
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
        vcfRecords.append(vcfDict)
        sys.stdout.write("done\n")

    

    # Read input file
    individual = 0                # 0 = father, 1 = mother, 2 = child
    childVCFRecords = vcfRecords[2]
    childFaSequence = faSequence[2]
    Variants = ({},{}) # indexed using the positions of the child
    childVariants = {}

    # Preprocessing the data for father and mother
    for filename in [args.f_input, args.m_input]:
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
            record = vcfRecords[individual][int(pos)]

            # Type cast to int
            pos = int(pos)
            c_pos = int(c_pos)

            # Add variant
            v = [record.REF]
            v.extend([x.sequence for x in record.ALT])
            Variants[individual][c_pos] = (pos, v)

            if childVCFRecords.has_key(c_pos):
                # Child has a variant on the same position in it's VCF                
                child_record = childVCFRecords[c_pos]

                cv = [child_record.REF]
                cv.extend([x.sequence for x in child_record.ALT])

                if childVariants.has_key(c_pos):
                    childVariants[c_pos][1].extend(cv)
                else:
                    childVariants[c_pos] = (c_pos, cv)

            else:
                # Search fasta

                # Correct for one-indexing
                fasta_c_pos = c_pos-1

                # Get variants
                cv = []
                for item in v:
                    cv.append(str(childFaSequence[fasta_c_pos:(fasta_c_pos+len(item))]))

                if childVariants.has_key(c_pos):
                    childVariants[c_pos][1].extend(cv)
                else:
                    childVariants[c_pos] = (c_pos, cv)
            
        
        sys.stdout.write("done\n")
        handle.close()
        individual += 1
        
    #Processing the child
    filename = args.c_input
    sys.stdout.write("Processing '"+filename+"'...")
    if filename[-2:] == 'gz':
        handle = gzip.open(filename, "r")
    else:
        handle = open(filename, "r")

    # Open output file
    ofilename = args.c_output
    if ofilename[-2:] == 'gz':
        ohandle = gzip.open(ofilename, "w")
    else:
        ohandle = open(ofilename, "w")
    vcf_template = vcf.Reader(filename=VCF_OUTPUT_TEMPLATE)
    vcf_writer = vcf.Writer(ohandle, vcf_template)


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

        record = childVCFRecords[int(pos)]

        #record.ALT.append(vcf.model._Substitution("A"))
        #print(record.ALT)

        vcf_writer.write_record(record)

    sys.stdout.write("done\n")
    ohandle.close()
    handle.close()


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

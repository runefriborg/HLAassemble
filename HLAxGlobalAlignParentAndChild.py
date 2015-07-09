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

VERSION="0.03"
UPDATED="2015-07-09"

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
                c0_pos, c1_pos, f0_pos, _ , m0_pos, _ = l
                
                # Typecast

                if parent == 'f':
                    if f0_pos != 'None':
                        p_pos = int(f0_pos)
                        c_pos = int(c0_pos)
                        self.content.append((c_pos, p_pos))
                else:
                    if m_pos != 'None':
                        m_pos = int(m_pos)
                        p_pos = int(m0_pos)
                        c_pos = int(c1_pos)
                        self.content.append((c_pos, p_pos))

                
            except ValueError:
                sys.stderr.write("Error!: Failing parsing line "+str(i)+" in '"+filename+"'!\nGot: "+str(l)+"\n")
                sys.exit(1)            

        sys.stdout.write("done\n")
        handle.close()
    
    def clean(self):
        """
        Removes crossed phased positions 
        """

    
        L = []
        
        last = 0
        for c_pos, p_pos in self.content:
            if p_pos > last:
                L.append((c_pos, p_pos))
            else:
                L.pop()
                print(str(last)+ " > " +str(p_pos))

            last = L[-1][1]

            
        print(len(self.content))
        print(len(L))

        self.content = L

class Alignment():
    def __init__(self, fasta, positions, parent, max_mem):
        
        self.parent = parent
        self.max_mem_mb = max_mem*1024
        self.fasta_id = (str(fasta[0].id), str(fasta[1].id))
        self.fasta = (str(fasta[0].val), str(fasta[1].val))
        self.segments = positions

    def process(self, output_prefix):
        
        sys.stdout.write("Main processing:\n")


        sys.stdout.write("Running global alignment on "+str(len(self.segments)+2)+" segments\n")
        start = self.segments[0]

        ohandle_fasta_child = open(output_prefix + '.consensus_with_'+self.parent+'.c.fa', "w")
        ohandle_fasta_parent = open(output_prefix + '.consensus_with_c.'+self.parent+'.fa', "w")
        sys.stdout.write("\tWriting " + output_prefix + '.consensus_with_'+self.parent+".c.fa...\n")
        sys.stdout.write("\tWriting " + output_prefix + '.consensus_with_c.'+self.parent+".fa...\n")


        # Add header
        ohandle_fasta_child.write(">"+self.fasta_id[0]+"\n")
        ohandle_fasta_parent.write(">"+self.fasta_id[1]+"\n")

        # Add fasta before first segment
        ohandle_fasta_child.write(self.fasta[0][:start[0]])
        ohandle_fasta_parent.write(self.fasta[1][:start[1]])


        for next_phased_pos in self.segments[1:]:
            child = self.fasta[0][start[0]:next_phased_pos[0]]
            parent = self.fasta[1][start[1]:next_phased_pos[1]]
            
            sys.stdout.write(str(start[0]))
            align_child, align_parent, score = self.nwalign(child, parent)
            print(" : Final score: " + str(score))
            
            for seq in align_child:
                ohandle_fasta_child.write(seq)

            for seq in align_parent:
                ohandle_fasta_parent.write(seq)
                
            start = next_phased_pos
        
        # Add fasta after last segment
        ohandle_fasta_child.write(self.fasta[0][start[0]:])
        ohandle_fasta_parent.write(self.fasta[1][start[1]:])
       
        ohandle_fasta_child.close()
        ohandle_fasta_parent.close()

        sys.stdout.write("Main processing completed.\n")        

    def nwalign(self, child, parent, l=0):

        result_child = []
        result_parent = []
        score = 0

        estimated_mem = ((15*len(child)*len(parent))/1024)/1024

        sys.stdout.write("..+"+str(len(child))+"|"+str(len(parent))+" mem: " +str(estimated_mem)+"MB")
        sys.stdout.flush()


        if estimated_mem > self.max_mem_mb:

            # First split on child
            # Splits both sequences based on the half lenght of the child sequence
            # Aligns from the start
            split = len(child)/2
            sys.stdout.write("\n\t"+"\t"*l)
            sub_result_child_1, sub_result_parent_1, score_1  = self.nwalign(child[:split], parent[:split], l=l+1)
            sys.stdout.write("\n\t"+"\t"*l)
            sub_result_child_2, sub_result_parent_2, score_2  = self.nwalign(child[split:], parent[split:], l=l+1)
            sub_score = score_1+score_2

            # Set initial best
            best_score = sub_score
            best_result_child = sub_result_child_1 + sub_result_child_2
            best_result_parent = sub_result_parent_1 + sub_result_parent_2

            # Second split on child
            # Splits both sequences based on the half lenght of the child sequence
            # Aligns from the end
            split = -1*(len(child)/2)
            sys.stdout.write("\n\t"+"\t"*l)
            sub_result_child_1, sub_result_parent_1, score_1  = self.nwalign(child[:split], parent[:split], l=l+1)
            sys.stdout.write("\n\t"+"\t"*l)
            sub_result_child_2, sub_result_parent_2, score_2  = self.nwalign(child[split:], parent[split:], l=l+1)
            sub_score = score_1+score_2

            # Update best
            if sub_score > best_score:
                best_score = sub_score
                best_result_child = sub_result_child_1 + sub_result_child_2
                best_result_parent = sub_result_parent_1 + sub_result_parent_2

            # First split on parent
            # Splits both sequences based on the half lenght of the parent sequence
            # Aligns from the start
            split = len(parent)/2
            sys.stdout.write("\n\t"+"\t"*l)
            sub_result_child_1, sub_result_parent_1, score_1  = self.nwalign(child[:split], parent[:split], l=l+1)
            sys.stdout.write("\n\t"+"\t"*l)
            sub_result_child_2, sub_result_parent_2, score_2  = self.nwalign(child[split:], parent[split:], l=l+1)
            sub_score = score_1+score_2

            # Update best
            if sub_score > best_score:
                best_score = sub_score
                best_result_child = sub_result_child_1 + sub_result_child_2
                best_result_parent = sub_result_parent_1 + sub_result_parent_2

            # Second split on parent
            # Splits both sequences based on the half lenght of the parent sequence
            # Aligns from the end
            split = -1*(len(parent)/2)
            sys.stdout.write("\n\t"+"\t"*l)
            sub_result_child_1, sub_result_parent_1, score_1  = self.nwalign(child[:split], parent[:split], l=l+1)
            sys.stdout.write("\n\t"+"\t"*l)
            sub_result_child_2, sub_result_parent_2, score_2  = self.nwalign(child[split:], parent[split:], l=l+1)
            sub_score = score_1+score_2

            # Update best
            if sub_score > best_score:
                best_score = sub_score
                best_result_child = sub_result_child_1 + sub_result_child_2
                best_result_parent = sub_result_parent_1 + sub_result_parent_2


            score = best_score            
            result_child.extend(best_result_child)
            result_parent.extend(best_result_parent)
        else:                    
        
            align_child, align_parent = nw.global_align(child, parent, gap_open=-6, gap_extend=-2, matrix='NUC.4.4')
        
            score = nw.score_alignment(align_child, align_parent, gap_open=-6, gap_extend=-2, matrix='NUC.4.4')
            
            result_child.extend(align_child)
            result_parent.extend(align_parent)

        return result_child, result_parent, score

##############################################################
####################### Main #################################

def main(args):

    # Read Fasta files (child and parent)
    faSequence = []
    for filename in [args.c_fa, args.parent_fa]:
        faSequence.append(Fasta(filename))

    phasePos = PhasedPositions(args.phase_pos, args.parent)
    print "CLean1"
    phasePos.clean()
    print "CLean2"
    phasePos.clean()

    align = Alignment(faSequence, phasePos.content, args.parent, max_mem=args.max_mem)
    align.process(output_prefix=args.c_output_prefix)
    
    
##############################################################
######################### Help ###############################

def usage():
    print("""HLAxGlobalAlignParentAndChild version """+VERSION+""" by Rune M. Friborg (updated """+UPDATED+""")
Usage:
  HLAxAssembleFastaFromPhasing \
    --parent-vcf=<file> --parent-fa=<fasta file from previous HLAx stage> --parent=m|f \
    --c-vcf=<file> --c-fa=<fasta file from previous HLAx stage> \
    --phase-pos=<output from previous HLAx stage> --c-output-prefix=<file prefix> \
    --max-mem=1024 (gigabyte)
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
        self.max_mem  = 1024

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
        opts, dirs = getopt.getopt(sys.argv[1:], "", ["help", "parent-vcf=", "parent-fa=", "parent=", "c-vcf=", "c-fa=", "phase-pos=", "c-output-prefix=", "max-mem="])
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
        elif o == "--max-mem":
            args.max_mem = int(a)
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

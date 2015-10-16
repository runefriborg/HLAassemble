#!/bin/bash

./HLAxCorrectConsensusFasta.py --fa-input 1009.haplotype.cf.consensus.formatted.fa  --var-pos=1009.cf.call.poo.genotype.v3.txt:2 --var-ref-len=1009.cf.call.poo.genotype.v3.txt:4 --var-new-seq=1009.cf.call.poo.genotype.v3.txt:3 --var-new-alt=1009.cf.call.poo.genotype.v3.txt:13 --fa-output=1009.haplotype.cf.consensus.formatted.corrected.fa --vcf-output=1009.haplotype.cf.consensus.formatted.corrected.vcf

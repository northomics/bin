#!/usr/bin/env python3

'''
    this is an in silico trypsin digestion program. The input is a fasta file which contains protein sequence to be digested, the output is a txt file which contains all trypsin digested peptides and corresponding protein accessions.
'''

import sys
import os
#import getopt
import argparse
from collections import defaultdict
from Bio import SeqIO
from tqdm import tqdm

def TRYPSIN(proseq,miss_cleavage):
    peptides=[]
    #peptides = defaultfict(list)
    cut_sites=[0]
    for i in range(0,len(proseq)-1):
        if proseq[i]=='K' and proseq[i+1]!='P':
            cut_sites.append(i+1)
        elif proseq[i]=='R' and proseq[i+1]!='P':
            cut_sites.append(i+1)
    
    if cut_sites[-1]!=len(proseq):
        cut_sites.append(len(proseq))

    if len(cut_sites)>2:
        if  miss_cleavage==0:
            for j in range(0,len(cut_sites)-1):
                peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])

        elif miss_cleavage==1:
            for j in range(0,len(cut_sites)-2):
                peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j+2]])
            
            peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])

        elif miss_cleavage==2:
            for j in range(0,len(cut_sites)-3):
                peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j+2]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j+3]])
            
            peptides.append(proseq[cut_sites[-3]:cut_sites[-2]])
            peptides.append(proseq[cut_sites[-3]:cut_sites[-1]])
            peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])
    else: #there is no trypsin site in the protein sequence
        peptides.append(proseq)
    return peptides


################  Comand-line arguments ################
#if len(sys.argv[1:])<=1:  ### Indicates that there are insufficient number of command-line arguments
#    print("Warning! wrong command, please read the manual in Readme.tx.")
#    print "Example: python trypsin.py --input input_filename --output output_filename --miss 1"
#else:
#    options, remainder = getopt.getopt(sys.argv[1:],'', ['input=',
#                                                         'miss=',
#                                                         'output='])
#    for opt, arg in options:
#        if opt == '--input': input_file=arg
#        elif opt == '--miss': n=int(arg)  #number of miss cleavage allowed
#        elif opt == '--output':output_file=arg
#        else:
#            print "Warning! Command-line argument: %s not recognized. Exiting..." % opt; sys.exit()
#
usage = "usage: %prog [options] -i input.fasta -c 2 -m 5 -a 50 -o output.fasta"

parser = argparse.ArgumentParser(description="in siloco trypsin digest")
parser.add_argument("-i", "--input", dest="input", required=True, metavar="input.fasta", help="fasta file of proteins to digest")
parser.add_argument("-o", "--output", dest="output", required=True, metavar="output.fasta", help="fasta file of peptide sequences")
parser.add_argument("-c", "--cleavages", dest="c", required=False, default = 2, type = int,  metavar="missed cleavages", help="missed trypsin cleavages")
parser.add_argument("-m", "--minsize", dest="minsize", required=False, default = 5, type = int,  metavar="minsize", help="minimum length of analyzed tryptic peptides")
parser.add_argument("-a", "--maxsize", dest="maxsize", required=False, default = 50, type = int,  metavar="maxsize", help="maximum length of analyzed tryptic peptides")
args = parser.parse_args()

handle=SeqIO.parse(args.input,'fasta')

## NOT deduplicating here, because we will deduplicate in the id_matchig step
#dedup_records = defaultdict(list)
dedup_records = []

for record in handle:
    proseq=str(record.seq)
    peptide_list=TRYPSIN(proseq,args.c)
    for peptide in peptide_list:
        if len(peptide) > args.minsize and len(peptide) < args.maxsize:
            dedup_records.append((str(peptide), record.id))
            # output.write("%s\t%s\n" % (record.id,peptide))
            #dedup_records[str(peptide)].append(record.id)

numberofpeps = len(dedup_records)

with open(args.output, 'w') as output:
#    for seq, ids in tqdm(dedup_records.items(), total=numberofpeps): 
    for seq in tqdm(dedup_records, total=numberofpeps): 
# Join the ids and write them out as the fasta
        #idjoin = ','.join(map(str, set(ids)))
        output.write(str(seq[0]) + "\t" + str(seq[1]) + "\n")
#        output.write(str(idjoin) + "\t" + str(seq) + "\n")


##
# check the annotatio for MH0389_GL0123916,DLM001_GL0023305,V1.FI31_GL0099642,MH0442_GL0198297
## see how this works. can you add your dedup code to this?

#for record in SeqIO.parse(args.input, "fasta"):
#    count = count + 1 
#    if count % 100000 == 0:
#        print("Reading sequence number " + str(count))
#    # Use the sequence as the key and then have a list of id's as the value
#    dedup_records[str(record.seq)].append(record.id)
#
#count = 0
#with open(args.output, 'w') as output:
#    for seq, ids in dedup_records.items():
#        count = count + 1
#        if count % 100000 == 0:
#            print("Writing sequence number " + str(count))
#        # Join the ids and write them out as the fasta
#        output.write(">" + "pep_id" + str(count) + " " + "Dup:" + str(len(ids))+"\n")
#        output.write(seq + "\n")
#

import os
import re
import strfuzzy_modified
import sys
import gzip
import resource
import pandas as pd
import numpy as np
import pdb
import collections
#import datatable


# TODO: add codecs block to make sure you are using the utf-8 for read_table
# import codecs
# BLOCKSIZE = 1048576 # or some other, desired size in bytes
# with codecs.open(sourceFileName, "r", "your-source-encoding") as sourceFile:
#     with codecs.open(targetFileName, "w", "utf-8") as targetFile:
#         while True:
#             contents = sourceFile.read(BLOCKSIZE)
#             if not contents:
#                 break
#             targetFile.write(contents)


# fuzzes
# left_primer_fuzz = int(sys.argv[2])
# right_primer_fuzz = int(sys.argv[3])

left_primer_fuzz = 1
right_primer_fuzz = 1

# if len(sys.argv) < 4:
#     raise SystemError("Error: You must specify arguments\n")
# file_fastq_R1 = sys.argv[1]

################################################################################
##input files
#file_fastq_R1="/home/magdalenabus/Magdalena/All_data_Sammed/Population_test_data/210422_M50235_0071_000000000-J9TF8_HIS_population/Data/Intensities/BaseCalls/8898_S12_L001_R1_001.fastq.gz"
file_fastq_R1="foo.txt"

file_primer="Updated_IDseek_SNP85_Primer_Seq_hg38_03232021_SM020922.txt"

mydata = pd.read_table(file_primer, delimiter="\t")

################################################################################


complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', \
              'R': 'R', 'Y': 'Y', 'S': 'S', 'W': 'W', \
              'K': 'K', 'M': 'M', 'B': 'B', 'D': 'D', \
              'H': 'H', 'V': 'V', 'N': 'N'}

def reverse_complement(seq):
    '''
    This function gives out a reverse complement of sequence
    @param seq - type of str
        The input seqeunce to construct a reverse complement of

    @return bases - type of str
        Reverse complement of input sequence
    '''
    bases = list(seq)
    bases = ''.join(complement[base] for base in reversed(bases))
    return bases



# def f(row):
#     val = row['Right primer']
#     if row['Left/right switch']==1:
#         val = reverse_complement(row['Right primer'])
#     return val

# works with above functions
# mydata['Right_recvomplement'] = mydata.apply(f, axis=1)

mydata['Right_revcomplement'] = mydata.apply(lambda x: \
reverse_complement(x['Right primer']) if (x['Left/right switch'] == 1) \
else x['Right primer'], axis=1)

mydata_dict = dict(zip(mydata['Left primer'],mydata['Right_revcomplement']))

#breakpoint()
numofLine = 0

TarPList = []
primernot = set()
primerin = set()
#pdb.set_trace()
with open(file_fastq_R1, "rt") as handleR1:
    while(True):
        r1= handleR1.readline()
        if not r1:
            break;
        # numofLine += 1
        # if numofLine == 4:
        #     numofLine = 0
        # if numofLine == 2:
        readR1 = r1.rstrip()
        for leftp, rightp in mydata_dict.items():
            leftp_idx = strfuzzy_modified.fuzzyFind(readR1, \
            leftp, fuzz = left_primer_fuzz)
            rightp_idx = strfuzzy_modified.fuzzyFind(readR1, \
            rightp, fuzz = right_primer_fuzz)
            # find all primers that are not found in the string, distinct 
            # of all those should be 85.
            if((leftp_idx == -1) and (rightp_idx == -1)):
                primernot.add(leftp)
            if((leftp_idx != -1) and (rightp_idx != -1)):
                primerin.add(leftp)
                tseq = readR1[leftp_idx[0]+len(leftp_idx[1]):rightp_idx[0]]
                # if(tseq=="TCCTTCGAGCTCCAGC"):
                #     pdb.set_trace()
                leftp_seq = leftp_idx[1]
                leftp_ham = leftp_idx[2]
                rightp_seq = rightp_idx[1]
                rightp_ham = rightp_idx[2]
                #pdb.set_trace()
                TarPList.append((leftp, leftp_seq, leftp_ham, \
                tseq, \
                rightp, rightp_seq, rightp_ham))  
                #TarPList.append((leftp_seq, tseq, rightp_seq))
                #TarPList.append(readR1)
                    
# print set to file
print(primernot)
print(primerin)

TarPListCount = collections.defaultdict(int)
for k in TarPList:
    TarPListCount[k] += 1
    
for k, v in TarPListCount.items():
    print("{}\t{}".format('\t'.join(map(str,k)),v))

# for k, v in TarPListCount.items():
#     print("{}\t{}".format(k,v))
#     
                
            #breakpoint()

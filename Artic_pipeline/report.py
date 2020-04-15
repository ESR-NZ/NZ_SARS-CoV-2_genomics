#!/usr/bin/env python

from Bio import SeqIO
import os
import re
from collections import Counter
import argparse

## Parse the args. Needs path to where the consensus fasta files are () 
parser = argparse.ArgumentParser()
parser.add_argument("--consensus_dir")

args = parser.parse_args()
CONSENSUN_DIR = args.consensus_dir
files_paths = [os.path.join(CONSENSUN_DIR, file) for file in os.listdir(CONSENSUN_DIR)] 


fasta_objs = [SeqIO.parse(fasta, 'fasta') for fasta in files_paths]


## Dictonary of data
sequences = {}
for fasta_obj in fasta_objs:
    for record in fasta_obj:
        metrics = {}
        ids=record.id.split('/')[0].split('_')[-1] 
        
        sequence = record.seq[54:-67]
        Ns_list = re.findall("N+", str(sequence))
        Num_gaps = len(Ns_list)
        gap_lenghts = [len(n) for n in Ns_list]
        hist = Counter(sorted(gap_lenghts))
        total_Ns = sequence.count('N')
        
        #metrics['sequence'] = sequence   
        metrics['seq_length'] = len(sequence)
        metrics['total_Ns'] = total_Ns
        metrics['num_gaps'] = len(Ns_list)
        metrics['hist'] = hist
        
        sequences[ids] = metrics

## Report
def report(metric_dict):
    
    f.write(f"Sequence Length: {metric_dict['seq_length']}\n")
    f.write(f"Number of gaps: {metric_dict['num_gaps']}\n")
    f.write(f"Total Ns: {metric_dict['total_Ns']}\n")
    f.write('\n')
    f.write('Gap histogram\n')
    for item in metric_dict['hist'].items():
        f.write(str(item[0]).ljust(5))
        f.write("#"*(item[1]))
        f.write('\n')
    
    f.write('\n\n\n')

    


##Write to file
f = open(os.path.join(CONSENSUN_DIR, "consensus_report.txt"), "a")
for k,v in sequences.items():
    f.write(k+'\n\n')
    report(v)
f.close() 

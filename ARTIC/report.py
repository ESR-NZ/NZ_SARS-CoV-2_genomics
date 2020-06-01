#!/usr/bin/env python

from Bio import SeqIO
import os, glob
import re
from collections import Counter
import argparse

## Parse the args. Needs path to where the consensus fasta files are () 
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--consensus_dir", help="Directory containing consensus genome files")
parser.add_argument("-b", "--bucket",        help="AWS bucket to upload genomes to", nargs='?', default=False)
parser.add_argument("-v", "--verbose",       help="Modify output verbosity", action="store_true")


args = parser.parse_args()
if not args.consensus_dir:
    parser.print_help()
    exit(0)

CONSENSUN_DIR = args.consensus_dir    
DEBUG = args.verbose

def GatherFiles(file_path, file_type="fasta"):
    file_list={}
    files_paths = [os.path.join(file_path, file) for file in glob.glob(file_path+"*."+file_type)]
    for fasta_file in files_paths:
        record = SeqIO.read(fasta_file, file_type)
        sid = record.id.split('/')[0].split('_')[-1]
        if DEBUG:
            print(record.id, sid)
        file_list[sid] = {}
        file_list[sid]["file"] = fasta_file
    return file_list

## Gather stats for sequences
def GatherStats(file_list, file_type="fasta"): 
    for sid in file_list:
        record = SeqIO.read(file_list[sid]["file"], file_type)
        if DEBUG:
            print(record)
        
        # Ignoring the first 54 and last 67 bases that are always N for SARS-CoV-2
        sequence = record.seq[54:-67] or ""

        # Gather N statistics
        Ns_list = re.findall("N+", str(sequence))
        Num_gaps = len(Ns_list)
        gap_lenghts = [len(n) for n in Ns_list]
        hist = Counter(sorted(gap_lenghts))
        total_Ns = sequence.count('N')
        
        if len(sequence) <= 1:
            print("[WARNING] empty sequence detected {}".format(sid))

        # Store in metric dictionary
        metrics = {} 
        #metrics['sequence'] = sequence   
        metrics['seq_length'] = len(sequence)
        metrics['total_Ns'] = total_Ns
        metrics['num_gaps'] = len(Ns_list)
        metrics['hist'] = hist
        metrics['pass_QC'] = total_Ns < 3000 and len(sequence) > 29000
        
        # Add to global dictionary
        file_list[sid]["metrics"] = metrics

## Report
def Report(metric_dict, fstream):
    fstream.write("""Sequence Length: {seq_length}
    Number of gaps: {num_gaps}
    Total Ns: {total_Ns}
    """.format(**metric_dict))
    
    fstream.write('Gap histogram\n')
    for item in metric_dict['hist'].items():
        fstream.write(str(item[0]).ljust(5))
        fstream.write("#"*(item[1]))
        fstream.write('\n')
    
    fstream.write('----\n')


## Upload genomes that pass QC
def Upload(file_list, bucket):
    if not bucket:
        print("Not uploading data")
        return 1
    
    #guix package -i awscli
    #os.system("guix package -i awscli")

    for sid in file_list:
        if DEBUG:
            print(sid, file_list[sid]["metrics"]["pass_QC"])
        if file_list[sid]["metrics"]["pass_QC"]:
            os.system("aws s3 cp {0} {1}".format(file_list[sid]["file"], bucket))
            os.system("cp {0} {1}".format(file_list[sid]["file"], "/NGS/active/VIR/SARS-CoV2/results/consensus/PASS_QC/"))
            if DEBUG:
                print("aws s3 cp {0} {1}".format(file_list[sid]["file"], bucket))
    #pip3 install nextstrain-cli --user
    #nextstrain remote upload {bucket} {file(s)}

seqstats = GatherFiles(CONSENSUN_DIR)
GatherStats(seqstats)

##Write to file
f = open(os.path.join(CONSENSUN_DIR, "consensus_report.txt"), "a")
print("Writing reports")
for k in seqstats:
    f.write('{} - {}:\n'.format(k, seqstats[k]["file"]))
    Report(seqstats[k]["metrics"], f)
f.close()

print("Uploading PASS genomes")
Upload(seqstats, args.bucket)

"""
A simple script for nominating protospacers based on cas and fasta file.
"""

import argparse
import yaml
from Bio import SeqIO

import sys
sys.path.append('code')

from crispr import CasTemplate
from collection import HitCollection


def main(cas, in_path, csv_path, top_n):
    
    with open(in_path) as handle:
        seqs = SeqIO.parse(handle, 'fasta')
        hit_collection = HitCollection.from_sequences(cas, seqs, top_N=top_n)
    
    
    with open(csv_path, 'w') as handle:
        hit_collection.to_csv(handle)
        
    yaml_path = csv_path.replace('.csv', '.yaml')
    with open(yaml_path, 'w') as handle:
        handle.write(f'min_count: {min(hit_collection.counts.values())}\n')
        handle.write(f'max_count: {max(hit_collection.counts.values())}\n')
        


def parse_args():
    parser = argparse.ArgumentParser(description='Nominate gRNAs from a fasta-file')
    parser.add_argument('pam', type=str, help='gRNA pattern')
    parser.add_argument('cas_type', type=str, help='Type II or V to indicate PAM/PS relationship')
    parser.add_argument('fasta_file', type=str, help='Path to the fasta-file input')
    parser.add_argument('csv_file', type=str, help='Path to the csv-file output of nominated gRNAs')
    parser.add_argument('--top', type=int, default=1000, help='Number of top nominations to return')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    
    cas = CasTemplate(pam_motif=args.pam,
                      typ = args.cas_type,
                      protospacer_length=20)
    
    main(cas, args.fasta_file, args.csv_file, args.top)
    

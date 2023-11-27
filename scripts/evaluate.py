"""
A simple script for evaluating nominated protospacers based on cas and fasta file.
"""

import argparse
import yaml

import sys
sys.path.append('code')

from crispr import CasTemplate
from collection import HitCollection


def main(cas, nominate_path, fasta_path, csv_path, mismatch):
    
    with open(nominate_path) as handle:
        hit_collection = HitCollection.from_csv(handle, cas=cas)
    
    res = hit_collection.offinder_search(fasta_path, mismatch = mismatch)
    
    res.to_csv(csv_path, index=False)
    yaml_path = csv_path.replace('.csv', '.yaml')
    
    detected_protospacers = len(res['protospacer'].unique())
    
    with open(yaml_path, 'w') as handle:
        handle.write(f'total_hits: {len(res.index)}\n')
        handle.write(f'detected_protospacers: {detected_protospacers}\n')
        


def parse_args():
    parser = argparse.ArgumentParser(description='Evaluate nominated gRNAs against a fasta-file')
    parser.add_argument('pam', type=str, help='gRNA pattern')
    parser.add_argument('cas_type', type=str, help='Type II or V to indicate PAM/PS relationship')
    parser.add_argument('nominate_file', type=str, help='Path to the nominated protospacers.')
    parser.add_argument('fasta_file', type=str, help='Path to the fasta-file search')
    parser.add_argument('csv_file', type=str, help='Path to the csv-file output of nominated gRNAs')
    parser.add_argument('--mismatch', type=int, default=3, help='Number of mismatches allowed')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    
    protospacer_start = -20 if args.cas_type == 'II' else len(args.pam)
    cas = CasTemplate(pam_motif=args.pam,
                      typ = args.cas_type,
                      protospacer_length=20)
    
    main(cas, args.nominate_file, args.fasta_file, args.csv_file, args.mismatch)
    
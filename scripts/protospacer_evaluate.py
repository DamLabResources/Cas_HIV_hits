"""
A script for evaluating one or more protospacers.
"""

import sys
sys.path.append('code')


import argparse

from crispr import CasTemplate
from pipeline import NDNFPipeline


def main(cas_type, pam, protospacers, reference, diversify=0, ontarget=None, offtarget=None):
    
    cas = CasTemplate(pam_motif = pam, typ=cas_type, 
                      protospacer_length=20)
    
    pipeline = NDNFPipeline.from_protospacers(cas, protospacers)
    
    if diversify > 0:
        pipeline = pipeline.diversify(diversify)
    
    mapping = pipeline.reference_map(reference, enforce_pam = False)
    
    if ontarget is not None:
        _, protospacer_hist_freq = pipeline.narrow(ontarget, return_freqs=True)
        mapping['on_target_freq'] = mapping['protospacer'].map(protospacer_hist_freq.get)
        
    if offtarget is not None:
        _, off_target_counts = pipeline.filter(offtarget, return_counts=True)
        mapping['off_target_counts'] = mapping['protospacer'].map(off_target_counts.get)
        
    return mapping


def parse_args():
    parser = argparse.ArgumentParser(description='Evaluate protospacers')
    
    parser.add_argument('cas_type', type=str, help='Type II or V to indicate PAM/PS relationship',
                        default = 'II')
    parser.add_argument('pam', type=str, help='PAM',
                        default = 'NGG')
    
    parser.add_argument('protospacers', type=str, help='Path to a nomination fasta dataset.',
                        nargs='+')
    
    parser.add_argument('--ontarget_file', type=str, help='Path to a testing fasta dataset.')
    parser.add_argument('--offtarget_file', type=str, help='Path to a offtarget fasta dataset.')
    
    parser.add_argument('--reference', type=str, help='Reference genome to use for mapping.',
                        default = None)
    
    parser.add_argument('--diversify', type=int, help='Add N random mutations to the provided protospacers before analysis.',
                        default = 0)
    
    return parser.parse_args()


if __name__ == '__main__':
    
    args = parse_args()
    
    mapping = main(args.cas_type, args.pam, args.protospacers, 
                   args.reference, diversify=args.diversify, 
                   ontarget=args.ontarget_file, offtarget=args.offtarget_file)
    
    mapping.to_csv(sys.stdout, index=False)
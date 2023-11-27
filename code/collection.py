"""
A set of classes for handling the processing of many sequences and collection highly conserved protospacers
"""

from collections import defaultdict
from dataclasses import dataclass
import regex as re
from crispr import Hit, reverse_complement
from tqdm import tqdm
from itertools import groupby
from tempfile import NamedTemporaryFile
from subprocess import check_call
import shlex
import pandas as pd
import numpy as np
import csv

def seqrecord2hitstream(cas, sequences):
    
    top_re = re.compile(cas.to_regexp())
    bot_re = re.compile(cas.to_regexp(rc=True))
    for seqR in tqdm(sequences):

        # Forward
        for match in top_re.finditer(str(seqR.seq), overlapped=True):
            gdict = match.groupdict()
            yield Hit(cas, 
                      protospacer=gdict['spacer'], 
                      pam=gdict['pam'], 
                      span=match.span(),
                      strand = '+',
                      origin=seqR.id)
            

        # Reverse
        for match in bot_re.finditer(str(seqR.seq), overlapped=True):
            gdict = match.groupdict()
            yield Hit(cas, 
                      protospacer=reverse_complement(gdict['spacer']), 
                      pam=reverse_complement(gdict['pam']), 
                      span=match.span(),
                      strand = '-',
                      origin=seqR.id)


class HitCollection(object):
    
    def __init__(self, cas, counts, representatives):
        
        self.cas = cas
        self.counts = counts
        self.representatives = representatives
        
    
    @staticmethod
    def from_stream(cas, hitstream, top_N = 1000, filters=None):
        
        counts = defaultdict(int)
        reps = {}
        
        for hit in hitstream:
            if filters is not None:
                if any(not f(hit) for f in filters): continue
            counts[hit.protospacer] += 1
            if hit.protospacer not in reps:
                reps[hit.protospacer] = hit
        
        order = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        
        for spacer, _ in order[top_N:]:
            counts.pop(spacer)
            reps.pop(spacer)
        
        return HitCollection(cas, counts, reps)
    
    @staticmethod
    def from_sequences(cas, sequences, top_N = 1000, filters='default'):
        
        if filters == 'default':
            filters = [run_length_filter, gc_content_filter]
        
        stream = seqrecord2hitstream(cas, sequences)
        return HitCollection.from_stream(cas, stream, top_N = top_N, filters=filters)
    
    @staticmethod
    def from_csv(handle, cas=None):
        
        reader = csv.DictReader(handle)
        counts, reps = {}, {}
        for row in reader:
            reps[row['protospacer']] = Hit(cas=cas,
                                           protospacer=row['protospacer'],
                                           pam=row['pam'],
                                           span=(int(row['span_start']), int(row['span_stop'])),
                                           strand = row['strand'],
                                           origin = row['origin'])
            counts[row['protospacer']] = int(row['counts'])
        
        return HitCollection(cas, counts, reps)
    
    def limit(self, protospacers):
        """Return a new collection that is limited to these protospacers."""
        
        counts = dict((key, self.counts.get(key, 0)) for key in protospacers)
        reps = dict((key, self.representatives[key]) for key in protospacers)
        
        return HitCollection(self.cas, counts, reps)
    
    def diversify(self, target_size, random_state=None):
        
        if random_state is None:
            random_state = np.random.RandomState()
        
        t = tqdm(total=target_size) # Initialise
        t.update(len(self.representatives))
        
        while len(self.representatives) < target_size:
            hit = random_state.choice(list(self.representatives.values()))
            new_hit = hit.diversify()
            if new_hit.protospacer not in self.representatives:
                self.representatives[new_hit.protospacer] = new_hit
                t.update(1)
        t.close()
        
        
    
    
    def to_rows(self):
        
        for proto, count in self.counts.items():
            info = {'protospacer': proto, 
                    'counts': count}
            rep = self.representatives[proto]
            info['pam'] = rep.pam
            info['span_start'] = rep.span[0]
            info['span_stop'] = rep.span[0]
            info['strand'] = rep.strand
            info['origin'] = rep.origin
            yield info
        
        
    
    def to_df(self):
        
        return pd.DataFrame(list(self.to_rows()))
    
    
    
    def to_csv(self, handle):
        
        fieldnames = ['protospacer', 'pam', 'span_start', 'span_stop', 'strand', 'origin', 'counts']
        writer = csv.DictWriter(handle, fieldnames)
        writer.writeheader()
        writer.writerows(self.to_rows())
            
    
    def offinder_search(self, search_path, dna_bulge=2, rna_bulge=1, mismatch=5):
        """Runs cas-offinder and returns a dataframe of results"""
        
        with NamedTemporaryFile(suffix='.txt', mode='w') as in_handle:
            self.to_cas_offinder(search_path, in_handle, 
                                 dna_bulge=dna_bulge,
                                 rna_bulge=rna_bulge,
                                 mismatch=mismatch)
            in_handle.flush()
            with NamedTemporaryFile(suffix='.tsv', mode='r') as out_handle:
                cmd = f'cas-offinder {in_handle.name} C0 {out_handle.name}'
                
                check_call(shlex.split(cmd))
                df = pd.read_csv(out_handle.name, sep='\t',header=None,
                                 names = ['protospacer', 'chrom', 'start', 
                                          'hit', 'strand', 'mismatches'])
                df['protospacer'] = df['protospacer'].map(lambda x: x.replace('N', ''))
                return df
                
    def offinder_counts(self, search_path, name='hits', one_hit_per_chrom=True, dna_bulge=2, rna_bulge=1, mismatch=5):
        """returns a series indexed by protospacer that indicates the number of counts found from the offinder search"""
        
        hits = self.offinder_search(search_path, dna_bulge=dna_bulge, rna_bulge=rna_bulge, mismatch=mismatch)
        
        if one_hit_per_chrom:
            protospacer_hits = hits.groupby('protospacer')['chrom'].nunique()
        else:
            hits['chromstart'] = hits.apply(lambda row: row['chrom']+str(row['start']), axis=1)
            protospacer_hits = hits.groupby('protospacer')['chromstart'].nunique()
        
        # missing ones have 0 hits
        protospacer_hits = protospacer_hits.reindex(list(self.representatives.keys())).fillna(0)
        protospacer_hits.name = name
        
        return protospacer_hits
    
    
    
    
    def to_cas_offinder(self, search_path, handle, dna_bulge=2, rna_bulge=1, mismatch=5):
        """Writes a cas-offinder compatable file with this collection.
        https://github.com/snugel/cas-offinder
        Example input file (DNA bulge size 2, RNA bulge size 1):
        /var/chromosomes/human_hg38
        NNNNNNNNNNNNNNNNNNNNNRG 2 1
        GGCCGACCTGTCGCTGACGCNNN 5
        CGCCAGCGTCAGCGACAGGTNNN 5
        ACGGCGCCAGCGTCAGCGACNNN 5
        GTCGCTGACGCTGGCGCCGTNNN 5
        """
        
        handle.write(search_path + '\n')
        
        if self.cas.typ == 'II':
            template = 'N'*self.cas.protospacer_length + self.cas.pam_motif
        else:
            template = self.cas.pam_motif + 'N'*self.cas.protospacer_length 
            
        handle.write(f'{template} {dna_bulge} {rna_bulge}\n')
        
        pam_suffix = 'N'*len(self.cas.pam_motif)
        for rep in self.representatives.values():
            if self.cas.typ == 'II':
                target = rep.protospacer + pam_suffix
            else:
                target = pam_suffix + rep.protospacer
            handle.write(f'{target} {mismatch}\n')
            
        
    
        


def gc_content(protospacer):
    return sum((s == 'C') | (s == 'G') for s in protospacer)/len(protospacer)

def gc_content_filter(hit, mx=0.6, mn=0.4):
    """ Returns True if the gc content is within [mn, mx]
    implying the protospacer should be kept"""
    
    gc = gc_content(hit.protospacer)
    return (gc>mn) & (gc<mx)

def run_length(protospacer):
    """Calculate the longest run"""
    
    longest = 0
    for _, run in groupby(protospacer):
        this = len(list(run))
        if this > longest:
            longest = this
    return longest


def run_length_filter(hit, max_run = 4):
    return run_length(hit.protospacer)<=max_run

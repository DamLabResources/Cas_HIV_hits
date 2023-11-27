from collection import defaultdict
from itertools import groupby
from copy import deepcopy

import seaborn as sns
import regex as re
from tqdm import tqdm
import numpy as np
from Bio import SeqIO

from crispr import Hit, CasTemplate, reverse_complement
from casoffinder import offinder_search
import visualize


class NDNFPipeline(object):
    
    def __init__(self, cas, protospacer_index):
        
        self.cas = cas
        self.protospacer_index = protospacer_index
        
        
    @staticmethod
    def from_protospacers(cas, protospacers):
        
        hit_stream = []
        for proto in protospacers:
            hit_stream.append(Hit(cas=cas,
                                  protospacer = proto,
                                  pam = cas.pam_motif,
                                  span = (None, None),
                                  strand = '+', origin = 'cmdline'))
        _, full_index = hitstream2index(hit_stream, [])
        return NDNFPipeline(cas, full_index)
                              
            
    
    @staticmethod
    def nominate(cas, sequences, top_n = 1000, filters = 'default',
                 make_figure = False, return_counts = False):
        """Nominate a set of protospacers from sequences."""
        
        if filters == 'default':
            filters = [run_length_filter, gc_content_filter]
            filters = []
        
        hit_stream = sequences2hitstream(cas, sequences)
        counts, full_index = hitstream2index(hit_stream, filters)
            
        order = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        nominated_index = dict((proto, full_index[proto]) for proto, _ in order[:top_n])
        
        if return_counts:
            
            return NDNFPipeline(cas, nominated_index), counts
        
        return NDNFPipeline(cas, nominated_index)
    
    def diversify(self, target_size=10_000, random_state=None,
                  make_figure = False):
        """Generate new protospacers through random mutation"""
        
        if random_state is None:
            random_state = np.random.RandomState()
        
        new_index = deepcopy(self.protospacer_index)
        
        # Create a progress bar
        t = tqdm(total=target_size) # Initialise
        t.update(len(new_index))
        
        # Run until we have created enough diversity
        while len(new_index) < target_size:
            # Pick a random hit
            hit = random_state.choice(list(new_index.values()))
            
            # Mutate it until new
            new_hit = hit.diversify()
            while new_hit.protospacer in new_index:
                new_hit = new_hit.diversify()
            
            # Add it
            new_index[new_hit.protospacer] = new_hit
            t.update(1)
        t.close()
        
        def num_muts(name):
            return sum(l == '+' for l in name)
        
        mut_counts = [num_muts(rep.origin) for rep in new_index.values()]
        
        if make_figure:
            fig = visualize.mutation_count_figure(mut_counts)
            return NDNFPipeline(self.cas, new_index), fig
        
        return NDNFPipeline(self.cas, new_index)
    
    def narrow(self, target_population_fasta, min_rate = 0.9, mismatch=3,
               make_figure=False, return_freqs = False):
        """Only keep protospacers that target at least {min_rate) sequences"""
        
        # Count how many hits there are
        with open(target_population_fasta) as handle:
            for nunique_chroms, _ in enumerate(SeqIO.parse(handle, 'fasta')):
                pass
        
        
        # Search for hits
        protospacer_hits = offinder_search(self.cas.pam_motif, self.cas.typ, self.cas.protospacer_length, 
                                           target_population_fasta, 
                                           list(self.protospacer_index.keys()),
                                           mismatch=mismatch)
        
        # Count how many hits each protospacer has
        top_mm_counts = protospacer_hits.groupby('protospacer')['chrom'].nunique()
        top_mm_counts = top_mm_counts.reindex(list(self.protospacer_index.keys())).fillna(0)
        
        # Only include ones over threshold
        freq = top_mm_counts/nunique_chroms
        valid_protos = freq[freq>min_rate].index
        
        index = dict((proto, self.protospacer_index[proto]) for proto in valid_protos)
        
        if return_freqs:
            return NDNFPipeline(self.cas, index), freq
            
        return NDNFPipeline(self.cas, index)
        
    
    def filter(self, off_target_fasta, mismatch=2, max_hits = 0, make_figure=False,
               return_counts = False):
        """Remove protospacers with hits in fasta."""
        
        # Search protospacers
        protospacer_hits = offinder_search(self.cas.pam_motif, self.cas.typ, self.cas.protospacer_length, 
                                           off_target_fasta, 
                                           list(self.protospacer_index.keys()),
                                           mismatch=mismatch)
        # Count ones that occur
        top_mm_counts = protospacer_hits.groupby('protospacer')['chrom'].nunique()
        
        # Reindex to include ones that did not (and set them to 0 hits)
        hits = top_mm_counts.reindex(list(self.protospacer_index.keys())).fillna(0)
        
        # Filter to allowed ones
        valid_protos = hits[hits<=max_hits].index
        
        index = dict((proto, self.protospacer_index[proto]) for proto in valid_protos)
        if return_counts:
            return NDNFPipeline(self.cas, index), hits
            
        return NDNFPipeline(self.cas, index)

    
    def reference_map(self, ref_fasta, mismatch=4, enforce_pam=True):
        
        pam = self.cas.pam_motif if enforce_pam else 'N'*len(self.cas.pam_motif)
        
        return offinder_search(pam, self.cas.typ, self.cas.protospacer_length, 
                               ref_fasta, 
                               list(self.protospacer_index.keys()),
                               mismatch=mismatch)

    
def full_pipeline(cas, sequences,
                  on_target_fasta,
                  off_target_fasta, 
                  on_target_mm = 3, on_target_minrate = 0.5,
                  off_target_mm = 3, off_target_allowed = 0,
                  nominate_n = 20_000,
                  diversify_n = 30_000,
                  save_path=None):
    
    nominated = NDNFPipeline.nominate(cas, sequences, top_n=nominate_n, make_figure=False)
    print('Nominated:', len(nominated.protospacer_index))
    diversed = nominated.diversify(diversify_n, make_figure=False)
    print('Diversified to:', len(diversed.protospacer_index))
    narrowed = diversed.narrow(on_target_fasta, min_rate=on_target_minrate, mismatch=on_target_mm)
    print('Narrowed to:', len(narrowed.protospacer_index), 'with a rate>=', on_target_minrate)
    if len(narrowed.protospacer_index) == 0:
        return narrowed
    safe = narrowed.filter(off_target_fasta, mismatch=off_target_mm, max_hits = off_target_allowed)
    print('Filtered to:', len(safe.protospacer_index), 'with a offtarget<=', off_target_allowed)
    
    if save_path is not None:
        with open(save_path, 'w') as handle:
            handle.write(f'num_broad_spec: {len(narrowed.protospacer_index)}\n')
            handle.write(f'num_safe: {len(safe.protospacer_index)}\n')
            if len(safe.protospacer_index):
                handle.write('safe:\n')
                for proto in safe.protospacer_index.values():
                    handle.write(f' - protospacer: {proto.protospacer}\n')
                    handle.write(f'   origin: {proto.origin}\n')
                    handle.write(f'   span: {proto.span}\n')
    
    return safe
    
    
    
        
def sequences2hitstream(cas, sequences):
    
    VALID_LETS = set(['A', 'C', 'G', 'T'])
    top_re = re.compile(cas.to_regexp())
    bot_re = re.compile(cas.to_regexp(rc=True))
    for seqR in tqdm(sequences):

        # Forward
        for match in top_re.finditer(str(seqR.seq), overlapped=True):
            gdict = match.groupdict()
            if all(l.upper() in VALID_LETS for l in gdict['spacer']):
                yield Hit(cas, 
                          protospacer=gdict['spacer'], 
                          pam=gdict['pam'], 
                          span=match.span(),
                          strand = '+',
                          origin=seqR.id)


        # Reverse
        for match in bot_re.finditer(str(seqR.seq), overlapped=True):
            gdict = match.groupdict()
            if all(l.upper() in VALID_LETS for l in gdict['spacer']):
                yield Hit(cas, 
                          protospacer=reverse_complement(gdict['spacer']), 
                          pam=reverse_complement(gdict['pam']), 
                          span=match.span(),
                          strand = '-',
                          origin=seqR.id)
            
def hitstream2index(hit_stream, filters):
    """Consume a stream of hits and create an index and counts"""
        
    counts = defaultdict(int)
    index = {}
    for hit in hit_stream:
        if filters is not None:
            if any(not f(hit) for f in filters): continue
        if hit.protospacer not in index:
            index[hit.protospacer] = hit
        counts[hit.protospacer] += 1
        
    return counts, index
            
            
            
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

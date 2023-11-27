"""
A place to hold the CRISPR classes used for searching
"""
from dataclasses import dataclass
import regex as re
import yaml
import numpy as np


@dataclass
class CasTemplate:
    """Class describing generic Cas.
    All positions are relative to the 5' end of the pam
    spCas9
    pam_motif: NGG
    typ: II
    protopacer_length: 20
    """
    pam_motif: str
    typ: str
    protospacer_length: int
    name: str = None
    
    @classmethod
    def from_yaml(cls, handle):
        
        info = yaml.full_load(handle)
        return cls(pam_motif=info['pam_motif'],
                   typ=int(info['typ']),
                   protospacer_length=int(info['protospacer_length']),
                   name=info.get('name'))
    
    def to_regexp(self, rc = False):
        """Convert this to a searchable regexp."""
        
        pam_reg = motif2regex(self.pam_motif, rc)
        pam_reg = f'(?<pam>{pam_reg})'
        
        ln = self.protospacer_length
        spacer_reg = '(?P<spacer>.{' + f'{ln},{ln}' + '})'
            
        items = [spacer_reg]
        items.append(pam_reg)
        
        if ((self.typ == 'II') & (not rc)) | (((self.typ == 'V') & rc)):
            items = items[::-1]
        
        return ''.join(items)
    
    
    def search_string(self, sequence, origin=None, top=True, bottom=True):
        
        
        if top:
            reg = re.compile(self.to_regexp())
            for match in reg.finditer(sequence, overlapped=True):
                gdict = match.groupdict()
                yield Hit(self, 
                          protospacer=gdict['spacer'], 
                          pam=gdict['pam'], 
                          span=match.span(),
                          strand = '+', origin=origin)
        if bottom:
            reg = re.compile(self.to_regexp(rc=True))
            for match in reg.finditer(sequence, overlapped=True):
                gdict = match.groupdict()
                yield Hit(self, 
                          protospacer=reverse_complement(gdict['spacer']), 
                          pam=reverse_complement(gdict['pam']), 
                          span=match.span(),
                          strand = '-', origin=origin)
                
        
@dataclass
class Hit:
    cas : CasTemplate
    protospacer : str
    pam : str
    span : (int, int)
    strand: str
    origin: str = None
    
    @property
    def gc_content(self):
        return sum(l == 'G' or l == 'C' for l in protospacer) 
    
    
    @property
    def max_run(self):    
        max_run = 0
        for letter, run in groupby(self.protospacer):
            ln = len(list(run))
            if ln > max_run:
                max_run = kn
        return max_run
    
    def diversify(self, random_state=None):
        
        if random_state is None:
            random_state = np.random.RandomState()
        
        mut_site = random_state.randint(len(self.protospacer))
        nucs = list(self.protospacer)
        target = nucs[mut_site]
        mut = random_state.choice(list(NT2MUT[nucs[mut_site]]))
        nucs[mut_site] = mut
        
        if self.origin is not None:
            origin = self.origin + f' + {target}{mut_site}{mut}'
        else:
            origin = f'{target}{mut_site}{mut}'
            
        
        return Hit(cas = self.cas,
                   protospacer = ''.join(nucs),
                   pam = self.pam,
                   span = self.span,
                   strand = self.strand,
                   origin = origin)
        
    
    
    
    
SpCas9 = CasTemplate(pam_motif = 'NGG',
                     typ = 'II',
                     protospacer_length = 20)


AacCas12b = CasTemplate(pam_motif = 'TTN',
                        typ = 'V',
                        protospacer_length = 20)


NT2MUT = {'A': 'CGT',
          'C': 'AGT',
          'G': 'ACT',
          'T': 'ACG'}



NT2RE = {
 'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'T',
    'M': '[AC]',
    'R': '[AG]',
    'W': '[AT]',
    'S': '[CG]',
    'Y': '[CT]',
    'K': '[GT]',
    'V': '[ACG]',
    'H': '[ACT]',
    'D': '[AGT]',
    'B': '[CGT]',
    'N': '[ACGT]'
}

NT2REV = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'U': 'A',
    'M': '[GT]',
    'R': '[CT]',
    'W': '[AT]',
    'S': '[CG]',
    'Y': '[GA]',
    'K': '[CA]',
    'V': '[CGT]',
    'H': '[ACT]',
    'D': '[AGT]',
    'B': '[ACG]',
    'N': '[ACGT]'
}





def motif2regex(motif, rc):
    
    if rc:
        return ''.join(NT2REV[m] for m in motif[::-1])
    else:
        return ''.join(NT2RE[m] for m in motif)
    

def reverse_complement(sequence):
    
    return ''.join(NT2REV[s] for s in sequence[::-1])
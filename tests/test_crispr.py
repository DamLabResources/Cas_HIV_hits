import sys
sys.path.append('code')

import crispr


def test_motif2regex():
    
    
    reg = crispr.motif2regex('NGG', False)
    expected = '[ACGT]GG'
    assert reg == expected, f'Expected {expected}, got {reg}' 
    
    reg = crispr.motif2regex('NGG', True)
    expected = 'CC[ACGT]'
    assert reg == expected, f'Expected {expected}, got {reg}' 
    
    
    
    
def test_spcas9_to_regexp():
    
    expected = '(?P<spacer>.{20,20})(?<pam>[ACGT]GG)'
    reg = crispr.SpCas9.to_regexp()
    assert reg == expected, f'Expected {expected}, got {reg}' 
    
    expected = '(?<pam>CC[ACGT])(?P<spacer>.{20,20})'
    reg = crispr.SpCas9.to_regexp(rc=True)
    assert reg == expected, f'Expected {expected}, got {reg}' 
    
    
    
def test_AacCas12b_to_regexp():
    
    expected = '(?<pam>TT[ACGT])(?P<spacer>.{20,20})'
    reg = crispr.AacCas12b.to_regexp()
    assert reg == expected, f'Expected {expected}, got {reg}' 
    
    expected = '(?P<spacer>.{20,20})(?<pam>[ACGT]AA)'
    reg = crispr.AacCas12b.to_regexp(rc=True)
    assert reg == expected, f'Expected {expected}, got {reg}' 
    

def test_reverse_complement():
    
    seq = 'AATTTGGGCCCC'
    expected = 'GGGGCCCAAATT'
    
    revcomp = crispr.reverse_complement(seq)
    
    assert revcomp == expected, f'Expected {expected}, got {revcomp}' 
    
    
def test_spcas9_search_sequence():
    
    sequence = 'AGTAGCCGATGATCGATAGCTAGCTAGCGGATGCTAGATAGCTAGCTAGCAT'
                       
    hits = list(crispr.SpCas9.search_string(sequence))
    
    assert len(hits) == 2, f'Expected 2 hits, got {len(hits)}'
    
    
    hit1 = crispr.Hit(cas = crispr.SpCas9,
                      protospacer='GATGATCGATAGCTAGCTAG',
                      pam = 'CGG',
                      span = (7, 30),
                      strand='+')
    
    assert hits[0] == hit1
    
    
    hit2 = crispr.Hit(cas = crispr.SpCas9,
                      protospacer= crispr.reverse_complement('ATGATCGATAGCTAGCTAGC'),
                      pam = 'CGG',
                      span = (5, 28),
                      strand = '-')
    
    assert hits[1] == hit2
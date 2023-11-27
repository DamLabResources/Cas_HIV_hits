import sys
sys.path.append('code')

import crispr
import collection


def test_collection_from_stream():
    
    
    seqA = 'A'*50 + 'TGAGCTAGCGATCGATCGAT' + 'AGG' + 'T'*50
    seqB = 'A'*90 + 'CGAGCATAGCTTAGCTAATC' + 'AGG' + 'T'*50
    
    stream = [seqA]*50 + [seqB]*100
    
    hits = collection.HitCollection.from_sequences(crispr.SpCas9, stream)
    
    print(hits.counts)
    assert len(hits.counts) == 2
    assert hits.counts['TGAGCTAGCGATCGATCGAT'] == 50
    assert hits.counts['CGAGCATAGCTTAGCTAATC'] == 100
    
    assert hits.representatives['TGAGCTAGCGATCGATCGAT'].span == (50, 73)
    assert hits.representatives['CGAGCATAGCTTAGCTAATC'].pam == 'AGG'
    
    

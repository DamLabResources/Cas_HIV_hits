import pandas as pd
from tempfile import NamedTemporaryFile
from subprocess import check_call
import shlex


def offinder_search(pam_motif, typ, protospacer_length, search_path, 
                    protospacers,
                    dna_bulge=2, rna_bulge=1, mismatch=5):
    """Runs cas-offinder and returns a dataframe of results"""

    with NamedTemporaryFile(suffix='.txt', mode='w') as in_handle:
        cas_offinder_template_writer(pam_motif, typ, protospacer_length,
                                     search_path, 
                                     protospacers,
                                     in_handle, 
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


def cas_offinder_template_writer(pam_motif, typ, protospacer_length,
                                 search_path, 
                                 protospacers,
                                 handle, dna_bulge=2, rna_bulge=1, mismatch=5):
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

    if typ == 'II':
        template = 'N'*protospacer_length + pam_motif
    else:
        template = pam_motif + 'N'*protospacer_length 

    handle.write(f'{template} {dna_bulge} {rna_bulge}\n')

    pam_suffix = 'N'*len(pam_motif)
    for spacer in protospacers:
        if typ == 'II':
            target = spacer + pam_suffix
        else:
            target = pam_suffix + spacer
        handle.write(f'{target} {mismatch}\n')

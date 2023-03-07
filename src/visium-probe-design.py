#! /usr/bin/env python3

'''
visium-probe-design.py

Design a set of ligation-based probes for the 10x Genomics Visium fixed RNA
assay.

Input is a FASTA of mRNA sequence targets. Guidelines for these sequeces are in
the README. 

Default output is a TSV of probe designs to be synthesized.

Specify `--output_fasta` to generate a FASTA of hybridization portions that can
be compared to genomic sequence.

Specify `--idt` to generate probes in a format for easy pasting into an IDT
ordering spreadsheet.

Example:

$ python src/visium-probe-design.py test/test.fa

A PDF from 10x Genomics with design guidelines is in `/doc`.
'''

__author__ = 'Jay Hesselberth <jay.hesselberth@cuanschutz.edu>'
__license__ = 'MIT'

from collections import Counter, defaultdict

import click
import rle

from Bio import SeqIO

PROBE_LHS = 'CCTTGGCACCCGAGAATTCCA'
PROBE_RHS = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
MOD_RHS   = '/5Phos/'

IDT_SCALE = '25nm'
IDT_PURIF = 'STD'

PROBE_LEN = 25
GC_MAX    = 0.73
GC_MIN    = 0.44
HOMOP_MAX = 3

ALLOWED_NUCS = set(['A','G','T','C']) 

@click.command()
@click.option('--output_fasta', default=False, is_flag=True)
@click.option('--idt', default=False, is_flag=True)
@click.argument('target_fasta')

def main(target_fasta, output_fasta, idt):

    # generate all possible probes and identify candidates
    cand_probes = defaultdict(list)
    for record in SeqIO.parse(open(target_fasta), format = 'fasta'):
        rc_seq = record.seq.reverse_complement()

        for pos, lhs, rhs in gen_probe_pairs(rc_seq):

            if not _check_lhs_seq(lhs):
                continue
            if not all(map(_check_gc, [lhs, rhs])):
                continue
            if not all(map(_check_homopolymer, [lhs, rhs])):
                continue

            cand_probes[record.id].append(ProbePair(lhs, rhs, pos))

    # evaluate candidate probes to ensure spacing between them 
    keep_probes = defaultdict(list)
    for id, probes in cand_probes.items():

        last_pos = None
        for probe in probes:
            # ensure probe pairs are at least 50 nt away from another
            # probe pair
            if last_pos and abs(probe.lhs_start - last_pos) < (PROBE_LEN * 2):
                continue

            keep_probes[id].append(probe)
            last_pos = probe.lhs_start

    # report the probes in TSV or FASTA format
    if not output_fasta:
        if idt:
            colnames = ['#id', 'sequence', 'scale', 'purification']
            print('\t'.join(colnames))
        else:
            colnames = ['#id', 'pos', 'hyb_region', 'probe']
            print('\t'.join(colnames))

    for id, probes in keep_probes.items():
        for probe in probes:

            if output_fasta:
                print(f'>{id}-{probe.lhs_start}-lhs\n{probe.lhs}')
                print(f'>{id}-{probe.lhs_start}-rhs\n{probe.rhs}')
            elif idt:
                fields_left = [id+'-'+str(probe.lhs_start)+'-lhs', PROBE_LHS + probe.lhs,
                               IDT_SCALE, IDT_PURIF]
                fields_right = [id+'-'+str(probe.lhs_start)+'-rhs', MOD_RHS + probe.rhs + PROBE_RHS,
                               IDT_SCALE, IDT_PURIF]

                print('\t'.join(map(str, fields_left)))
                print('\t'.join(map(str, fields_right)))
            else:
                fields_left = [id+'-'+str(probe.lhs_start)+'-lhs', probe.lhs_start,
                  probe.lhs,
                  PROBE_LHS + probe.lhs]
                fields_right = [id+'-'+str(probe.lhs_start)+'-rhs', probe.lhs_start,
                  probe.rhs,
                  MOD_RHS + probe.rhs + PROBE_RHS]

                print('\t'.join(map(str, fields_left)))
                print('\t'.join(map(str, fields_right)))


class ProbePair():
    def __init__(self, lhs, rhs, pos):
        self.lhs = lhs
        self.rhs = rhs
        # coordinate of the start of the lhs on the reference
        self.lhs_start = pos

def gen_probe_pairs(seq):

    for i in range(len(seq)):
        subseq = seq[i:i+(PROBE_LEN * 2)]
        if len(subseq) < PROBE_LEN * 2: break

        lhs = subseq[:PROBE_LEN]
        rhs = subseq[-PROBE_LEN:]
                
        yield (i, lhs, rhs)

def _check_gc(seq):
    '''
    min / max come from 10x spec sheet

    >>> seq1 = 'GGGGGGGGG'
    >>> _check_gc(seq1)
    False

    >>> seq2 = 'GGGGGGGGGGAAAAAAAAAAA'
    >>> _check_gc(seq2)
    True

    >>> seq3 = 'NNNNN'
    >>> _check_gc(seq3)
    False
    '''
    ct = Counter(seq)
    if not all([i in ALLOWED_NUCS for i in ct.keys()]):
        return False

    gc = (ct['G'] + ct['C']) / ct.total()
    return gc >= GC_MIN and gc <= GC_MAX

def _check_homopolymer(seq):
    '''
    use rle to check homopolymer stretches.
    return: False if any stretches greater then max_len

    >>> seq1 = 'AAAATTTTCCCCGGGG'
    >>> _check_homopolymer(seq1)
    False

    >>> seq2 = 'ATGCTAGCTAGTCGAT'
    >>> _check_homopolymer(seq2)
    True
    '''
    hps = rle.encode(seq)
    return all([i <= HOMOP_MAX for i in hps[1]])

def _check_lhs_seq(seq):
    '''
    At the ligation junction (where the two probes
    meet), the 3' end of the LHS probe should be a
    thymine.

    >>> seq1 = 'ATGTT'
    >>> _check_lhs_seq(seq1)
    True

    >>> seq2 = 'ATGTC'
    >>> _check_lhs_seq(seq2)
    False
    '''
    if seq[-1] == 'T':
        return True
    return False

if __name__ == '__main__':
    main()


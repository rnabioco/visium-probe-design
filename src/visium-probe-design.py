#GGGGGGGGG! /usr/bin/env python3

'''
visium--probe-design.py

Design a set of ligation-based probes for the 10x Genomics Visium fixed RNA assay.

Input is a FASTA of mRNA sequence targets.

A PDF from 10x with design guidelines is in `/doc`.
'''

__author__ = 'Jay Hesselberth <jay.hesselberth@cuanschutz.edu>'
__licence__ = 'MIT'

import pdb

from collections import Counter, defaultdict

import click
import rle

from Bio import SeqIO

PROBE_LHS = 'CCTTGGCACCCGAGAATTCCA'
PROBE_RHS = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
MOD_RHS   = '/5Phos/-'

PROBE_LEN = 25
GC_MAX    = 0.73
GC_MIN    = 0.44
HOMOP_MAX = 3

ALLOWED_NUCS = set(['A','G','T','C']) 

@click.command()
@click.argument('target_fasta')
def main(target_fasta):
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

    keep_probes = defaultdict(list)
    for id, probes in cand_probes.items():

        last_pos = None
        for probe in probes:
            if not last_pos:
                last_pos = probe.lhs_start
                continue

            if abs(probe.lhs_start - last_pos) < (PROBE_LEN * 2): continue

            keep_probes[id].append(probe)
            last_pos = probe.lhs_start

    colnames = ['#id', 'pos', 'hyb_lhs', 'hyb_rhs', 'probe_lhs', 'probe_rhs']
    print('\t'.join(colnames))

    for id, probes in keep_probes.items():
        for probe in probes:
            fields = [id, probe.lhs_start,
                      probe.lhs, probe.rhs,
                      PROBE_LHS + probe.lhs,
                      MOD_RHS + probe.rhs + PROBE_RHS]        
            print('\t'.join(map(str, fields)))

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
    return: True if any stretches greater then max_len

    >>> seq1 = 'AAAATTTTCCCCGGGG'
    >>> _check_homopolymer(seq1)
    False

    >>> seq2 = 'ATGCTAGCTAGTCGAT'
    >>> _check_homopolymer(seq2)
    True
    '''
    hps = rle.encode(seq)
    return not any([i > HOMOP_MAX for i in hps[1]])

def _check_lhs_seq(seq):
    '''
    At the ligation junction (where the two probes
    meet), the 3' end of the LHS probe should be a
    thymine.

    >>> seq1 = 'ATGTT'
    >>> _check_lhs_seq(seq1)
    True

    >>> seq1 = 'ATGTC'
    >>> _check_lhs_seq(seq1)
    False
    '''
    if seq[-1] == 'T':
        return True
    return False

if __name__ == '__main__':
    main()


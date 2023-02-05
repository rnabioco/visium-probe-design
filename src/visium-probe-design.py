#! /usr/bin/env python3

'''
visium--probe-design.py

Design a set of ligation-based probes for the 10x Genomics Visium fixed RNA assay.

Input is a FASTA of mRNA sequence targets.

A PDF from 10x with design guidelines is in `/doc`.
'''

__author__ = 'Jay Hesselberth <jay.hesselberth@cuanschutz.edu>'
__licence__ = 'MIT'

from collections import Counter, defaultdict
from itertools import combinations

import click
import rle

import pdb

from Bio import SeqIO

PROBE_LHS = 'CCTTGGCACCCGAGAATTCCA'
PROBE_RHS = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
MOD_RHS   = '/5Phos/-'

PROBE_LEN = 25
GC_MAX    = 73
GC_MIN    = 44
HOMOP_MAX = 4

@click.command()
@click.argument('target_fasta')
def main():
    cand_probes = defaultdict(list)

    for record in SeqIO.parse(open(target_fasta)):
        for pos, target_seq in gen_subseq(record.seq):

            lhs = target_seq[:PROBE_LEN]
            rhs = target_seq[-PROBE_LEN:]
                
            if not _check_lhs_seq(lhs):
                continue
            if not all(map(_check_gc, [lhs, rhs])):
                continue
            if not all(map(_check_homopolymer, [lhs, rhs])):
                continue

            cand_pairs[record.id].append(ProbePair(lhs, rhs, pos))

    for id, probes in cand_probes.items():
        last_pos = None
        keep_probes = defaultdict(list)

        for probe in probes:
            if not last_pos:
                last_pos = probe.lhs_start
                continue

            if abs(probe.lhs_start - last_pos) < (PROBE_LEN * 2): continue

    colnames = ['#id', 'hyb_lhs', 'hyb_rhs', 'probe_lhs', 'probe_rhs']
    print('\t'.join(colnames))

    for id, probes in keep_probes:
        for lhs, rhs in probes:
            fields = [id, lhs, rhs,
                      PROBE_LHS + lhs,
                      MOD_RHS + rhs + PROBE_RHS]        
            print('\t'.join(fields))

class ProbePair():
    def __init__(lhs, rhs, pos):
        self.lhs = lhs
        self.rhs = rhs
        # coordinate of the start of the lhs on the reference
        self.lhs_start = pos

def gen_subseq(seq):
    for i in len(seq):
        subseq = seq[i:i+(PROBE_LEN * 2)]
        if len(subseq) < PROBE_LEN * 2: break

        yield (i, subseq)

def _check_gc(seq):
    ''' min / max come from 10x spec sheet '''
    ct = Counter(seq)
    gc = ct['G'] + ct['C'] / ct.total()
    return gc >= GC_MIN and gc <= GC_MAX

def _check_homopolymer(seq):
    ''' use rle to check homopolymer stretches.
        return: True if any stretches greater then max_len
    '''
    hps = rle.encode(seq)
    return any([i > HOMOP_MAX for i in hps[1]])

def _check_lhs_seq(seq):
    ''' At the ligation junction (where the two probes
    meet), the 3' end of the LHS probe should be a
    thymine.'''
    if seq[-1] == 'T':
        return True
    return False

def check_overlaps(seq):
    pass

if __name__ == '__main__':
    main()

# Design probe pairs for the Visium Fixed RNA assay

Input is a FASTA file containing target sequences. Assording to the
design guidelines, these ideally should be coding regions, not UTRs.

In addition, input sequences should be pre-screened for repeat elements.

To generate a sample list of probes:

```bash
python src/visium-probe-design.py test/test.fa
```

## Generating target sequences

The UCSC Genome Browser can to used to generate target sequences that
conform to the 10x Genomics specifications.

1. Navigate to the locus of interest and click an isoform.
1. Under the section "Sequence and Links to Tools and Databases", click
"Genomic Sequence at the top left.
1. Uncheck the 5/3 UTR, Introns, and Up/Downatream buttons.
1. Check "CDS Exons"
1. Check the "Mask repeats" button at the bottom, selecting "to N".

Click submit. This will generate a FASTA sequence that can be pasted into 
a list of mRNA targets.

### Checking probe sequences

Once targets are generated, you can inspect the probes using the UCSC
browser. 

Generate the probe set using `--output_fasta`. Paste the resulting probe
(which only contains the hybridization portions) into the UCSC
Blat tool, and inspect the loci for probes. The probes should only
hybridize to coding exonic regions.

## Ordering probes

The `probe_lhs` and `probe_rhs` pairs can be ordered from e.g. IDT. The
sequences specific the 5'-phosphate modification on the RHS required for 
probe ligation.

## Testing

To run the tests:

```bash
python -m doctest src/visium-design-probes.py
```

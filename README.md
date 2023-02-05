# visium-design-probes.py

Design probe pairs for the Visium Fixed RNA assay.

Input is a FASTA file containing target sequences. Assording to the
design guidelines, these ideally should be coding regions, not UTRs.

In addition, input sequences should be pre-screened for repeat elements.

To generate a sample list of probes:

```bash
python src/visium-probe-design.py test/test.fa
```

To run the tests:

```bash
python -m doctest src/visium-design-probes.py
```

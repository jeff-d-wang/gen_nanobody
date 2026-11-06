# Test fixture: 1MEL (cAb-Lys3 VHH + hen lysozyme)

Smoke-test input for `run_mpnn_cdr.sh`, so the pipeline runs before C's real targets land.
Chains relabelled from the RCSB entry: VHH `A` becomes `H`, lysozyme `L` becomes `T`, with
the second copy and all heteroatoms dropped. CDR spans are carried as
`REMARK PDBinfo-LABEL` lines, the same way `chothia2HLT.py` writes them.

```
src/evaluate/run_mpnn_cdr.sh src/evaluate/testdata/hlt
```

40 designed positions out of 132.

The CDR spans here are anchor-derived, not Chothia. They were read off the conserved
cysteine pair, the conserved tryptophan, and the WGQG motif, which gives a fixture without
needing an install. Real targets must come from `chothia2HLT.py`, which applies Chothia
H1 26-32, H2 52-56, H3 95-102, the definition A's RFdiffusion uses. Do not quote recovery
numbers from this fixture in the report.

## What a good run looks like

From the Nov 4 run at seed 37, 16 sequences:

- No framework drift in any sample. Every position outside the CDRs matches the native.
- Our pooled CDR recovery matches ProteinMPNN's own `seq_recovery` exactly on all 16
  samples, which is the cross-check that the indexing is right.
- cAb-Lys3 has an interloop disulfide, so CDR1 and CDR3 each hold one cysteine.
  ProteinMPNN kept both in 16 of 16 designs.

If CDR3 recovery comes back near 1.0 on real targets, suspect the mask rather than the
model.

# Strand-aware hyperediting competition

Introduced in development version `1.1.1-dev.10`; not a release or qualification
claim. This contract changes no model, score constant, trigger density or
primary coverage threshold.

## Coordinate convention

The input read has one original orientation. An alignment candidate on the
minus strand is scored using the reverse complement of that input. Therefore
the same conversion class has complementary names in the two candidates'
reference-oriented coordinates:

| Signature candidate | Competing candidate | Effective direction |
| --- | --- | --- |
| same alignment strand | same | retain A>G or T>C |
| opposite alignment strand | opposite | complement A>G to T>C or T>C to A>G |

Continue choosing the first existing signature-producing candidate. Preserve
that candidate's alignment strand together with its direction. Derive the
effective direction for every candidate before both genomic re-extension and
common rescoring. A competing candidate need not independently reach the
signature threshold to be scored in the same original-read alphabet.

This is not an unrestricted per-locus choice of whichever exemption wins:
two candidates on the same strand must use the same direction. Candidates
with evidence for incompatible original-read domains still use the existing
signature selection rule; resolving that ambiguity is outside this fix.

## Reproduction and scope

Synthetic two-locus fixtures have a strict leader with seven directional
substitutions on 150 bases and a reverse-oriented competitor with nine
substitutions plus one ordinary mismatch. Only the latter reaches the
unchanged 5% signature trigger. Before this correction, its reference-oriented
exemption is incorrectly applied to the other orientation, scoring 276 versus
258. The strand-consistent scores are 276 versus 286.

Four tests in `crates/map/src/align.rs` exercise:

1. Both trigger strands and both directions, without genomic re-extension.
2. The same four cases including genomic re-extension.
3. A clipped 129M21S genomic competitor recovers 150M only with complemented
   re-extension; reverting that call alone makes this test fail.
4. A same-strand negative control: a different directional spectrum cannot
   select its own exemption.

The two opposite-strand tests fail on the previous implementation before the
behavior patch. Fixtures contain generated sequences only, no human reads.

Competitor-supported human development paths motivated the investigation but
do not enter mapper decisions. New complete mapping, independent scoring,
body/ALT/tail deltas and confidence review belong to the consuming research
workspace. Passing these mechanism tests does not establish locus uniqueness,
calibrated MAPQ, or independent generalization.

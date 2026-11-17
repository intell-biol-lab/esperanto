# Experimental library-supported short-exon candidates

`map --experimental-short-exon-candidates` is disabled by default, including
through the older flow APIs. It supplements candidate discovery; it does not
change the ordinary chain gate, sequence scores, MAPQ formula, model or pair
anchor requirements.
In two-pass mapping it runs only in the final pass; discovery retains its
original candidate policy.

A clipped read can contain a 39-base exon and a 12-base exon, with a longer
contiguous remote homolog. The short exon has no 15-base seed, while the
longer exon's minimizers can span only 37 bases, below the ordinary 40-point
chain gate. Lowering that gate alone is insufficient: the existing tail
rescue postprocesses only the best ordinary extension, which can still be the
remote homolog. A synthetic regression reproduces both the original remote
winner and recovery of the library-supported splice.

When all ordinary genomic placements cover less than 90% of the query,
the option considers additional chains scoring at least 25 and below the
ordinary chain threshold. A chain endpoint must be within the existing
split-library window of a library junction. The primary index and its
ordinary occurrence cap remain unchanged. Complete seed lists share a
20,000-anchor limit; an actual anchor truncation propagates the existing
search-limited MAPQ0 rule, retaining earlier limitation evidence.

Up to eight supplemental chains receive a separate bounded continuous
extension budget. Each proposal must pass the existing split-rescue sequence
gain and strict score improvement requirements, and every emitted junction
must occur in the current library. The library may include both annotation
and junctions learned during discovery; library presence is not truth.
Library search uses reference-coordinate tail direction on either read
orientation, preserving the junction's transcript strand. Any existing clip
at the opposite outer end is retained in the full-query CIGAR; the original
retry coverage and seam requirements apply to the body plus proposed tail.
Unspliced failed proposals are discarded. The supplemental pass does not
consume the original extension budget or remove original candidates.

Accepted splice paths join the unchanged common rescoring, eligibility,
deduplication and ambiguity competition. There is no genotype lookup, target
gene exception, forced coordinate, mate-based relocation or confidence floor.
Annotation can be wrong, homologous paths may remain unresolved, and bounded
candidate generation is not exhaustive. Synthetic tests establish the
mechanism, not human read origin, confidence calibration or general mapper
superiority. Native diagnostic replay, independent interface audits and
per-read adverse accounting remain required before adoption.

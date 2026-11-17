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

After common strict rescoring, when all positive ordinary placements from both
the genomic and L1 layers cover less than 90% of the informative query,
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

## Rejected first native candidate

Candidate `6b08c75972fd8bb385b5a3956a91fce86f3823b8` passed the software
gate and recovered the intended synthetic splice. In a post hoc HG001 replay,
252 of 261 selected hotspot cases placed the same source base at STAR's
RPS27 coordinate, and 248 matched STAR's full supporting alignment. Both
281-pair BAMs passed the independent full-sequence/quality/SAM contract.
This is a mechanism result, not per-read origin truth. One previously mapped
read became unmapped; the small replay differs from the full run's junction
discovery context for that template.

Seven reused simulated development cells (872,459 pairs) rejected general
adoption: 59 ancestral-exact reads gained versus 19 lost; 5,493 correct bases
gained versus 451 lost; 17 correctly projected RNA-event bases gained versus
11 lost; 16 new false-positive junctions and 16 new false-positive deletions.
Net improvements cannot replace these gross loss gates. Every input remains
in the independent scorer denominator. The option remains disabled by default
and must not be described as a qualified improvement.

One observed adverse mechanism is terminal mismatch clipping: a supplementary
partial splice can outscore an existing full transcript path after discarding
an edited terminal base. A new discovery proposal must account for the
existing transcript and genomic explanations together, with adverse
source-base and event projection controls. Simply lowering the chain gate,
forcing the mate locus, or changing the downstream evaluator would not resolve
this failure.

Private evidence is retained outside Git under
`experiments/esp-short-exon-library-20261117-attempt01` in the research
workspace, including every changed/adverse score row and the rejected gate.

## Coverage qualification repair

The first supplement checked `informative_query_coverage` before genomic
placements had been rescored, when `strict.aligned` was still zero. It also
omitted L1 candidates from that check. The supplement now runs after the
ordinary candidates from both layers have received common strict scores,
before coverage-class deduplication. A positive path covering at least 90%
of the read suppresses supplemental discovery. N columns do not provide
informative coverage. The genomic minimizers and original search-limitation
evidence are retained; extra candidates still receive identical rescoring.

Generated regressions reproduce replacement of an existing high-coverage
genomic path and loss of a terminal A-to-G base from a complete L1 projection
on both query strands. This is a discovery-eligibility repair, not a change
to match, mismatch, clipping, gap or splice scores. False tail extensions and
pair effects on incomplete pools require native all-input evaluation.

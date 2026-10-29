# Known scientific qualification issues

Dev.10 corrects an orientation error in hyperediting re-arbitration: A>G/T>C
must complement for candidates on opposite alignment strands to represent one
original-read conversion domain. This does not resolve conflicting signature
domains, same-base terminal splice ambiguity, short-body confidence, or
independent held-out qualification. See the
[strand-aware editing-domain contract](strand-aware-edit-domain-contract.md).

Optional `map --candidate-diagnostics` retains legal positive candidates below
primary coverage for inspection. This adds observability, not recall recovery:
96/98 bp candidates in 150 bp reads remain below the unchanged 99 bp minimum.
The sidecar is captured before collapsed rescue and must be joined with final
BAM for final placement status. See
[the executable diagnostics contract](subthreshold-candidate-diagnostics.md).

## Internal soft clipping in primary CIGARs

Status at organization adoption (2026-10-29): **unfixed**. The
`fix/internal-softclip-cigar` development candidate implements a repair;
integration and scientific qualification must be assessed separately.

The inherited mapper emitted primary alignments with internal soft clipping in
the human small development pilot:

| Seed | Invalid primary records |
| --- | ---: |
| 20261028 | 51 |
| 20261029 | 35 |

Representative CIGARs:

```text
10M1S44595N139M
133M5S26110N10M2S
5M1D4M1D8M4S11643N5S124M
```

Frozen ESP and research rescue-on/off variants all failed strict SAM interface
qualification. Synthetic regressions now demonstrate that split rescue joins
two independently local extensions without anchoring their inner ends; the
same problem affects run-through tail extension. A separate caller defect
retains the old start after right-tail reextension. The candidate checks the
seams before accepting candidates and updates that start. See the executable
[rescue contract](cigar-rescue-contract.md).

The research workspace retains the original BAM records, FASTQ sequence/quality,
truth provenance and independent evaluation results. Those human read fixtures
are not imported into this public foundation. A future regression fixture must
be synthetic or have explicitly verified redistribution rights.

A fix must preserve query/reference consumption and correct genomic placement,
not merely replace internal S with I or drop failing records without accounting.
Add focused unit tests, retain source truth / record denominators, rerun both
seeds under strict evaluation and review changed valid alignments as well.

The publication source commit is
`d37fda36c16a78d0c403e0b287bcae3a14702192`. The fork starting commit
`3d3c3adfae57f8f8164d1dc5dd5f988f8a18c865` has the same mapping source.
Foundational formatting and repository changes do not resolve this defect.
The repair is a conservative candidate-eligibility change, not an anchored
dynamic-programming implementation: it can reject a local optimum even when a
different boundary-constrained alignment could be valid. Unit tests establish
the construction behavior; named human reruns and independent scoring establish
only the qualification and effects observed on those development inputs.

## Rescue arbitration regressions found during dev.2 review

Full-seed, sequence-filtered diagnostics identified a weak novel-tail proposal
that seeded a second-pass junction prior, a perfect contiguous read rewritten
solely by that prior, and a repeat-derived segmented chain whose successful
splice suppressed the correct continuous alternative. The dev.3 candidate
adds strict proposal-versus-primary acceptance with the existing rewrite cost,
requires actual sequence evidence gain, and generates continuous segment
alternatives alongside splices. Pure synthetic regressions fail on dev.2
before this change. No human read data are included in the public repository.

The research workspace retains exact instrumented patches and full-seed
receipts; instrumentation was checked against every uninstrumented primary.
Dev.3 qualification, residual regressions and adoption must be assessed on its
own named outputs, independently of the dev.2 interface result. This remains
development software; no general mapping-accuracy or biological claim follows
from these mechanism tests.

## Coverage eligibility and fallback found during dev.3 review

Added continuous competitors exposed coverage refusal after destructive
same-locus dedup: a shorter higher score could remove a full candidate, and
rejecting the score leader ended arbitration without considering a fallback.
The dev.4 candidate separates coverage classes during dedup, selects the first
positive full-coverage candidate, and retains different-locus short evidence
in the MAPQ margin. Reference deletions do not count as query coverage.
Full nonpositive candidates remain available for existing editing rescoring,
but cannot become primary without qualifying in the active domain.

The executable [candidate contract](candidate-eligibility-contract.md) and
eight synthetic mechanism tests define these rules. They do not require all
previously mapped low-confidence reads to be restored. Qualification,
per-read body/ALT changes and any remaining wrong placements require a new
complete human development replay and independent review.

## Fast intron proposals suppressing generic candidates (dev.4 review)

Ten refused human pilot reads exposed candidate loss before arbitration. The
fast intron path searched all anchors, then its caller rewrote only N length
using another chain's junction and returned early. Common scoring rejected
that bad explanation, but an independently supported generic splice candidate
was never offered. Dev.5 removes this geometry rewrite, validates CIGAR/declared
junction agreement and continues generic generation. The independent scorer
and coverage threshold remain unchanged.

A generic candidate can still cover only the last two of three true exons:
the missing first exon lacks chain anchors in this case. Recovery of that
partial body does not establish full read recall, general accuracy, or MAPQ
calibration. Human replay results belong in the separate research report.

Dev.5 alone is not a freeze candidate: seam rejection exposed that the fast
constructor already combined refined genomic coordinates with a midpoint
query split. Some supported fast proposals disappeared when the generic
chain could not reproduce them. Dev.6 reconstructs from the selected fragment
pair using pseudo-reference extension and consistent stitching. Historical
dev.5 output remains a failed scientific qualification candidate even though
its interface and unit gates passed.

## Neutral endpoint and read/transcript orientation competition (dev.6 review)

The `catalog_10@512` regression included three clipped body bases whose strict
score was exactly zero (two matches, one T>C). Local traceback/end trimming
removed them on a score tie. Right-end clipped ties at identical coordinates
could also win same-locus dedup by generator order. Dev.7 retains neutral
informative diagonal flanks and prefers informative coverage after score and
ordinary mismatch ties. Negative clips and ambiguous columns remain excluded.

Read-aware junction refinement previously consulted only the read alignment
strand, while common rescoring accepted annotation/motifs on either RNA
strand. This could slide the opposite mate 1-2 bases at repeat boundaries.
Dev.7 considers both transcript strands in that same refinement search.
These synthetic mechanism fixes require independent full-seed/event review;
they do not establish global recall or calibrated MAPQ. Historical dev.6
NO-GO and its lost event remain preserved research evidence.

## Non-N technical tails and endpoint tie consistency (dev.7 review)

Dev.7's zero-score flank scan admitted ordinary mismatches, unlike common
arbitration's preference for fewer ordinary mismatches before coverage. This
could absorb a neutral technical flank solely to add query coverage. Dev.8
stops at such mismatches, retaining only matches and A>G/T>C columns up to the
longest neutral endpoint before that barrier. No score, model or eligibility
threshold changes.

The original development review found 499 reads adding nonbody bases without
a correct-body gain; 472 had new ordinary-mismatch columns and 27 only had
matching/edit-direction columns. These are diagnostic cohorts, not a promised
count of corrected reads. Positive-score accidental technical alignments,
technical sequences identical to terminal edit evidence, incomplete mapping
and population MAPQ calibration remain unresolved. Five new synthetic tests
include explicit positive-score and identical-input limits. Full independent
replay must determine the actual per-read/body/event costs of the change.

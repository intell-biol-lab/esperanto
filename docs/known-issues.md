# Known scientific qualification issues

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

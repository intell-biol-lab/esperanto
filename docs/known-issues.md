# Known scientific qualification issues

## Internal soft clipping in primary CIGARs

Status at organization adoption (2026-10-29): **unfixed**.

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
qualification. These observations do not prove which mapper construction path
is responsible. Trace the original read, block alignment and splice projection
before deciding the repair.

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

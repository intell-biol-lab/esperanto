# map — splice-aware read alignment

Dev.19 retains full-query exact structural alternatives at one locus, displays
a deterministic representative and preserves MAPQ0 ambiguity regardless of
annotation rewards. See [the exact-path contract](../../exact-path-competition.md).
`map --no-exact-path-competition` independently restores dev18 policy.

Dev.18 tests additional context proposals against a constant-base background
in the same strict score domain. `map --no-context-background` independently
restores the dev17 policy. See [the contract](../../context-background-contract.md).
This changes proposal qualification only; human adoption remains separate.

The mapper aligns single-end or paired-end FASTQ reads to a `.paidx` genomic index and writes BAM alignments. RNA alignment uses minimizer seeds, chaining, splice-aware extension, and paired-end evidence.

Dev.17 adds bounded ordinary high-frequency context proposals after failure of
both original chain rounds, while preserving old edge/splice/indel candidates.
The additional scan/proposal cost and conservative MAPQ0 policy are specified
in [the repeat-context contract](../../repeat-context-contract.md).
`map --no-repeat-context` is a component ablation; converted rescue disables it.
Human qualification and scientific adoption remain separate.

In dev.16, an ordinary RNA query with no chain from either initial round can
complete its query-edge seeds within the existing mid-round budgets. Successful
rounds, the reference sketch, DNA and editing-variant paths are unchanged.
See [the bounded query-edge contract](../../query-edge-seeds.md).

## Inputs

- Plain or gzip-compressed FASTQ with synchronized paired-end records.
- A `.paidx` built from the matching reference; default minimizer parameters are k=15 and w=5.
- An optional GTF junction library.
- An optional complete transcriptome L1 bundle and `.tidx` sidecar built from the same reference.

`AlignConfig::rna_default()` supplies the production parameters. Primary alignment uses strict four-letter matching. Junction discovery and a second alignment pass support unannotated splice junctions. Annotation and discovered junctions are combined deterministically.

The paired-endpoint development proposal infers conservative read-through
intervals before alignment/cache interning. RNA defaults enable it; DNA and
single reads are unchanged. `map --no-fragment-endpoints` disables it.
Supported effective intervals use the existing 0.66 fraction while preserving
raw sequence/quality/output coordinates and declaring `XF:i:T` evidence.
This is a new conditional eligibility/preprocessing policy, not scientific
adoption; see [the endpoint contract](../../paired-fragment-endpoints.md).

## Coordinates and orientation

Alignment positions and reference intervals use 0-based coordinates. CIGAR operations are `M`, `I`, `D`, `N`, and `S`; skipped introns use `N`. BAM sequence for a mapped minus-strand record is reverse-complemented and its qualities are reversed. Unmapped reads retain their original orientation.

## Placement evidence

Candidate placement ranking considers alignment score, alternative locations, splice structure, and paired-end agreement. MAPQ is bounded by 60. Ambiguous placements remain representable with low MAPQ; downstream observations require MAPQ at least 20.

The aligner can recover a read using an anchored mate. Alignment rescue and its `RE` tags are read-level evidence and do not bypass the downstream model's filters. Decoy reference contigs use the `decoy_` prefix and compete during alignment.

## Collapsed-read rescue

When `.ag.paidx` and `.tc.paidx` are present beside the primary index, the pipeline can realign unmapped and insufficiently covered reads using separate A-to-G and T-to-C transformations. These indexes use k=31 and w=10 and exclude decoy contigs.

Placements must pass verification against the true four-letter reference. Required evidence includes an ungapped alignment, at most two non-editing mismatches, editing density at least 5%, editing-class purity above 60% (above 80% for reads of at most 60 bases), high-quality edited bases, and internal cluster geometry. Paired-end placements additionally require mate concordance. Accepted records carry `RE:Z:collapsed`; rescue outputs are sealed only after the merged BAM is complete.

Valid XF-bound primary records, including unmapped ones, do not re-enter
raw-read rescue or stale-rescue replacement. The complete original unmapped
pool is preserved except for actual rescued identities, including after resume.

## Outputs

`raw.bam` is unsorted alignment output. `unmapped.fq.gz` contains unresolved reads. `align_qc.json` records total and mapped reads, proper pairs, insert-size summaries, MAPQ histogram, junction counts, read-rescue counts, and elapsed time. `discovered_junctions.bed` reports the junction-discovery result. The flow crate sorts and indexes BAM output for discovery and scoring.

## Index format

The `.paidx` magic is `PAIDXfmt` and the format version is 1. The file stores seed parameters, reference sequences and contig metadata, minimizer keys, and positions. Loading validates the file structure and integrity information before alignment.

## Errors and determinism

Malformed FASTQ, inconsistent mate records, invalid indexes, and I/O errors are reported with input context. Placement decisions and record order are deterministic. Worker scheduling does not choose between tied placements or alter scientific result ordering.

# Full-query exact structural alternatives

Dev.19 RNA alignment retains distinct full-query exact mappings that overlap
at one locus instead of allowing annotation priors to erase them. This is a
development mechanism; software CI does not establish scientific adoption.

A placement is exact-full only when all query bases match unambiguous reference
bases, no substitutions/indels/clips occur, and its CIGAR contains positive M/N
operations. Adjacent M operations are coalesced for mapping identity. An L1 and
genomic copy of the same mapping remains a duplicate.

Ordinary scoring and its per-junction priors are unchanged. If the ordinary
eligible leader is exact-full and another exact-full mapping overlaps its
reference M footprint by at least 0.5 on the same contig/strand, both remain
observable. Among that leader's compatible local paths, primary representation
uses fewer introns, then the existing genomic/L1 origin rank, position and
normalized mapping. This deterministic choice is not proof of the source
isoform. It is limited to exact-full, sequence-indistinguishable paths;
imperfect and partially aligned paths keep the previous arbitration.

The selected placement retains its actual common score, but receives durable
`exact_path_ambiguous` evidence and MAPQ0. Cache/mate re-arbitration recomputes
the evidence; explicit confidence overrides cannot defeat it. EK reports
`mapq_src=exact_path_tie` and a tied policy margin (`chain=sub`, `dq=0`).
These margin fields represent sequence equivalence, not two raw annotation
scores. Raw scores remain available in candidate diagnostics. Existing limited
seed-search evidence takes precedence and still forces MAPQ0.

This does not delete short exons. A 1/2bp terminal exon that supplies a
discriminating reference base can retain its supported spliced placement.
Sequence-identical true short-exon and continuous sources cannot be resolved
from this read alone; both receive the same ambiguity treatment.

`map --no-exact-path-competition` disables retention, representative selection
and ambiguity policy, restoring dev18 at the same head. RNA defaults enable
the policy; DNA defaults disable it. Prior flow APIs remain wrappers around
`map_stage_with_exact_path_options`.

No new candidate-generation, DP or posting budget is added. Additional retained
alternatives were already generated and strictly scored. There is no complete
search or globally calibrated probability claim. The 0.66 coverage rule,
read-through intervals, repeat-context backgrounds, collapsed rescue and
model weights remain separate mechanisms.

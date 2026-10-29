# Bounded ordinary repeat-context proposals (dev.17)

This is an unadopted development proposal, not a qualified release or a claim
of mapper superiority. It uses the unchanged ordinary index and scoring model.

## Admission and work

Ordinary RNA may attempt context verification only when both original chain
rounds have no chain. Preserve query-edge completion and all its existing
splice/indel/continuous candidates. A perfect whole-query ordinary candidate
already reaches the score ceiling and bypasses the additional search. DNA,
editing extension variants, queries above256bases and converted rescue do not
use this proposal. RNA defaults enable it; `map --no-repeat-context` or
`AlignConfig.repeat_context=false` independently disables it.

Query all valid canonical k-mers at index k. For lists above5000occurrences,
sort by(count,query position,canonical kmer,strand) and admit only whole lists
that fit an ADDITIONAL20000posting scan bound. Oversized lists are skipped, not
sampled by coordinates. Different query starts remain different hypotheses.
Reference intervals must contain the entire query; do not pad outside contigs.

Deduplicate full-query diagonals by(contig,start,strand), then find a strictly
scored local ungapped interval containing the queried exact seed. N scores0
and does not buy coverage. Equal-score choices prefer fewer ordinary mismatches,
then more informative bases: neutral edit/match blocks may remain, ordinary
neutral tails cannot buy coverage. Emit legal M/S proposals, with no fictitious
chain score. Retain at most8additional proposals, prioritizing the existing
primary-coverage class before the usual score/mm/coverage order.

The candidate origin is explicitly `repeat_context`,not chain-derived
`genomic`. It is common-rescored if an editing domain is selected but does not
enter genomic DP reextension. This preserves the additional-family work bound
and avoids inventing chain evidence; stronger existing genomic/L1 candidates
still win exact common-order ties.

The old8chain/two-DP scheduler is unchanged; these8ungapped proposals and up to
20000*256verified diagonal columns are EXTRA work. Common strict rescoring,
coverage-class dedup and final arbitration consider old and new proposals
together. This is not an ungapped proxy used to delete splice/indel paths.

## Incomplete confidence and durable evidence

An optional `seed_budget.repeat_context` object records query length,
eligible/scanned/omitted seed lists,eligible/scanned postings,unique intervals,
and positive/retained/omitted proposals. The old three seed-budget anchor
counters retain their original meaning; they may all be0for a context-only
search. Absent context data preserves the prior serialized object.

Every attempted context search is conservatively limited in this development
version, regardless of whether its selected lists fit: ungapped contexts do not
establish globally complete indel/splice search or calibrated MAPQ. The existing
limited flag forces MAPQ0for any winner,including old-edge/L1/mate winners.
Empty candidate caches carry it as well. Existing EK/raw-BAM,rescue and restart
policies apply; a later override cannot restore confidence. Diagnostic pools
copy the same typed budget. XF endpoints still exclude raw adapter rescue.

## Controls and limitations

Native-index generated witnesses cover sparse AG and ordinary transversions on
both strands; identical copies,whole-list overflow,empty-cache mate fallback,
other-mode/ablation/perfect-edge bypass,full cached-candidate parity,and final
one/two-pass paired BAM raw sequence/quality/XF/diagnostic consistency are tested.
The local kernel is checked against exhaustive subinterval score/mm/coverage
ordering on generated short diagonals. Positive-score technical alignment
remains an explicit indistinguishability limit,not a successful negative.

Synthetic success does not establish human qualification or superiority.
Freeze a new exact-source candidate and full replay before adoption. Existing
terminal ALT losses,technical absorption,independent source/region/family
holdout,equal-information competitors,MAPQ calibration and resource comparison
remain separate requirements.

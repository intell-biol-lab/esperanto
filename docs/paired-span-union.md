# Paired transcript span and shared introns

Dev.21 measures the pair's enclosing reference interval and removes the union
of both CIGARs' `N` intervals. Shared,partially overlapping,nested or touching
intervals contribute their total covered length once. `M` and `D` consume
reference; `I` and `S` do not. Mate order does not change the magnitude.

For `50M2000N30M` at100 and `30M2000N50M` at120,the enclosing interval is
[100,2200),the shared skip is[150,2150),and the transcript span is100.
Previously the same skip was summed twice,returning zero after saturation.
That could understate insert summaries and incorrectly admit an excessive
span as proper on generated controls.

The same-contig/opposite-strand rule and1000-base inclusive proper span limit
remain. Paired record TLEN keeps ESP's inherited transcript-aware convention;
this document does not redefine standard SAM TLEN semantics. Proper pairs
receive opposite signed values,while the pipeline continues emitting0 for
discordant/unplaced pairs. Sign handling and mate coordinate semantics are
unchanged.

This fixes arithmetic,not candidate search,score,MAPQ or coverage. Synthetic
geometry tests include an independent small-coordinate mask and native
one/two-pass spliced BAM with both strand orders and nonuniform quality.
Scientific qualification,full replay and independent generalization remain
separate from software checks.

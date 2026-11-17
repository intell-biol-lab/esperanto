# Ordinary context background competition (dev.18)

Additional high-frequency `repeat_context` proposals must strictly beat the
best local score of the entire effective query against any constant A,C,G,T
reference. Use the existing strict match/mismatch and N=0 domain. A tie rejects
the additional proposal. Compute the maximum over all query subintervals for
each base in O(4n), once per attempted context search.

This is a candidate specificity policy, not a calibrated statistical test or
proof of biological/technical origin. A whole read is not blacklisted for its
base composition. Stronger full low-complexity genomic matches can survive;
identical copies remain ambiguous. A weak real local match may also lose to
the background, which must be assessed in the independent complete replay.

No ordinary chain/edge/L1/splice/indel/mate candidate is deleted by this gate.
Only the additional family's proposals are filtered, before its8candidate
capacity selection. Its posting/DP bounds,common scores,0.66coverage policy
and conservative limited-search MAPQ0 remain unchanged. The optional context
budget and an empty candidate cache remain durable after every rejection.
Converted rescue disables this family; this gate is not a converted-domain
test.

`AlignConfig.context_background` defaults true for RNA,false for DNA.
`map --no-context-background` independently disables this gate and reproduces
the dev17 proposal policy. `--no-repeat-context` still disables the whole
additional family. Prior flow option wrappers remain available.

When enabled, `seed_budget.repeat_context.background` records:

- `score`: the strict constant-base local maximum for the effective query;
- `rejected_candidates`: positive,exact-seed-containing distinct diagonals
  whose best local scores do not exceed that maximum.

`positive_candidates` still counts emitted positive proposals before capacity
truncation: positive=retained+omitted. Background rejections have their own
counter. Repeated seeds on the same rejected diagonal are counted once.
When disabled the optional object is absent, preserving old diagnostic shape.

Generated native-index tests cover both-strand partial low-complexity refusal,
full unique and full repeated low-complexity positives,ambiguity,independent
ablation,counter/cache/mate limitations,one/two-pass paired BAM rawSEQ/QUAL/XF,
background ties and an independent exhaustive local-score oracle. These are
mechanism tests,not a claim of human qualification or mapper superiority.
The existing stronger splice/indel and sparse-edit/transversion controls remain.

Constant-base backgrounds do not model all adapters,tandem repeats or complex
artifacts. Previously observed terminal body/ALT losses and other technical
absorption remain separate blockers. Full independent development replay and
source/region/family holdout,equalQC/information comparisons,MAPQ calibration
and resource measurements are required before broader claims.

# Terminal background anchor qualification

dev.29 extends the existing ordinary repeat-context background policy. A
constant terminal run can supply almost an entire boundary-crossing seed;
a few incidental flank matches can then exceed the constant-background score.
The new qualification requires independent flank support rather than a larger
score margin.

For each effective query, identify maximal terminal A/C/G/T runs whose length
is at least the index seed length `k`. A context proposal must still strictly
beat the existing constant-base background. If either terminal run qualifies,
its selected ungapped local interval must also contain an exact, unambiguous
`k`-base anchor wholly outside both qualifying runs. Verify this anchor in the
reference interval already fetched; do not search another posting list or
extend/change the selected alignment. Without a qualifying terminal run the
existing seed-containing-interval test suffices.

The rule applies only when `context_background` is enabled. The
`--no-context-background` ablation disables both background checks. Ordinary
chain, splice, indel and mate generation, scores, coverage, model, search
budgets and MAPQ policies are unchanged. Weak-mate proposals remain opt-in.
No read is blacklisted globally: an independently supported full biological
sequence remains eligible through this or other existing candidate paths.

Diagnostic `repeat_context.background` retains the original `score` and
`rejected_candidates` (score-rejected) fields and adds:

- `terminal_prefix_bases`, `terminal_suffix_bases`: excluded lengths in the
  effective query's input orientation; internally swap on the minus strand.
- `rejected_terminal_candidates`: seed-containing, positive, score-qualified
  diagonals rejected for missing an independent flank anchor.

Count each diagonal once, before the candidate-retention cap. Positive
candidate counts include only proposals passing both checks. Rejected
proposals do not erase the existing incomplete-search evidence.

Generated native-index regressions cover the positive-score four-base flank
failure, both orientations and scan budgets, `k-1` versus `k` flanks, long
genuine biological flanks, remote identical competitors and background-off
reproduction. Unit controls cover all constant bases, variable index k,
both-ended runs, ambiguous/nonmatching anchors and clipped-away support.

This is a conservative eligibility policy. A true biological run supported by
only a short matching flank can be indistinguishable from a technical tail and
receive the same refusal. Positive-scoring technical alignments with full
independent flank support, interrupted tails, adapters and nonterminal repeats
remain possible. Neither rejection nor acceptance proves biological origin,
complete search or calibrated uniqueness. Human-development qualification,
gross losses and held-out evaluation are separate gates.

# Affine gap scoring domain

Extension and common candidate rescoring use the same positive gap parameters.
A run of `n >= 1` inserted or deleted bases costs
`gap_open + n * gap_ext`; one base therefore costs open plus extension.
Ordinary defaults are open=4 and extension=2. The existing collapsed-rescue
configuration uses 400/400. Neither parameter set changes in dev.24.

Before dev.24, extension charged open only on entry and extension on each
subsequent base. Common candidate rescoring already charged open+n*extension.
A selected extension could consequently lose one extension penalty per gap
when rescored, after an equally good or better ungapped proposal had vanished.
For example, a synthetic 120-base exact body followed by a one-base reference
deletion and three matches scored 242 in extension but 240 in common scoring.
Clipping the tail scored 240 throughout.

The shared core now charges open+extension on entry and extension on
continuation. Scalar i16/i32, AVX2 full/banded, traceback and end trimming all
consume that core. Existing tie conventions still apply: neutral informative
diagonal flanks may be retained; neutral gapped flanks receive no new tie rule.
No editing direction, coverage rule, splice prior, MAPQ policy, model or
independent evaluation semantics change. This is a consistency correction, not
an optional rescue heuristic.

Regression coverage includes central insertions/deletions of lengths 1,2,3,7;
reverse complements, offsets, Off/Ag/Tc domains, hinted/full extension and
execution-path parity; agreement with common rescoring; native true-indel
mapping; terminal-gap competition and the high-penalty rescue configuration.
Production parameter sets fit the existing packed accumulator representation.
These tests do not certify arbitrary user-constructed scoring ranges.

Scientific qualification, per-read adverse changes, confidence effects and
independent generalization require their own named replays. Passing these
mechanism tests does not establish mapper superiority or authorize adoption.

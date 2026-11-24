# Experimental repeat-context budgets

dev.28 exposes two independent `map` options:

- `--experimental-context-postings N` (default 20000, positive u32).
- `--experimental-context-candidates N` (default 8, positive u32).

The first limits scanned postings, retaining a seed's whole reference hit list
or skipping it. The second limits retained ordinary ungapped context proposals.
It does not limit all mapper candidates. Both are independent of the dev.27
`--experimental-weak-mate-proposals` flag, which remains disabled by default.

Defaults preserve dev.26 output. Existing flow entry points retain 20000/8.
The new flow wrapper accepts budgets explicitly and rejects zero values.
Limits apply each time context generation runs, including pass1; diagnostics
record actual scans and omissions, and scientific runs must compare discovered
junction libraries. Weak mate still runs only in pass2.

Changing budgets never changes the conservative context MAPQ0 policy. Even a
fully scanned set of ungapped seed-based proposals does not establish complete
search over gapped, spliced or other omitted explanations. No threshold, model,
ordinary anchor budget, extension budget, coverage rule or scoring is changed.

The development ablation predefines the 2x2x2 matrix of postings 20000/80000,
retention 8/32 and weak mate off/on. Only fresh development sources are used.
All cells require strict SAM and independent column qualification, with body
losses, ALT losses and added technical query positions recorded separately
from gains. Preserve unsuccessful cells. Default and weak-only runs must match
the frozen dev.26 and dev.27 final records respectively. Diagnostic-free timed
runs must match their diagnostic counterparts.

These experimental controls are not an optimized default or an adoption.
Do not merge, release, change model weights or use held-out results to tune them.

# Score-tied mate-window experiment

`--experimental-mate-window-score-ties` is disabled by default. Existing flow
entry points keep it disabled. It applies to paired final mapping, with the
existing confident, unrescued mate requirement; discovery keeps the parent
policy.

The additional tie case requires an already placed original read on the same
contig as the mate. The supporting mate must be unrescued, have no search-limit
flag, meet the existing pileup MAPQ gate, and match its reference base-for-base
over one full-length unspliced match operation. A score margin driven by a
splice prior or mismatches cannot enable the exception. The coordinate-only
low-level mate API keeps the old policy because it lacks this support evidence.

Previously, the first eligible candidate had to overlap the proposed mate
window before the mate's reach filter could run. A remote repeat copy sorted
first could veto an equally scoring exact full-length window candidate.

The experiment also admits a retained candidate matching the window's contig,
start, strand and CIGAR when its strict score equals the global eligible leader.
The window must consist of one full-length match operation, match every query
base, and use ordinary scoring. A lower score, terminal clipping or an editing
proposal cannot qualify through this exception. Existing reach filtering and
arbitration then run. Every newly admitted result is rescued and forced to
MAPQ0. Search limits and existing coverage requirements are preserved.

This is conditional coordinate selection using mate evidence. A confidently
misplaced mate can still select the wrong repeat copy. Nearby ties and splice
alternatives remain unresolved by this change. Unit tests establish the veto
mechanism, default preservation and refusal of a lower-scoring proposal; they
do not establish scientific gain or generalization.

An initial four-cell development pilot of the wider tie policy gained 2,149
exact reads and lost ten; it failed its predeclared zero correct-base loss gate.
The narrower support policy is a new candidate requiring a fresh experiment.
These restrictions cannot rule out every wrongly placed exact supporting mate.

Development comparison must use the previous stronger weak-mate candidate,
with gross read, pair, source and base losses, junction TP/FP and fixed-site
no-edit ALT gates. Keep models, thresholds, truth and the independent evaluator
unchanged. Preserve failed experiments. No full BAM qualification or adoption
follows from enabling the option.

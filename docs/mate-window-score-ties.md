# Score-tied mate-window experiment

`--experimental-mate-window-score-ties` is disabled by default. Existing flow
entry points keep it disabled. It applies to paired final mapping, with the
existing confident, unrescued mate requirement; discovery keeps the parent
policy.

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

Development comparison must use the previous stronger weak-mate candidate,
with gross read, pair, source and base losses, junction TP/FP and fixed-site
no-edit ALT gates. Keep models, thresholds, truth and the independent evaluator
unchanged. Preserve failed experiments. No full BAM qualification or adoption
follows from enabling the option.

# Compact molecular feature windows

The native27 reader previously built template-to-position maps and then copied
observations and full template names into position-to-template maps. High-depth
RNA regions could exceed the analysis cgroup limit. This change interns the
read-group + normalized template name into a checked u32 ID, merges directly into
position maps, and removes unresolved conflicts only after all reads have merged.
It does not cap depth, discard duplicates, alter BQ/MQ, or change mate-conflict
resolution. IDs are local identities, not hashes; unrelated templates do not
collide. IDs stay consistent across all positions in a window, preserving reuse.

Feature export groups are limited to401bp (one distinct target's +/-200 window)
to avoid connected-site window expansion. This trades repeated indexed reads for
lower peak residency. Memory still scales with the evidence in a single window;
this is NOT an absolute constant-memory guarantee. `score_sites` grouping is not
changed by this patch. Floating-point reduction order can change, requiring
NaN-aware numerical parity testing; neither tests nor a build certify biology.

Regression tests cover compact position keys, RG separation, mate suffixes,
equal-quality conflict exclusion and later higher-quality conflict resolution.
Real depth/RSS and baseline feature parity are recorded outside the repository
in the ESP research experiment; no human read data is committed here.

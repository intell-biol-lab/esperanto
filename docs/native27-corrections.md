# Native27 numerical corrections (dev.42)

This is an experimental scoring/extraction correction, not a new model or an
adopted release. Base: `c045225eb9419296600f7d45e629f77bc8b5e9c4`.

The same 27-column schema, uniform editing-fraction prior, symmetric BQ error
model, filters, template merge and model weights/threshold are retained.

## Reproductions and changes

1. The <=63-observation likelihood is a polynomial of degree <=63, so GL32
   integrates it exactly in exact arithmetic. Summing in linear space can
   underflow when many high-BQ third-base observations contribute tiny constant
   factors. Sum the weighted node likelihoods with log-add-exp instead.
2. For >63 observations the integrator splits a mode-centered interval. The
   largest derivative over *both* endpoints is not a local mass width:
   it compresses the boundary intervals while fixed quadrature misses mass in
   the remainder. Use slope at the actual boundary combined with square-root
   curvature there. Curvature matters when the boundary derivative is nearly0.
3. Evaluate `(1-p)*P(base|REF) + p*P(base|ALT)` directly. Subtracting almost equal
   numbers at p=1 loses accuracy in the reference-base endpoint when BQ is high.
   Derivative and curvature use the identical stable endpoint mixture.

Frozen failures include depth10,000/ALT3/BQ30 and depth63/third-base-only/BQ93.
The former's independently checked BF pair is
`(-8.254418713486338,79997.38405204924)` rather than
`(-12.83825901833373,79992.80021174355)`. The latter is `(0,0)`, not `(-inf,-inf)`.
Closed-form pure REF/ALT tests cover63/64/1000/10000 observations and
BQ13/20/30/60/93. Mixed-BQ and nearly flat-boundary oracle cases are included.

These are numerical fixes, not evidence that the inherited classifier is
calibrated with the corrected inputs. Real-BAM feature impact and independent
label admission must precede training or model efficacy claims. Keep old
feature datasets and binaries; do not relabel old columns as corrected.

The equal-BQ same-base template metadata order dependence is a **separate**
semantic correction; it is not addressed by dev.42. Existing low-depth
third-base edge failures are not evidence that already accepted data contains
infinite values. Fixed-node quadrature still needs independently quantified
precision/runtime over any expanded operational domain.

## Deterministic template metadata (dev.43; separate stacked change)

At equal maximum BQ, reduce auxiliary fields for same-base observations as
well as different-base conflicts. Use maximum finite MAPQ (NaN if all unknown),
maximum read-end fraction, minimum softclip/gap-nearness fractions, and
orientation consensus (`-1` when inconsistent). These operators already applied
to the conflicting-base branch; extending them to same-base ties makes the
metadata reduction commutative and associative. The result is explicitly
**composite template evidence**, not metadata of a selected physical read.

Base/BQ are unchanged: unequal-BQ observations keep the highest-BQ observation,
equal-BQ different bases conflict, and only a strictly higher BQ can resolve
an existing conflict. Counts, VAF, likelihoods, neighborhood counts, sequence
and annotations must remain unchanged relative to dev.42. Potentially changed
columns are ALT mean MAPQ/end-fraction/softclip/gap-nearness and strand bias.
Unknown orientation is not assigned to either direction by the strand statistic.

This is a semantic change, kept as a separate commit/build/PR and not claimed
to match the frozen model's training distribution. Real-BAM impact must be
reported separately from numerical corrections. No old Parquet or model
is overwritten and no model is fit or scored as part of this correction.

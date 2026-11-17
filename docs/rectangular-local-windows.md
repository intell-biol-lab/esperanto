# Rectangular full-matrix local alignment

Local alignment permits a reference window shorter than the query. Stitched
splice pseudo-references naturally have this shape when a raw read has clips.
The old full-matrix scalar and AVX2 entry points asserted `m<=n` in debug
builds, aborting otherwise valid splice pipeline controls.

The assertion is not an algorithmic precondition. In scalar full-matrix mode
the row center is clamped into [1,n], so half-width n covers all columns even
when m>n. The AVX2 sweep uses bounds for the full [1,m]×[1,n] rectangle. Remove
the two unjustified debug assertions only. Do not widen bands, change scoring,
dispatch, recurrence, traceback, tie behavior or DP budgets.

Literal short-window controls verify late-query matches for both accumulator
widths, public dispatch and AVX2. Randomized scalar/AVX2 parity and small
16-lane boundary shapes now include m>n. The pre-fix release scalar literal
control already passed; this is a debug-contract correction, NOT a claimed
release alignment recovery. Human disabled-endpoint replay must still verify
parent equivalence rather than assume it.

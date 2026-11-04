# Extension priority after limited repeat search

Development proposal: `1.1.1-dev.13`. This rule does not establish scientific
adoption, population accuracy, or calibrated confidence.

When the mid-occurrence seed round omits seed lists to respect its existing
20,000-anchor budget, its higher chain scores must not displace all extension
opportunities from a nonempty complete first round. Genomic extension therefore
uses that first round as its priority round. It keeps that round's chain order,
locus deduplication, best-chain constructors and decisive fast-lane eligibility.
Other rounds fill remaining shortlist slots in their existing order. If the
first round is absent or empty, the existing highest-best-chain round is used.
Complete repeat search keeps the previous highest-best-chain policy.

The final shortlist remains at most eight chains, and all selected chains share
two banded-DP fallback calls. Cheap existing extension paths remain available.
The complete first round can itself fill the shortlist or consume the entire
DP budget. This policy preserves its extension opportunity; it does not ensure
that every new repeat locus is extended or that every correct path is found.
Nor does first-round completeness certify a correct genomic placement.

Newly extended paths compete using the existing sequence/splice score,
eligibility, footprint deduplication and arbitration. There is no score bonus
for the first round, and a stronger repeat placement can win. Scores, primary
coverage (0.66), model weights and downstream filters do not change. Limited
search still propagates through cache, pairing and rescue and forces MAPQ 0.

Synthetic tests isolate scheduling with ordinary sequence, sparse first-round
anchors and stronger-chain-score repeat competitors. They cover both read
strands and outer clip orientations, short N tails and longer non-N tails,
strict-score competition, missing/empty first rounds, exact-repeat ambiguity,
the eight-chain limit and shared two-DP budget. They do not establish performance
on a representative population, independently unseen families or real samples.

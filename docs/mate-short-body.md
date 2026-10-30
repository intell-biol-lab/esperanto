# Observable paired short-body qualification

Dev.22 retains ordinary RNA coverage 0.66. Only an otherwise unplaced paired
read can use this separate qualification. `map --no-mate-short-body` disables
it independently; DNA, single reads and coordinate-only mate APIs do not use it.
The rejected exact-path competition experiment stays off by default.

The pool is the existing scored placements plus one existing bounded mate-window
proposal, rescored in the strict four-letter domain. Qualify before exact-path
deduplication; two distinct passing contig/start/strand/CIGAR paths refuse.
The prior candidate generator can still have omitted or collapsed alternatives.
No complete genomic search or global uniqueness is claimed.

A candidate must be positive but below ordinary coverage, have a legal exact
M/N body and one terminal canonical pure homopolymer S >=30. The adjacent body
base must differ from the tail. Do not trim matching boundaries to make them
eligible. Matching biological terminal runs and coincidental technical matches
remain unknown. The body is >=31 bases, has >=3 base types, and no base >80%.
Native minimizers inside M blocks must supply >=2 disjoint query intervals,
each with one primary-index posting agreeing with contig, coordinate and strand.
This is sketch support relative to the supplied index.

The complete current mate must be nonrescued, MAPQ>=20, same contig and opposite
strand, with union-intron transcript span <=1000 inclusive. Geometry is evaluated
in the effective-query frame before restoring raw sequence/quality/CIGAR.
Pair qualification is never sequence-cached. An accepted rescue cannot
bootstrap a second rescued mate.

Every result has MAPQ0 and RE:Z:mate_short_body. Exact indexed support is not a
calibrated probability, and a missing true source with identical represented
observations is unidentifiable. These placements do not enter the usual MAPQ20
pileup. They do not change two-pass discovery and do not reenter collapsed
rescue solely due to raw clipping. Existing ordinary mapped placements remain
unchanged except legitimate pair flags/coordinates/span when a mate recovers.

Candidate diagnostics keep ordinary coverage eligibility separate. The
mate_short_body pool adds a short_body decision, distinct passing path count
and disjoint seed (query start, canonical key) pairs. Query positions are in
raw alignment orientation, including minus-strand prefix lifting for supported
effective fragments. This schema does not claim all alternatives were searched.

Native tests cover both chains, mate orders, one/two pass, nonuniform qualities,
cache isolation, strict confidence/span bounds, four tail bases, prefix/suffix,
unobservable boundaries, spliced/indel/ambiguous paths and duplicated index
support. Development qualification and holdout/competitive advantage require
separate experiment results; software tests alone do not establish them.

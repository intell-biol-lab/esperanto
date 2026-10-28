# Splice rescue construction contract

Applies to the `1.1.1-dev.2` development candidate. Coordinates below are
zero-based and half-open. A junction skips `[start, end)` on the forward
reference; reverse-strand reads are already reverse-complemented before these
checks. Transcript strand does not reverse genomic CIGAR order.

For every accepted split or run-through rescue:

1. Positive-length `M`, `I` and `S` operations consume exactly the input read
   length. `S` may occur only at the two outer CIGAR ends.
2. Walking the CIGAR from the returned alignment position reaches `start`
   immediately before the introduced `N`; that `N` ends at `end`.
3. A locally extended left segment reaches the query split and donor with no
   unaligned inner suffix. A right segment begins at the query split and
   acceptor with no unaligned inner prefix or hidden reference offset.
4. When a primary body is retained, its actual CIGAR reference consumption,
   including indels, agrees with the declared junction boundary. A suffix
   retains its original reference position. Query offsets alone do not imply
   reference offsets.
5. Invalid candidates are rejected during selection. Channel A ranks valid
   retries; channel B tries scan candidates in descending score with stable
   ties; channel C and run-through try subsequent candidates after rejection.
   When no rescue qualifies the original placement remains available.
6. Right-tail reextension returns a new genomic start. Both the short-tail
   mismatch gate and the accepted placement use that start.

Outer soft clips and valid insertions/deletions are allowed. No internal soft
clip is silently converted into an insertion, no BAM records are discarded to
pass qualification, and no independent scoring semantics are changed.
This eligibility repair does not add boundary-constrained alignment search.

Executable synthetic checks live in `split.rs`, `span.rs` and `align.rs`:

```sh
cargo test --locked -p esperanto-map --lib
sh scripts/check.sh
```

The tests cover clipped seams, reference gaps at either boundary, exact
splices, outer clips, indels on both sides, strand-symmetric seam checks,
right-rescue position propagation and fallback to the original primary.
All fixtures are generated from synthetic sequences. Human qualification
requires new output directories, exact source/binary/input/reference identities
and strict independent scoring of both development seeds; software CI alone
does not establish scientific qualification.

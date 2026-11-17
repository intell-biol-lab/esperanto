# Scoped directional rescue indexes

`Index<'a>` and `IndexHit<'a>` now borrow tables/reference from their backing
storage. `Index::build`, `index_io::load` and `decode` retain the legacy
process-lifetime storage interface (`Index<'static>`). `OwnedIndex::open` shares
the same strict paidx decoder but owns its mmap and decoded packed reference;
`index()` lends a zero-copy view tied to that owner. No static-lifetime cast or
manual unmapping of borrowed memory is used. A compile-fail example protects
against dropping the owner while continuing to access its view.

The rescue flow preserves the four pass order AG, TC, TC, AG. Since dev.33 it
reuses one immutable TC owner across the two adjacent TC passes; each pass
still builds its own aligners, transforms its reads and verifies its placements.
After Rayon workers join and produce owned placements, the last pass using an
index releases its owner before verification. Consequently AG storage is gone
before TC loads, and TC is gone before the final AG loads. Ordinary Rust
unwinding/early returns release storage as well.
Standalone rescue owns its fallback four-letter reference for its whole call.
Read order, candidate retention, transforms, search budgets, gates, scores and
MAPQ are unchanged. Primary-map legacy storage remains process-lifetime; this
change does not claim to remove every historical leak or bound read buffering.

The mmap file must remain immutable during a run. Owned and legacy loaders both
check the same format/table invariants; checked table-size and offset arithmetic
also rejects integer overflow. This is a resource-lifetime change, with no paidx
format or index construction change. Parallel workers borrow shared immutable
views as before. Mmap unmapping releases process mappings/RSS; it does not evict
those file pages from the operating system's shared file cache.

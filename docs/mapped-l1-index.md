# Demand-paged L1 index access

Development version dev.30 keeps the L1 transcript index format and query
semantics unchanged. `L1Index` now opens `MappedTidx` instead of reading the whole
file and decoding both entry tables into separate heap arrays.

The mapped reader reuses the owning reader's layout validation: header fields,
supported versions, checked section sizes, trailing-byte rejection, UTF-8 names
and monotonic name/bucket offsets. The reader retains exactly the same existing
entry validation scope; this change does not add or remove key/payload checks.
Keys and payload fields are decoded explicitly from little-endian bytes, without
requiring aligned table starts or casting bytes to Rust tuples.

`lookup_hits` returns a borrowed view with the exact original hit count and
order. The engine's existing hit cap, budget counters and candidate selection
consume this view unchanged. Existing `Tidx` trait implementations automatically
adapt their slice lookup through the new method's default implementation.

The old owning `esperanto_tidx::Tidx` reader and slice API remain available.
`MappedTidx::lookup`, including the legacy `L1Index` trait slice interface,
materializes the payload table once when explicitly requested. Normal engine
mapping uses `lookup_hits` and never takes that compatibility path. Clones share
mapping and compatibility-table ownership through `Arc`; they do not duplicate
the large tables.

## Mapping lifetime and safety

The mapping belongs to the index object, and hit views borrow it. Dropping the
last owner releases the mapping. There are no new leaked mappings or static
borrow casts. The single unsafe mapping construction operation is isolated in
`tidx::mapped`; all layout parsing, lookup and iteration use safe Rust.

Index files must remain immutable while open, matching the genome index's mmap
precondition. Replacing a path with a new inode does not change an existing
mapping; mutating or truncating an open mapped inode is unsupported.

RSS can still include file-backed resident pages. Demand paging avoids eagerly
touching and copying whole tables; it does not guarantee that every workload
has a small resident working set, and cold storage can introduce page faults.
Report wall/CPU time, peak RSS and workload/cache conditions together.

## Required qualification

Reader tests cover v1/v2, unaligned sections, duplicate hits, empty buckets,
misses, full-width k=32, legacy views, concurrent queries, cloned owner lifetime,
and format-error parity. Human development qualification must additionally
compare complete alignment records and tags, candidate evidence, all scientific
counts and independent scores with matched dev.29 output.

No model, threshold, search budget, candidate count, scoring, MAPQ, reference or
ISA change is part of this implementation.

# Verified annotation junction cache

Mapping streams GTF records through a 1 MiB input buffer. Transcript grouping,
coordinate parsing, lossy UTF-8, line numbering, support counts and unknown contig
handling follow the previous whole-file parser. Input line/transcript sizes still
bound parser memory; this is not a fixed-memory transcript grouping algorithm.

By default map uses `$XDG_CACHE_HOME/esperanto/annotation-v1`, or
`$HOME/.cache/esperanto/annotation-v1`. `ESP_ANNOTATION_CACHE=/path` overrides it;
`ESP_ANNOTATION_CACHE=off` (or empty) bypasses cache I/O. If neither user directory
is available caching is disabled. Cache errors fall back to parsing; the stderr
outcome/timing does not enter scientific QC or alter mapping parameters.

Each invocation hashes the entire GTF. Keys include that digest, the reference
FASTA SHA256, ordered contig names/lengths, and parser/schema version. Each hit
also verifies the stored key/schema, canonical payload SHA256, sorted unique
junctions, positive counts and valid intervals, then reconstructs the acceptor
index. Initial parsing hashes consumed bytes again before publication. Metadata
alone is never used. Input/reference files must remain immutable during a run,
as required by the mmap index readers. Digests detect accidental corruption;
they do not authenticate a cache modified by an adversarial writer.

Only annotation junctions and support are cached. Pass-1 sample discoveries and
merged/sample-specific libraries are never cached. Entries are disposable and
limited to 128 MiB; larger entries use the normal parser. Writers create unique
exclusive temporary files, sync, then atomically rename; simultaneous writers
may publish the same deterministic content. Failed writes remove their temporary
file. No reference/index directory is modified. Cache pruning is user-managed;
old schema/source entries can be deleted at any time between runs.

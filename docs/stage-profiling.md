# Optional stage wall-time profiling

Set `ESP_PROFILE_STAGES=1` before starting the mapper to write stderr records:

```text
[esp-profile] flow.load.primary 0.850000
[esp-profile] primary.pe.read 0.250000
```

The process caches the setting at the first timer start. Other values disable
profiling; disabled guards do not read the clock or emit records. Guards occur
at stage/chunk boundaries, never once per read. Logging failures are ignored.
The records are engineering diagnostics, outside scientific statistics and SAM.

Top-level map components are `flow.load.primary`, `flow.annotation`,
`flow.load.l1`, `flow.primary.total`, `flow.primary.metadata` and `rescue.total`.
Their containing `flow.map.total` interval is inclusive. Process teardown can
continue after it ends; use full command wall time for user-visible cost.

Primary SE/PE phases are `primary.{se,pe}.read`, `.prepare`, `.search`, `.output`,
and `.finish`; PE adds `.pairs` for concrete pair arbitration. Each chunk emits
a separate interval, including a final EOF read. Two-pass `.2pass.pass1` and
`.2pass.pass2` are inclusive; library merging lies between them. Output timings
include record creation, compression-queue backpressure and writes, while
`.finish` waits for final sink completion. Compression can overlap computation,
so these intervals are not CPU-time attribution.

Rescue phases are `.prepare`, `.index`, `.search`, `.release`, `.verify`,
`.adjudicate` and `.output`. Four index intervals preserve pass order AG, TC,
TC, AG; the second TC interval measures borrowing the retained index. The owner
is retained only across that adjacency. Rescue output includes the existing
merge/rewrite and report updates; early empty-pool/empty-result exits omit
phases that did not execute.

Aggregate repeated names by summing their intervals, but never add an inclusive
parent to its children. Timings are wall seconds, rounded to six decimals.
Run profiling separately from ordinary performance comparisons, and verify
ordered alignment/evidence equality when enabling it.

# RUN_RECEIPT.json — the existing schema, as a contract

**Status:** descriptive, not aspirational. Everything below is read out of the
code that emits the file today. Nothing here is invented.

**Why this document exists.** Phase 1 adds a second writer, in `top.cpp`, so that
`--redock` stops producing runs with no engine identity. The moment there are two
writers there is a risk of two dialects, and every downstream parser then has to
sniff which one it is holding. The way to avoid that is to write the existing
shape down first and make the new writer conform to it — not to design a nicer
schema and leave 474 files speaking the old one.

**Sources.**

| What | Where |
|---|---|
| Field list, order, JSON encoding | `LIB/RunReceipt.cpp:53-84` (`build_run_receipt_json`) |
| Struct + types + defaults | `LIB/RunReceipt.h:19-47` |
| Field population | `LIB/DatasetRunner.cpp:5341-5396` |
| Nested `protocol_config` | `LIB/ProtocolConfig.cpp:345-438` (`ProtocolConfig::to_json`) |
| Hash helpers | `LIB/DatasetRunnerProvenance.cpp:44-61` |

---

## 1. Top-level fields

22 keys, emitted in exactly this order. `schema_version` is `1`
(`kRunReceiptSchemaVersion`, `RunReceipt.h:47`).

| # | Key | JSON type | Populated from |
|---|---|---|---|
| 1 | `schema_version` | int | compile-time constant `1` |
| 2 | `run_id` | string | `report.dataset_name`; if empty, `basename(config.output_dir)`; if still empty, `"dataset_run"` |
| 3 | `started_utc` | string | `utc_now_iso8601()` — `%Y-%m-%dT%H:%M:%SZ`, second resolution, no fractional part |
| 4 | `output` | string | `config.output_dir` (absolute) |
| 5 | `dataset` | string | `report.dataset_name` |
| 6 | `mode` | string | `receipt_mode_label`, e.g. `defined-cleft-redock` |
| 7 | `temperature_K` | number | `config.temperature`, cast to `double` |
| 8 | `pop` | int | `config.ga_population` |
| 9 | `gen` | int | `config.ga_generations` |
| 10 | `restarts` | int | `max(1, protocol_cfg_.restarts)` — never 0 |
| 11 | `seed_base` | int (uint64) | `protocol_cfg_.seed_base` |
| 12 | `seed_elitism` | **int 0/1** | derived; see §3 |
| 13 | `matrix_path` | string | resolved scoring-matrix path |
| 14 | `matrix_md5` | string | `md5 -q`, fallback `md5sum`; `""` on failure |
| 15 | `matrix_sha256` | string | `shasum -a 256`, fallback `sha256sum`; `""` on failure |
| 16 | `binary_path` | string | `flexaidds_bin` |
| 17 | `binary_sha256` | string | as above, over the binary |
| 18 | `runner_path` | string | `_NSGetExecutablePath` (Apple) / `/proc/self/exe` (Linux); `""` elsewhere |
| 19 | `runner_sha256` | string | as above, over `runner_path` |
| 20 | `git_commit` | string | `git rev-parse HEAD` via argv exec; `""` when it fails |
| 21 | `oracle_site_dir` | string | `protocol_cfg_.oracle_site_dir` |
| 22 | `oracle_site_dir_set` | **bool** | `!oracle_site_dir.empty()` |
| 23 | `protocol_config` | object | raw `ProtocolConfig::to_json()` text, 60 keys |

Verified against a production receipt
(`astex85_defcleft_claim_20260807_172821/run/RUN_RECEIPT.json`): all 22 keys
present, `protocol_config` an object of 60 keys.

---

## 2. Encoding rules a second writer must copy exactly

These are the parts that are easy to "improve" and thereby break.

**2.1 Two booleans, two encodings.** `seed_elitism` is emitted as `1`/`0`
(`RunReceipt.cpp:69`, ternary to int). `oracle_site_dir_set` is emitted as
`true`/`false` (`:79`). This is inconsistent, and it is the contract. A strict
parser written against `seed_elitism: 0` breaks the day a second writer emits
`false`. **Do not normalise this in Phase 1.** If it is to be fixed, it is a
`schema_version` bump with both writers changed together.

**2.2 Fixed 6-decimal floats everywhere.** The stream is configured
`o.setf(std::ios::fixed); o.precision(6)` in *both* `build_run_receipt_json`
(`:55-56`) and `ProtocolConfig::to_json` (`:347-348`). `temperature_K` therefore
serialises as `300.000000`, not `300`. A writer using default `ostream`
formatting produces numerically equal but textually different output, which
defeats byte-comparison of receipts across arms.

**2.3 Mixed pretty/compact style.** The top level is pretty-printed with
two-space indent and newlines; `protocol_config` is spliced in as compact,
space-free object text on a single line. That asymmetry is a consequence of
embedding pre-rendered JSON and is part of the on-disk shape.

**2.4 Escaping is partial.** `json_escape` (`RunReceipt.cpp:20-34`) handles only
`\` `"` `\n` `\r` `\t`. Other C0 control bytes pass through raw, which would
produce invalid JSON. In practice every string field is a path or an identifier,
so this has not bitten — but a second writer must not assume the helper is a
general-purpose escaper.

**2.5 Trailing newline.** `write_run_receipt` emits `body << "\n"`
(`:97`). The JSON text itself carries no trailing newline; the writer adds one.

**2.6 One call, two files.** `write_run_receipt(..., also_write_provenance_json)`
also writes a legacy slim `provenance.json` (9 keys: `dataset`, `matrix_path`,
`matrix_md5`, `matrix_sha256`, `binary_path`, `binary_sha256`, `git_commit`,
`oracle_site_dir`, `oracle_site_dir_set`, `protocol_config`). `DatasetRunner`
passes `true` (`:5396`). A Phase 1 writer must decide deliberately whether it is
also on the hook for that file; silently not writing it changes what older tools
find.

---

## 3. `seed_elitism` is derived, not copied

`DatasetRunner.cpp:5341-5347`:

```
receipt_seed_elitism = protocol_cfg_.seed_elitism
  ORACLE_CEILING                              -> forced true
  DEFINED_CLEFT_REDOCK | AUTONOMOUS | UNSET   -> forced false
```

The receipt records the **effective** value for the mode, not the configured
one. A second writer that copies `protocol_cfg_.seed_elitism` straight through
will emit a field with the same name and a different meaning — the worst
possible failure, because it is silent and only shows up as an inexplicable
cross-arm difference.

---

## 4. Known-weak fields

- **`git_commit` is frequently empty.** It shells out to `git rev-parse HEAD`
  with the runner's cwd. Drivers `cd` into the arm directory before launching, so
  unless that path is inside a checkout the call fails and the field is `""`.
  Observed empty in the production receipt sampled above. `run_t13_twotarget.sh`
  works around this by reading the commit the binary was *stamped* with at build
  time (`CMakeLists.txt:54-61`, `FLEXAIDS_GIT_COMMIT`) and refusing to launch an
  unstamped binary. That is the more reliable source and Phase 1 should prefer it.
- **`""` is the universal hash failure value.** `provenance_file_md5` /
  `provenance_file_sha256` return `""` when the path is empty, missing, or fails
  `is_safe_exec_path`. Empty is not "no hash requested"; it is "hashing failed".
  Consumers must distinguish those.
- **`runner_path`/`runner_sha256` are empty on any platform that is not Apple or
  Linux** — the `#if` has no `#else`.

---

## 5. The receipt is written *before* the docking runs

`write_run_receipt` is called at `DatasetRunner.cpp:5395`. The "Docking N
entries" loop starts at `:5405`. So `RUN_RECEIPT.json` is a **statement of
intent**, complete before a single pose exists.

Consequences, and they are the whole argument for the rest of Phase 0:

- The receipt contains no return code, no counts, no completion time, no
  outcome of any kind. It cannot.
- A receipt therefore proves a run was *configured*, never that it *finished*.
- `started_utc` is honestly named. There is no `finished_utc` and adding one
  would require a second write at the end — which is a different feature from
  Phase 1, not a free extra.
- This is precisely why `DONE` earns its place in the tier tables, and why the
  `KIND` sidecar carries a separate `status=` field. See `CONVENTION.md`.

---

## 6. The dialect collision already exists

Census over `~/flexaidds_results`, 2026-08-28:

```
RUN_RECEIPT.json files          476
containing "schema_version"     474   <- engine dialect, this document
not containing it                 2   <- hand-written by a shell driver
```

The two exceptions are under `ga1jd0_20260828_012729/`. They share a filename
with the engine dialect and nothing else — no `schema_version`, no
`protocol_config`, no `started_utc`. Instead: `schema: "1jd0_ga_wal400_v1"`,
`frozen_utc`, `engine_id`, `conditions`, `env_pins`, `assertions`, `supersedes`.

Both dialects are useful. The problem is only that they are not distinguishable
without opening the file and guessing. Two options:

1. **Reserve the filename.** `RUN_RECEIPT.json` means the engine dialect;
   campaign pre-registration goes in a differently named file (these two are
   effectively `PREREGISTRATION.json` — which is what
   `run_t13_twotarget.sh:151` already calls the same idea in text form).
2. **Require a discriminator.** Every writer emits a first key that identifies
   the dialect. Cheap for new files, but does not fix the 476 already on disk.

### Recommendation — not yet actioned

**Option 1**, on the reasoning that a rename fixes a census *retroactively*
while a discriminator key only fixes files written from now on. The 476 receipts
on disk are the corpus every existing analysis script walks; a fix that leaves
them ambiguous has not fixed the problem those scripts have. Two renames buy
correctness for all 476.

**Deliberately not done in this commit.** Both files live under
`ga1jd0_20260828_012729/`, which is Grok Bot's territory and adjacent to a live
batch. Renaming another seat's provenance artifacts is a coordination call
between LP and Grok Bot, not a change to slip into a documentation commit.
Recorded here so the decision is made on purpose rather than by whoever touches
the directory next.

Note also that the collision **predates this work**. It is the state of the
corpus as found on 2026-08-28, not something Phase 0 or Phase 1 introduces.
Phase 1 would be the third writer into an existing two-dialect situation, which
is an argument for settling the naming before Phase 1 lands, not after.

**Intersection, if both must coexist.** The keys present in both dialects are
`binary_sha256`, `runner_sha256`, `matrix_md5`, `pop`, `gen`, `restarts`, `mode`,
`seed_base`. A tool that reads only these works against either — which is a
reasonable floor for a census script, and a bad ceiling for anything else.

---

## 7. Contract for the Phase 1 `top.cpp` writer

1. Call `flexaids::write_run_receipt`. Do not hand-roll JSON.
2. Populate `RunReceiptInput` completely; leave no field defaulted-by-accident.
   `RunReceipt.h` defaults (`pop 1000`, `gen 2000`, `restarts 5`,
   `seed_elitism true`, `temperature_K 300.0`) are plausible-looking values that
   will silently misreport a `--redock` run that used anything else.
3. Reproduce the `seed_elitism` derivation of §3, or deliberately record why the
   single-target path differs.
4. Do not change `schema_version` unless keys change incompatibly — in which
   case both writers change in the same commit.
5. Prefer the build-stamped `FLEXAIDS_GIT_COMMIT` over shelling out to `git`.
6. Decide explicitly on `provenance.json` (§2.6).
7. `--redock` currently writes no receipt at all. Any receipt is an improvement;
   a receipt in a second dialect is not.

# stills_process → ImageSequence: prototype overview

Branch `stills_process_imagesequence` (dials + dxtbx). A prototype, **not intended for
merge** — its purpose is to surface the concrete changes, new objects, and sharp edges
a production effort will hit.

## What the branch does

`dials.stills_process` currently converts every imported image into a plain `ImageSet`
in `do_import()`. This branch keeps the data as an `ImageSequence` for the entire
pipeline. Two consequences follow:

- **Smaller output.** Experiments from one source file can share a single
  imageset/beam/detector object, so the `.expt` file stores one copy instead of one per
  frame.
- **Faster runtime.** Per-frame caching of the spotfinder and indexer becomes possible,
  because the cached objects now see a stable imageset identity across frames.

The rest of this document is the concrete list of what changed and why.

## Comparison to prior approach

A prior branch, `development_dsp_sequence_test`, explored the same core insight: remove the
`ImageSequence → ImageSet` conversion block in `Indexer.from_parameters()` that nulled scan
and goniometer for stills. That block is deleted here too. The approaches diverge in four ways:

**Direction of conversion.** The prior branch accepted stills arriving as `ImageSet` and
converted them *forward* to `ImageSequence` in `stills_process.py` via a bespoke
`convert_stills_to_sequences()` function. This branch keeps data as `ImageSequence` from the
format class layer — the conversion never happens.

**Stills detection.** The prior branch forced `stills.indexer = stills` as a phil default and
added `not experiment.is_still()` guards in the reflection predictor only. This branch
establishes a single contract everywhere: `expt.scan is None or expt.scan.is_still()`.

**Z-coordinate / bbox.** The prior branch had commented-out z-coordinate fixes in
`find_spots()` and `integrate()` — the bbox-incorrect problem was unresolved. This branch
fixes it upstream: per-frame slices are rebuilt with `Scan((1, 1))` in `do_import()`, so
no absolute frame range reaches the indexer.

**Scope.** The prior branch touched 5 files with no XFEL support or caching changes. This
branch adds the full XFEL architecture (`XFELBeam`, scan wavelength property,
`FormatXFEL` mixin) and the spotfinder and indexer caches.

## 1. The pivot: `do_import()` no longer converts to `ImageSet`

**Before:** `do_import()` called `ExperimentListFactory.from_filenames()` and then
`ImageSetFactory.imageset_from_anyset(experiment.imageset)`. That conversion rebuilt the
image data and handed back fresh, per-frame beam/detector model **copies**.

**After:** for multi-image formats it calls `format_class.get_imageset(as_sequence=True)`
directly, yielding an `ImageSequence` whose beam/detector are shared-identity Python
objects. The per-frame loop slices `imageset[i:i+1]` and keeps each slice an
`ImageSequence`.

Per-frame slices are `imageset[i:i+1]`, then rebuilt with a fresh `Scan((1,1))` to
normalize `array_range` to `(0,1)`. Every `ImageSequence` slice is rebuilt regardless
of whether the original scan was a still or a rotation — rotation CBF frames processed
via `convert_sequences_to_stills` carry absolute scan ranges (e.g., frame 2 has
`image_range=(2,2)` → `array_range=(1,2)`) that would confuse the integrator just as
much as XFEL absolute ranges.

**Why this was the blocker.** The old conversion destroyed *object identity*. The JSON
serializer deduplicates models by identity — same Python object written once. The
spotfinder mask cache (section 5) also keys on the imageset object. With per-frame
copies, neither could work. Removing the conversion is the root change; everything else
follows from it.

## 2. New / changed dxtbx objects

### `XFELBeam` — lightweight guard/marker type

**`XFELBeam`** (C++, extends `Beam`, in `model/beam.h` / `boost_python/beam.cc`) is a
**stateless** marker type: it carries direction, divergence, and polarization but
**no wavelength**. `get_wavelength()`, `set_wavelength()`, `get_s0()`, `set_s0()`, and
all scan-point methods deliberately **throw** — there is no single wavelength. Serializes
with `__id__: "xfel"`; `BeamFactory.from_dict` routes that id back to it.

The physical distinction matters: `PolychromaticBeam` (Laue/ToF) models a source with an
intrinsic spectral bandwidth and stores `wavelength_range`. `XFELBeam` models sequential
monochromatic pulses with varying energies — the per-shot energies are *measurements*,
stored in the scan properties table rather than in the beam model. Conflating the two is
incorrect.

A Python helper `get_monochromatic_beam(wavelength)` is injected onto `XFELBeam`
(in `model/__init__.py`). Given a wavelength in Ångström it returns a full monochromatic
`Beam` with the `XFELBeam`'s direction/divergence/polarization — this is called
internally by the `ImageSequence.get_beam(i)` override and is not a public consumer API.

### `ImageSequence.get_beam(index)` override — public per-frame beam API

Instead of a separate `XFELImageSequence` subclass, per-frame beam access is handled by
a Python override injected onto the existing `ImageSequence` class (in `imageset.py`).
When the imageset's beam is an `XFELBeam` **and** the scan has a `"wavelength"` property,
`get_beam(i)` reads `scan.get_property("wavelength")[i]` and dispatches to
`xfel_beam.get_monochromatic_beam(wl)`. All other cases — including `get_beam()` with no
index, non-XFEL beams, and non-XFEL imagesets — fall through unchanged.

This means callers in `stills_process.do_import` and `dials.import` that already call
`imageset.get_beam(i)` work without modification. The logic is encapsulated in dxtbx.

### Per-frame wavelengths in scan properties

The dxtbx scan properties table (`scan.set_property("wavelength", flex.double([...]))`)
stores per-image arbitrary properties. `FormatXFEL.get_imageset()` writes one entry per
frame when building an XFEL imageset. The property is read in three places:

1. `ImageSequence.get_beam(i)` — dispatches to `get_monochromatic_beam`
2. `stills_process.do_import` — stamps a single-entry `"wavelength"` onto each per-frame
   scan so the spotfinder can detect XFEL data on the in-flight path
3. `_rebuild_shared_imageset_output` — stamps per-frame wavelength onto each output scan
   so the property survives JSON round-trip

### `FormatXFEL` mixin + `FormatNXmxXFEL` / `FormatXTCXFEL`

The mixin's `get_imageset()` lets the parent format build the raw `ImageSet`, reads
per-frame wavelengths via a subclass `get_wavelengths()` (returning `list[float]`,
Ångström), builds an `XFELBeam` from the first frame's direction/divergence, and returns
a regular `ImageSequence` with that `XFELBeam` and a zero-oscillation `Scan` carrying the
`"wavelength"` property.

For NXmx, `get_wavelengths()` calls `CachedWavelengthBeamFactory.get_wavelengths()` (in
`dxtbx/nexus/__init__.py`), which reads the 1-D `incident_wavelength` HDF5 dataset.

`FormatNXmxXFEL.understand()` distinguishes genuine per-pulse XFEL data from
constant-wavelength files: it requires `wl.ndim > 0 and wl.size > 1`. This excludes both
scalar wavelengths (`ndim=0`) and 1-element arrays (`size=1`), which both represent a
fixed wavelength for all frames — those files use regular `FormatNXmx` instead.

`FormatXTCXFEL` delegates `understand()` to `FormatXTC` since XTC streams are always
XFEL stills with per-event wavelengths.

## 3. Cross-cutting change: how a "still" is detected

Once stills stay as `ImageSequence`, the old test `isinstance(imageset, ImageSequence)`
no longer separates stills from rotations — **both are now `ImageSequence`**. Every such
check had to move to the scan:

> a still is `scan is None or scan.is_still()` (zero-oscillation).

Touched: `Indexer.from_parameters` (stills/sequence routing),
`indexer._xyzcal_mm_to_px`, `reflection_predictor` (static vs scan-varying predictor
dispatch), `integrator._determine_max_memory_needed`, and `SpotFinderFactory` (via a new
`is_stills` argument). It also let an entire block be **deleted** from
`Indexer.from_parameters` — the code that used to convert `ImageSequence`→`ImageSet` and
null scan/goniometer is now unnecessary.

XFEL detection in the spotfinder (`finder.py`) is a separate concern: it uses
`scan.has_property("wavelength")` rather than any beam isinstance check — per-frame
imagesets carry monochromatic beams built by `do_import`, so an `XFELBeam` check would
always be False on the in-flight path.

This is the widest blast radius in the branch and the change most in need of a
pipeline-wide decision: is scan-based detection the contract we want everywhere?

## 4. Output: shared models per source file

`_rebuild_shared_imageset_output()` (new, in `stills_process.py`) runs before writing.
It groups experiments by source file and rebuilds them so every experiment from one file
points at **one** `ImageSequence` / beam / detector object, and multi-lattice
experiments from the same frame share **one** `Scan` object. This is the change that
actually shrinks the file. It is a no-op for rotation data. A helper `_frame_index()`
(new) recovers a frame's 0-based index, handling both representations — the per-frame
imageset from `do_import` and the full shared imageset after a JSON round-trip.

`combine_all_ranks` (new phil option) merges the per-rank (MPI) or per-worker
(multiprocessing) composite files into a single combined file with shared models. The
multiprocessing path is handled by `_combine_multiprocessing_outputs()` (new). As part
of this, `extend_with_bookkeeping` — which remaps reflection `id` and
`experiment_identifiers` when merging — was lifted out of the MPI `finalize` branch into
a module-level function so both the MPI and multiprocessing paths reuse it.

### 4a. Compact scan serialization

Even after `_rebuild_shared_imageset_output()` runs, the output still contains one
`Scan` JSON object per integrated image. Although individual scan properties are small,
the per-object overhead accumulates: a 764-experiment file carries 764 scan dicts totaling
285 KB (8.5 % of 3.3 MB).

`ExperimentList.to_dict/as_json` gains a `compact_stills_scans=True` keyword (in
`src/dxtbx/model/__init__.py`). When set, and all scans are single-frame stills
(oscillation width == 0), the N scan dicts are merged into **one** consolidated object:

```json
{"__stills_consolidated": true, "batch_offset": 0,
 "frame_numbers": [575, 821, ...],
 "properties": {"wavelength": [1.3087, 1.3078, ...]},
 "valid_image_ranges": {}}
```

All-zero property arrays (epochs, exposure_time, oscillation, oscillation_width) are
omitted and reconstructed on load. Each experiment gains `"scan_point": i` to identify
its position in the consolidated arrays; all point to `"scan": 0`.

On the read side, `ExperimentListDict._expand_consolidated_scans()` (in
`src/dxtbx/model/experiment_list.py`) runs in `__init__` before any model extraction.
It expands the consolidated object back to N per-frame scan dicts, so the rest of
`decode()` sees the standard format unchanged. The method is a no-op for old files.

All four composite output `as_json()` calls in `stills_process.py` pass
`compact_stills_scans=True`. Measured reduction on a 764-experiment file: **−267 KB
(8.1 %)**.

**Backward compatibility:** new reader + old file is transparent. Old reader + new file
fails (no `image_range` in the scan dict), but old readers are never exposed to the
new format unless the writer opts in via the flag.

## 5. Performance caching (unblocked by section 1)

**Spotfinder** (`finder.py` + `stills_process.py`). `ExtractPixelsFromImage` splits its
mask into a *static* component (lookup mask + untrusted rectangles, frame-independent)
and a *dynamic* per-frame component (`imageset.get_mask(index)`). A new
`update_imageset()` swaps in the next frame's data without rebuilding. `ExtractSpots`
caches the `ExtractPixelsFromImage` object and reuses it across stills frames. The
`Processor` itself also caches the configured `SpotFinderFactory` in
`self.spot_finder_factory` (built once by the new `get_spot_finder_factory()`), so the
factory is not reconstructed per frame.

**Indexer** (`stills_process.py`). `Processor` caches the `Indexer` object
(`self.idxr`, `self.idxr_method_list`, `self.idxr_known_crystal_models`); a new
`update_indexer()` resets per-frame state instead of reconstructing. `copy.deepcopy` of
the params object per image was replaced with direct mutate-and-restore of two fields:
`refinement.parameterisation.scan_varying` (forced `False`) and
`basis_vector_combinations.max_refine`. Note the latter is now pinned to `5` for stills
indexing — a behavioral change, not just a refactor — and restored afterward.

Both caches depend on section 1 — they only work because the imageset identity is now
stable across frames.

## 6. Surprises / sharp edges

- The **C++ `ImageSequence` requires a non-null scan**, but stills have none. Two
  separate code paths work around this, and they should not be conflated:
  (a) for XFEL data, `FormatXFEL.get_imageset()` builds a **real** zero-oscillation
  `Scan((1, n), (0.0, 0.0))` directly — there is no error to catch;
  (b) for non-XFEL multi-image formats, `format_class.get_imageset(as_sequence=True)`
  raises a `RuntimeError` (its `get_scan()` returns `None`), so `do_import` catches that,
  loads the data as an `ImageSet`, and wraps it in an `ImageSequence` with a
  **synthetic** zero-oscillation scan.
- **`load_models()` re-reads the synthetic still scan** back onto the experiment, so the
  scan/goniometer must be re-nulled after *every* `load_models()` call (both `do_import`
  and `do_work`).
- A **per-frame slice carries the absolute frame range** — XFEL frame 3913 gives
  `array_range=(3912,3913)`, and a rotation CBF frame 2 gives `array_range=(1,2)` even
  after `convert_sequences_to_stills=True`. Both conflict with the integrator's
  expectation of 0-based z coordinates. Every `ImageSequence` slice is rebuilt with
  `Scan((1,1))` in `do_import`, not just XFEL stills.
- **`Format.get_imageset` asserts `goniometer is not None`** in the sequence path —
  which fires when reading back output `.expt` with `check_format=False` because
  stills have `goniometer=None`. Fixed in dxtbx `Format.py` (~line 441): the assertion
  now only fires for non-still rotation sequences
  (`if goniometer is None and not (scan is not None and scan.is_still())`).
- **`FormatXTC.understand()` returned `True` unconditionally** when psana was importable,
  because the return value of `_get_datasource()` was discarded. This caused
  `FormatXTCXFEL` (the concrete registry class) to match arbitrary temp files in tests.
  Fixed: `return bool(ds)`.
- **`imageset.get_beam()` with no index returns `XFELBeam`** for XFEL imagesets — which
  throws on `get_s0()`. Callers that need a monochromatic beam must call `get_beam(i)`.
- **XFEL detection via scan property, not isinstance.** The spotfinder checks
  `scan.has_property("wavelength")` to suppress caching for XFEL data. Checking
  `isinstance(imageset.get_beam(), XFELBeam)` would silently fail on the in-flight path
  because per-frame imagesets have monochromatic beams, not `XFELBeam`.
- **`FormatNXmxXFEL` vs `FormatNXmx` selection.** `FormatNXmxXFEL.understand()` requires
  `wl.size > 1`. Files with a 1-element `incident_wavelength` array (constant-energy LCLS
  runs) use `FormatNXmx` instead — they have no per-pulse wavelength variation to capture.
- `easy_mp` tried to **pickle C++ threshold/indexer objects** in the result queue;
  `do_work` now returns `None` when `finalize=True` since output is already on disk.

## 7. Limitations & open questions (feedback wanted)

- Should `XFELBeam` become a permanent, supported dxtbx model?
- Is scan-based stills detection the right pipeline-wide contract?
- `test_pseudo_scan` now passes and asserts `ImageSequence` output with `goniometer=None`
  and still scan after round-trip. Dedicated unit tests for `XFELBeam`, scan
  `"wavelength"` property, and `ImageSequence.get_beam(i)` are not yet written.
- Downstream consumers (cctbx.xfel, merging) have not been run against the new output.
- `combine_all_ranks` output semantics are not yet settled.
- No performance numbers yet — file-size and runtime benchmarking is a next step.
- Crystal section (67 % of a typical file, ~2.2 MB) is the dominant remaining file-size
  cost; scan consolidation (section 4a) addressed the 8 % scan share. Crystals require a
  separate investigation.

## 8. How it would land

The work decomposes into reviewable pieces: **1 dxtbx PR** (the new objects +
`XFELBeam` + compact scan serialization) and **3 dependent dials PRs** — correctness
transition (including `compact_stills_scans=True` in composite output), then spotfinder
caching, then indexer caching. The dxtbx PR must merge first so dials can import the new
types. Note the two-step rollout for compact scans: the dxtbx read-side change
(`_expand_consolidated_scans`) must be deployed before the write-side flag is enabled in
stills_process.

## 9. Image viewer: reflection overlays for composite still output

`dials.image_viewer` loads the composite `.expt`/`.refl` files produced by section 4.
Without additional fixes, the reflection overlays (spots, centres of mass, predictions,
hkl labels, integrated shoeboxes) appeared on the first frame only; all later frames were
blank. Two independent root causes were found and fixed (commit `875d1a8cc`).

### 9.1 Root cause: `_identifiers_for_frame` binary search on unsorted list

`SpotFrame._identifiers_for_frame` (`spotfinder_frame.py`) selects, for the current
frame, the experiment identifiers whose per-experiment scan starts at that frame. It
then calls `select_on_experiment_identifiers()` on the composite reflection table to
extract that frame's rows. The original implementation binary-searched the experiment
list, with a docstring asserting "sorted ascending by frame index, guaranteed by the
shared ImageSequence ordering." This was false for two reasons:

- `_rebuild_shared_imageset_output` (section 4) preserves the input experiment order —
  it does not sort by absolute frame index.
- `_combine_multiprocessing_outputs` appends per-worker composite files end-to-end.
  Worker 0 processed frames {0, N, 2N, …}, worker 1 processed {1, N+1, 2N+1, …}, so
  the merged list interleaves frames across strides rather than presenting them in
  globally ascending order.

Binary search on an unsorted list returns incorrect or empty results for all frames
except the one that happens to be first. `identifiers=[]` → zero rows selected →
no overlay drawn.

**Fix:** replaced the binary search with an order-independent lookup:

```python
def _identifiers_for_frame(self, i_frame):
    frame_map = {}
    for elist in self.experiments:
        for expt in elist:
            if expt.scan is None:
                continue
            frame_map.setdefault(expt.scan.get_array_range()[0], []).append(expt.identifier)
    return frame_map.get(i_frame, [])
```

O(N_expt) per call, no ordering assumption.

### 9.2 Root cause: shared imageset scan loses start offset on JSON round-trip

`_rebuild_shared_imageset_output` builds the shared `ImageSequence` with a scan whose
`image_range` spans the full set of integrated frames — e.g., `Scan((6, 1449))` when
the earliest integrated frame is frame 5. After JSON round-trip the imageset's own scan
comes back as `Scan((0, 1443))`: the absolute start offset is discarded. (The
`imageset.indices()` list — e.g., `[5, 6, …, 1447]` — does survive round-trip via the
imageset's data block, but the scan's `array_range` does not.)

The frame index computation in `get_spotfinder_data` originally was:
```python
i_frame += imageset.get_scan().get_array_range()[0]   # returns 0, not 5
```
so chooser position 0 mapped to `i_frame=0`, while every per-experiment scan starts at
frame 5. `_identifiers_for_frame(0)` returned `[]` even after the dict fix.

**Fix:** new helper `_still_scan_frame_offset(imageset)`. Per-experiment scans are
individually serialized and keep their correct absolute frame index through round-trip.
The helper finds `min(expt.scan.get_array_range()[0] for expt in all_experiments_for_imageset)`
and caches the result per imageset (keyed by `id(imageset)`). `get_spotfinder_data` uses
this offset for the `viewing_still_scans` code path:
```python
if self.viewing_still_scans:
    i_frame += self._still_scan_frame_offset(imageset)
```

### 9.3 Per-overlay colours for composite output

The original colour selection indexed `prediction_colours[reflection["id"]]` — a
90-entry list. In the old `viewing_stills` path, `id` was reset to 0 per frame and
never exceeded bounds. In the new `viewing_still_scans` path, `id` is the composite
table's experiment id (0–763 for 764 experiments), which crashes and also produces
frame-varying colours.

Fix: for `viewing_still_scans`, integrated shoeboxes use a fixed purple (`#984ea3`) and
predictions a fixed green (`#4daf4a`), independent of experiment id. Non-still modes
keep the original `prediction_colours[id]` behaviour. The `all_pix` colour sites retain
a modulo guard as crash protection for rare multi-lattice frames.

### 9.4 Per-frame beam (`slip_viewer/frame.py`)

`chooser_wrapper.get_beam()` was changed from `image_set.get_beam()` to
`image_set.get_beam(self.index)` to return the per-frame monochromatic beam. The beam
centre and resolution ring drawing use the same pattern.

### 9.5 Performance of frame-flipping with large composite files

The per-frame overlay computation has two O(N) costs that scale poorly:

| Path | Complexity | Notes (764-expt / 484,603-refl test case) |
|------|-----------|-------------------------------------------|
| bbox / ctr-mass / shoebox overlays | O(N_refl) scan per flip | `_reflection_overlay_data` iterates all reflections |
| predictions | O(N_expt × N_refl) per flip | outer loop over all experiments, inner scan-and-select |

The predictions loop performs roughly 3.7 × 10^8 operations per frame flip in the
test case. The natural fix is to precompute a `{abs_frame_index: flex.size_t(row_indices)}`
dict at load time, reducing both paths to O(rows-on-frame) ≈ 630 rows per flip. This
has not been implemented.

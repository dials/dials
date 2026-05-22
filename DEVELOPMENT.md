# Development Context: `stills_process_imagesequence` branch

**Branch:** `stills_process_imagesequence` in both dials and dxtbx  
**Status as of 2026-05-21:** Architecture refactored per upstream dxtbx maintainer review (Waterman, McDonagh). `XFELImageSequence` removed; per-frame wavelengths moved to scan properties table. `XFELBeam` retained as lightweight guard/marker type (no wavelength array). Build passes; all smoke tests pass.  
**2026-05-21 (scan consolidation):** `ExperimentList.to_dict/as_json` gains `compact_stills_scans=True` option; stills composite output now writes one consolidated scan JSON object instead of N, saving ~8 % of file size (−266 KB on a 764-experiment test case).  
**2026-05-21 (bug fixes):** Per-frame scan rebuild in `do_import` extended to all `ImageSequence` slices (not only XFEL stills) — rotation CBF frames processed via `convert_sequences_to_stills` also carry absolute scan ranges that would confuse the integrator. `FormatXTC.understand()` fixed to `return bool(ds)` instead of `True`. `Format.get_imageset` goniometer assertion relaxed for still sequences. `test_pseudo_scan` fixed and strengthened with round-trip assertions.  
**2026-05-22 (image_viewer):** Reflection overlays (spots, centres of mass, predictions, hkl, integrated shoeboxes) now draw correctly on all frames of composite still output in `dials.image_viewer`. Two independent root causes fixed; per-frame beam corrected for XFEL. Commit `875d1a8cc`. See "image_viewer changes" section below.  
**2026-05-22 (combine_experiments):** `dials.combine_experiments` now works with stills-as-ImageSequence output. Four fixes across both repos; see "combine_experiments changes" section below.
**2026-05-22 (frame ordering):** `stills_process` and `dials.combine_experiments` composite output now writes experiments (and paired reflections) in ascending frame order. New helper `_sort_experiments_by_frame()` sorts by `(source path, frame index)` and remaps reflection-table `id` and `experiment_identifiers` to match. Called in both the multiprocessing path (`_combine_multiprocessing_outputs`, after `_rebuild_shared_imageset_output`) and the MPI composite-stride path (`finalize`, after the rebuild loop, for all four stages). Previously, experiments were written in worker/rank encounter order.

---

## What this branch does (one paragraph)

`dials.stills_process` normally converts every imported image into a plain `ImageSet` (one object per frame). This branch keeps data as an `ImageSequence` for the entire pipeline. Two payoffs: (1) experiments from one source file share a single imageset/beam/detector, so `.expt` files shrink; (2) the spotfinder and indexer can be cached across frames because imageset identity is now stable. The pivot is `do_import()` no longer calling `imageset_from_anyset()` — that conversion was destroying object identity and blocking both the JSON deduplication and the caches. For XFEL data specifically, per-frame wavelengths are stored in `scan.get_property("wavelength")` (the dxtbx scan properties table), and the shared beam is an `XFELBeam` (a guard/marker type that throws on `get_wavelength()`/`get_s0()`). `ImageSequence.get_beam(i)` dispatches to `XFELBeam.get_monochromatic_beam(wl)` to return a per-frame monochromatic `Beam`.

---

## Repo layout

| Repo | Path |
|------|------|
| dials | `/pscratch/sd/d/dwmoreau/cctbx_05_14_26/cctbx/modules/dials` |
| dxtbx | `/pscratch/sd/d/dwmoreau/cctbx_05_14_26/cctbx/modules/dxtbx` |

Run Python as: `libtbx.python` (not bare `python3`)

---

## dxtbx changes

### New/changed types

**`XFELBeam`** — C++ class inheriting `Beam` (lightweight guard/marker, no wavelength array)
- Files: `src/dxtbx/model/beam.h` (~lines 880–1045), `src/dxtbx/model/boost_python/beam.cc` (~lines 423–715)
- Analogous to `PolychromaticBeam` but **no** `wavelength_range_` — XFEL pulses are individually monochromatic; per-shot energies are measurements in scan properties, not a beam property
- Guards: `get_wavelength()`, `set_wavelength()`, `get_s0()`, `set_s0()`, all scan-point methods **throw** — intentional
- Constructors: default, `(direction, divergence, sigma_divergence)`, full with polarization/flux/transmission/probe
- `operator==` and `is_similar_to` compare direction + divergence + polarization (no wavelength comparison)
- Serializes with `"__id__": "xfel"`; no `wavelengths` key in dict
- Factory: `BeamFactory.make_xfel_beam(direction, divergence, sigma_divergence)` in `src/dxtbx/model/beam.py`
- Python helper `get_monochromatic_beam(wavelength)` injected in `src/dxtbx/model/__init__.py` — returns a monochromatic `Beam` with given wavelength; called by `ImageSequence.get_beam(i)`

**`XFELImageSequence` — REMOVED**
- Replaced by `ImageSequence.get_beam(index)` Python override in `src/dxtbx/imageset.py`
- When imageset beam is `XFELBeam` and scan has `"wavelength"` property: `get_beam(i)` → `xfel_beam.get_monochromatic_beam(wl[i])`
- All other cases (index=None, non-XFELBeam): returns shared beam unchanged

### Format classes

- **`FormatXFEL` mixin** (`src/dxtbx/format/FormatXFEL.py`): `get_imageset()` builds regular `ImageSequence` with `XFELBeam` (no wavelength) + zero-oscillation `Scan` with `"wavelength"` property. Subclasses implement `get_wavelengths()`. **Wavelength subset fix:** when the raw imageset covers only a subset of the source file's frames (e.g. loading a composite stills `.expt` output), wavelengths are now indexed by `raw_iset.indices()` rather than assumed to span the full file.
- **`FormatNXmxXFEL`** (`src/dxtbx/format/FormatNXmx.py`, ~lines 180–200): unchanged. `get_wavelengths()` delegates to `CachedWavelengthBeamFactory.get_wavelengths()`.
- **`FormatXTCXFEL`** (`src/dxtbx/format/FormatXTC.py`, ~lines 952–967): unchanged. `get_wavelengths()` calls `self.get_beam(i).get_wavelength()` per frame. **`FormatXTC.understand()` now returns `bool(ds)`** instead of `True` — previously the return value of `_get_datasource()` was discarded, causing `FormatXTCXFEL` to match arbitrary files whenever psana was importable.
- **`CachedWavelengthBeamFactory.get_wavelengths()`** (`src/dxtbx/nexus/__init__.py`, ~lines 82–93): renamed from `make_xfel_beam()`. Returns `list[float]` (Ångström) instead of `XFELBeam`.

---

## dials changes (10 commits)

### `src/dials/command_line/stills_process.py` — largest file (~592 ins / ~35 del)

**New helpers:**
- `_make_stills_sequence()` (~line 359): Builds an `ImageSequence` with zero-oscillation `Scan` and no goniometer — shared by `do_import()`, its error fallback, and `_rebuild_shared_imageset_output()`.
- `_frame_index(expt)` (~line 508): Recovers 0-based frame index from either a per-frame imageset (`imageset.indices()`) or a full-sequence imageset (`scan.get_array_range()`). Needed because output rebuilding uses different representations than import.
- `extend_with_bookkeeping()` (~line 492): Module-level function (lifted out of MPI `finalize`) that remaps reflection `id` and `experiment_identifiers` when merging experiment lists. Used by both MPI and multiprocessing paths.
- `_sort_experiments_by_frame(expts, refls)` (~line 535): Sorts experiments by `(imageset.paths()[0], _frame_index(expt))` and remaps reflection-table `id` and `experiment_identifiers` dict to match. Returns inputs unchanged if already sorted or if `expts` is empty. Called after `_rebuild_shared_imageset_output` in both output paths.

**`do_import()` refactor (~lines 379–489):**  
For multi-image formats: calls `format_class.get_imageset(as_sequence=True)` — gets a shared-identity `ImageSequence`. If that raises `RuntimeError` (format's `get_scan()` returns None), wraps the `ImageSet` in a synthetic `ImageSequence` via `_make_stills_sequence()`. Per-frame slices are `imageset[i:i+1]`, then rebuilt with a fresh `Scan((1,1))` to normalize `array_range` to `(0,1)`. Every `ImageSequence` slice is rebuilt regardless of whether the original scan was a still or a rotation — rotation CBF frames processed via `convert_sequences_to_stills` carry absolute scan ranges (e.g., frame 2 has `image_range=(2,2)` → `array_range=(1,2)`) that would confuse the integrator just as much as XFEL absolute ranges.

**Key detail:** `load_models()` re-reads the synthetic still scan back onto the experiment, so scan/goniometer must be re-nulled after *every* `load_models()` call (in both `do_import` and `do_work`).

**Output rebuilding:**
- `_rebuild_shared_imageset_output()` (~line 530): Groups experiments by source file; re-opens source to recover full `ImageSetData`; creates one `Scan((fi+1, fi+1))` per unique indexed frame; multi-lattice frames from the same frame share one `Scan` object. Guards against rotation data (returns unchanged). This is what actually shrinks the file.
- `_combine_multiprocessing_outputs()` (~line 615): Merges per-worker composites, calls `_rebuild_shared_imageset_output()`, rewrites as worker-0 file, deletes intermediates.
- `combine_all_ranks` (new phil): Triggers merge at end. MPI path sets `composite_stride = comm.size` so rank 0 gets all frames.

**Caching in `Processor`:**  
- `self.spot_finder_factory` (~line 1259): Built once by `get_spot_finder_factory()`; reused per frame.
- `self.idxr`, `self.idxr_method_list`, `self.idxr_known_crystal_models` (~line 1263): Indexer cached; `update_indexer()` (~line 1569) resets per-frame state instead of reconstructing. `copy.deepcopy(all_params)` removed — replaced with direct mutate-and-restore of `refinement.parameterisation.scan_varying` (forced False) and `basis_vector_combinations.max_refine` (pinned to 5, then restored). **The `max_refine=5` pin is a behavioral change, not just a refactor.**

### `src/dials/algorithms/spot_finding/finder.py` (~142 ins / ~135 del)

- `ExtractPixelsFromImage`: Added `is_stills` param; renamed `self.mask` → `self.image_mask` (static component only); `__call__()` merges static mask with dynamic `imageset.get_mask(index)` per call.
- New `update_imageset()` method: Swaps in next frame's data without rebuilding. This is how caching works — static mask survives, dynamic mask is re-read.
- `ExtractSpots` and `SpotFinder._find_spots()`: Cache `ExtractPixelsFromImage` in `self.function`; call `update_imageset()` on subsequent frames instead of constructing new.

### `src/dials/algorithms/spot_finding/factory.py` (~17 ins / ~15 del)

`SpotFinderFactory.from_parameters()` now accepts explicit `is_stills` parameter. When True, sets `no_shoeboxes_2d = True` directly rather than checking `isinstance(imageset, ImageSequence)`.

### `src/dials/algorithms/indexing/indexer.py` (~45 ins / ~60 del)

- Stills detection changed from `isinstance(expt.imageset, ImageSequence)` → `expt.scan is None or expt.scan.is_still()`. Applies to `from_parameters()` routing and `_xyzcal_mm_to_px()` guard.
- Deleted: the block that used to convert `ImageSequence → ImageSet` and null scan/goniometer for stills — no longer needed.

### `src/dials/algorithms/indexing/stills_indexer.py` (~67 ins / ~52 del)

- `self.hardcoded_phil` moved to `__init__()` (was re-parsed per call to `choose_best_orientation_matrix()`).
- `copy.deepcopy(self.all_params)` removed from `choose_best_orientation_matrix()` — params are not mutated here; mutation+restore happens in `Processor.index()`.
- `experiment_list_for_crystal()`: Loop changed from iterating imagesets to iterating experiments, so beam comes from `expt.beam` (wavelength-adjusted for XFEL) not `imageset.get_beam()` (base shared beam).

### `src/dials/algorithms/integration/integrator.py` (~28 ins / ~23 del)

`_determine_max_memory_needed()`: Prefers experiment's scan over imageset's scan for stills detection. Condition now checks `expt_scan is None or expt_scan.is_still()`.

### `src/dials/algorithms/spot_prediction/reflection_predictor.py` (~10 ins / ~18 del)

Stills detection: `isinstance(experiment.imageset, ImageSequence)` → `experiment.scan is not None and not experiment.scan.is_still()`. Removed `from dxtbx.imageset import ImageSequence` import.

---

## Cross-cutting contract change

**Old:** stills detected by `isinstance(imageset, ImageSequence)` — False for `ImageSet`, True for `ImageSequence`.  
**New:** stills detected by `expt.scan is None or expt.scan.is_still()` — works regardless of imageset type.  

Every place this changes is a potential regression site if any code path was missed. Touched: indexer, integrator, spotfinder factory, reflection predictor, stills_process itself.

---

## Sharp edges / gotchas

1. **C++ `ImageSequence` requires a non-null scan.** Two separate workarounds — don't conflate them:
   - XFEL formats: `FormatXFEL.get_imageset()` builds a zero-oscillation `Scan((1, n), (0.0, 0.0))` with `"wavelength"` property — no error to catch.
   - Non-XFEL multi-image formats: `get_imageset(as_sequence=True)` raises `RuntimeError`, so `do_import` catches it, loads as `ImageSet`, wraps with synthetic scan via `_make_stills_sequence()`.

2. **Per-frame scan range is absolute** for any `ImageSequence` slice — XFEL frame 3913 has `array_range=(3912,3913)`, and rotation CBF frame 2 has `array_range=(1,2)`. The integrator expects 0-based z. Every `ImageSequence` slice is rebuilt with `Scan((1,1))` in `do_import` regardless of whether the original scan was a still or a rotation.

3. **`Format.get_imageset` goniometer assertion is relaxed for still sequences.** The original `assert goniometer is not None` in the sequence path now only fires for non-still (rotation) sequences. This is required so output `.expt` files (where `expt.goniometer=None` for stills) can be read back with `check_format=False`. File: `dxtbx/src/dxtbx/format/Format.py` (~line 441).

4. **`load_models()` re-reads the synthetic still scan** back onto the experiment. Scan/goniometer must be re-nulled after every `load_models()` call.

5. **`ImageSequence.get_beam()` with no index returns `XFELBeam`** for XFEL data. Callers that need a monochromatic beam must call `get_beam(i)` — the injected override dispatches to `XFELBeam.get_monochromatic_beam(wl[i])`.

6. **Spotfinder XFEL detection uses `scan.has_property("wavelength")`**, not `isinstance(beam, XFELBeam)` — because per-frame imagesets carry monochromatic beams built by `do_import`, not XFELBeam.

7. **`easy_mp` tried to pickle C++ threshold/indexer objects** via the result queue. Fixed: `do_work` returns `None` when `finalize=True` since output is already on disk.

8. **`max_refine=5` pin** in the indexer cache path is a behavioral change. It was previously derived from params; now it's hardcoded to 5 for stills and restored afterward.

---

## `.expt` file size — scan consolidation (2026-05-21)

### Problem

In composite output mode, each integrated image produces its own `Scan` JSON
object even though scans differ only in `image_range` (one integer) and
`wavelength` (one float). For a 764-experiment file this costs 285 KB (8.5 %
of a 3.3 MB file) for data that could be stored in two compact arrays.

### Solution: `compact_stills_scans=True`

**Write side — `ExperimentList.to_dict/as_json` in `src/dxtbx/model/__init__.py`**

When `compact_stills_scans=True` and every scan is a single-frame still
(`image_range[0] == image_range[1]`, oscillation width == 0):

- The N scan dicts are replaced by **one** consolidated object:
  ```json
  {"__stills_consolidated": true, "batch_offset": 0,
   "frame_numbers": [575, 821, ...],
   "properties": {"wavelength": [1.3087, 1.3078, ...]},
   "valid_image_ranges": {}}
  ```
- All-zero property arrays (epochs, exposure_time, oscillation,
  oscillation_width) are omitted; they are reconstructed as zeros on load.
- Each experiment gains `"scan_point": i` and `"scan": 0`.

**Read side — `ExperimentListDict._expand_consolidated_scans()` in `src/dxtbx/model/experiment_list.py`**

Called in `__init__` before `_extract_models` runs. Expands the consolidated
object back to N per-frame scan dicts (no-op for files without
`__stills_consolidated`). The rest of `decode()` sees the standard format
unchanged.

**Caller — `src/dials/command_line/stills_process.py`**

All four composite output `as_json()` calls pass `compact_stills_scans=True`.

### Numbers (764-experiment test case)

| | Size |
|---|---|
| Original | 3,269 KB |
| Compact | 3,002 KB |
| Savings | **−267 KB (8.1 %)** |

### Backward compatibility

- New reader + old file: **transparent** (`_expand_consolidated_scans` exits immediately).
- Old reader + new file: **fails** (`ScanFactory.from_dict` gets a dict with no `image_range`). Safe because `compact_stills_scans=False` is the default; old readers are never exposed to the new format unless the writer explicitly opts in.
- Rollout: deploy dxtbx change (reader) before flipping the writer flag in stills_process.

### Remaining file-size opportunities

| Section | Size | % | Notes |
|---|---|---|---|
| crystal | 2.24 MB | 67 % | Largest target; not yet addressed |
| scan | ~19 KB | < 1 % | After consolidation |
| detector | 306 KB | 9 % | Shared, stored once — already optimal |
| profile | 87 KB | 2.6 % | Potential future target |

---

## Open questions (status 2026-05-21)

| # | Question | Status |
|---|----------|--------|
| 1 | Should `XFELBeam` become a permanent supported dxtbx model? | Unresolved — needs group discussion; architecture now matches maintainer feedback |
| 2 | Is scan-based stills detection the right contract everywhere? | Unresolved — widest blast radius |
| 3 | Tests for `XFELBeam` / scan `"wavelength"` property / `ImageSequence.get_beam(i)` | Not written yet |
| 4 | Downstream consumers (cctbx.xfel, merging) against new output | Not tested yet |
| 5 | `combine_all_ranks` output semantics | Not settled |
| 6 | Performance benchmarking (file size + runtime) | Not done yet |

---

## How it would ship (if the prototype validates)

1. **dxtbx PR** — `XFELBeam` C++ model (lightweight guard type) + `ImageSequence.get_beam(index)` Python override + `FormatXFEL` mixin + `FormatNXmxXFEL`/`FormatXTCXFEL` + `BeamFactory.make_xfel_beam()` + `CachedWavelengthBeamFactory.get_wavelengths()` + `compact_stills_scans` scan consolidation (`ExperimentList.to_dict/as_json` + `ExperimentListDict._expand_consolidated_scans`). Must land first.
2. **dials PR 1 (correctness)** — `do_import()` pivot + scan-based stills detection everywhere + `_rebuild_shared_imageset_output()` + `_combine_multiprocessing_outputs()` + `compact_stills_scans=True` in composite output.
3. **dials PR 2 (spotfinder caching)** — `ExtractPixelsFromImage` static/dynamic split + `update_imageset()` + `ExtractSpots`/`SpotFinder` caching.
4. **dials PR 3 (indexer caching)** — `Processor` indexer cache + `update_indexer()` + `hardcoded_phil` move + deepcopy removal.

PRs 2 and 3 depend on PR 1 (stable imageset identity is what makes caching valid).

---

## Quick file reference

### dxtbx
| File | What's there |
|------|-------------|
| `src/dxtbx/model/beam.h` | `XFELBeam` C++ class (lightweight guard, no wavelength array) |
| `src/dxtbx/model/boost_python/beam.cc` | Python bindings + pickle + to_dict/from_dict |
| `src/dxtbx/model/beam.py` | `BeamFactory.make_xfel_beam(direction, ...)`, routing in `from_dict()` |
| `src/dxtbx/model/__init__.py` | `get_monochromatic_beam(wavelength)` injection on `XFELBeam`; `ExperimentList.to_dict/as_json` `compact_stills_scans` flag |
| `src/dxtbx/imageset.py` | `ImageSequence.get_beam(index)` override (XFELBeam dispatch) |
| `src/dxtbx/format/FormatXFEL.py` | `FormatXFEL` mixin — builds ImageSequence with XFELBeam + scan wavelength property; wavelength subset fix for composite output |
| `src/dxtbx/format/FormatMultiImage.py` | `get_imageset` model-fill guards (`and format_instance is not None`) — fixes crash when `check_format=False` and goniometer=None |
| `src/dxtbx/format/FormatNXmx.py` | `FormatNXmxXFEL` class |
| `src/dxtbx/format/FormatXTC.py` | `FormatXTCXFEL` class |
| `src/dxtbx/nexus/__init__.py` | `CachedWavelengthBeamFactory.get_wavelengths()` |
| `src/dxtbx/model/experiment_list.py` | `ExperimentListDict._expand_consolidated_scans()` — expands compact scan format on load; `decode()` preserves imageset scan with wavelength property |

### dials
| File | What's there |
|------|-------------|
| `src/dials/command_line/stills_process.py` | Everything: do_import, helpers, caching, output rebuild |
| `src/dials/algorithms/spot_finding/finder.py` | Static/dynamic mask split, update_imageset, SpotFinder cache |
| `src/dials/algorithms/spot_finding/factory.py` | `is_stills` param for SpotFinderFactory |
| `src/dials/algorithms/indexing/indexer.py` | Stills detection refactor, deleted ImageSet conversion |
| `src/dials/algorithms/indexing/stills_indexer.py` | hardcoded_phil cache, deepcopy removal, experiment loop fix |
| `src/dials/algorithms/integration/integrator.py` | Memory calc stills detection fix |
| `src/dials/algorithms/spot_prediction/reflection_predictor.py` | Stills dispatch fix |

### dials/combine_experiments
| File | What's there |
|------|-------------|
| `src/dials/util/combine_experiments.py` | `CombineWithReference.__call__` scan logic — stills keep per-frame scan |
| `src/dials/command_line/combine_experiments.py` | `_save_only_experiments`, `_save_experiments_and_reflections` — `compact_stills_scans=True` on write; `_sort_experiments_and_reflections` rewired to sort stills automatically |

### dials/util/image_viewer
| File | What's there |
|------|-------------|
| `src/dials/util/image_viewer/spotfinder_frame.py` | `_identifiers_for_frame` (order-independent dict); `_still_scan_frame_offset` (derives imageset start offset from per-experiment scans); `get_spotfinder_data` frame-offset fix; fixed per-overlay colours for `viewing_still_scans` |
| `src/dials/util/image_viewer/slip_viewer/frame.py` | `chooser_wrapper.get_beam()` → `image_set.get_beam(self.index)` (per-frame beam); `get_beam_center_px` and resolution-rings use same pattern |

---

## image_viewer: reflection overlays for composite still output (2026-05-22)

**Commit:** `875d1a8cc`

Composite `.expt`/`.refl` output from this branch uses a shared `ImageSequence` for N stills,
each experiment carrying a per-experiment still scan. When loaded in `dials.image_viewer`,
the reflection overlays (spots, centres of mass, predictions, hkl, shoeboxes) appeared on
frame 0 only; all later frames were blank. Two independent root causes:

### Root cause 1: `_identifiers_for_frame` binary search on unsorted list

`_identifiers_for_frame` maps an absolute frame index → experiment identifiers whose
per-experiment scan starts at that frame. The original implementation binary-searched the
experiment list, with a docstring asserting "sorted ascending by frame index, guaranteed by
the shared ImageSequence ordering." This was false:

- `_rebuild_shared_imageset_output` preserves input order, does not sort by frame.
- `_combine_multiprocessing_outputs` appends worker composites end-to-end → frames are
  interleaved (worker 0's strided frames, then worker 1's…), not globally ascending.

Binary search on an unsorted list gives wrong/empty results for every frame except the first.
`identifiers=[]` → `select_on_experiment_identifiers([])` → no rows selected → no overlay.

**Fix:** replaced binary search with an order-independent `{frame_index: [identifiers]}` dict
built by scanning all experiments:

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

### Root cause 2: imageset scan loses start offset on JSON round-trip

`_rebuild_shared_imageset_output` builds the shared imageset with `Scan((min_f+1, max_f+1))`
(e.g., `Scan((6, 1449))` for frames 5–1447). After JSON round-trip, the imageset's own scan
comes back as `Scan((0, 1443))` — offset discarded. The per-experiment scans (individually
serialized) keep the correct absolute frame indices through round-trip.

The frame computation in `get_spotfinder_data` was:
```python
i_frame += imageset.get_scan().get_array_range()[0]  # 0, not 5
```
so chooser position 0 → `i_frame=0`, but per-experiment scans start at 5 → `identifiers=[]`.

**Fix:** new helper `_still_scan_frame_offset(imageset)` derives the absolute offset from the
minimum `scan.get_array_range()[0]` across all experiments that share the imageset. Cached per
imageset by `id()`. `get_spotfinder_data` uses this offset for the `viewing_still_scans` path:
```python
if self.viewing_still_scans:
    i_frame += self._still_scan_frame_offset(imageset)
```

### Root cause 3: `prediction_colours` IndexError / frame-varying colors

`prediction_colours` has 90 entries. With 764 experiments, `reflection["id"]` ranges 0–763 and
crashes `prediction_colours[reflection["id"]]`. In the old `viewing_stills` path, `id` was reset
to 0 per frame, so it never exceeded bounds; in the new `viewing_still_scans` path it is the
composite-table id.

**Fix:** for `viewing_still_scans`, integrated shoeboxes use a fixed purple (`#984ea3`) and
predictions a fixed green (`#4daf4a`), independent of experiment id. Non-still modes keep the
original `prediction_colours[id]` behaviour. The `all_pix` colour sites keep a modulo guard as
crash protection for multi-lattice frames.

### Performance note (not implemented)

The per-frame overlay computation scales poorly with large composite files:

| Path | Complexity | Test case |
|------|-----------|-----------|
| bbox / ctr-mass / shoebox | O(N_refl) full scan per flip | 484,603 reflections |
| predictions | O(N_expt × N_refl) per flip | 764 × 484,603 ≈ 3.7 × 10^8 |

A precomputed `{abs_frame: flex.size_t(row_indices)}` dict built at load time would reduce
both paths to O(rows-on-frame) per flip (~630 rows for the test case). Not implemented.

---

## combine_experiments: stills-as-ImageSequence support (2026-05-22)

**Commits:** dxtbx `e3133a6e`, dials `9494e3a0b`

`dials.combine_experiments` crashed and then lost per-frame wavelengths when processing
integrated output from this branch. Four independent root causes across two repos.

### Root cause 1: `FormatMultiImage.get_imageset` crash (dxtbx)

`combine_experiments` uses `check_format=False`, so `format_instance=None` inside
`FormatMultiImage.get_imageset`. Stills experiments have `goniometer=None` (no goniometer
in the JSON). The model-fill block called `format_instance.get_goniometer()` unconditionally
when `goniometer is None`, causing `AttributeError: 'NoneType' object has no attribute
'get_goniometer'`.

**Fix:** added `and format_instance is not None` to all four model-fill conditions (beam,
detector, goniometer, scan) at `FormatMultiImage.py` ~line 300. This matches the identical
guard already present in `Format.get_imageset` (added in an earlier commit).

### Root cause 2: `FormatXFEL.get_imageset` wavelength mismatch (dxtbx)

When loading a composite stills `.expt` file, the imageset covers only a subset of the
source file's frames. `get_wavelengths()` returns wavelengths for the full file, but the
imageset length is the subset size — the `Scan` was built with the full wavelength list,
making its property length inconsistent with the imageset length.

**Fix:** `FormatXFEL.get_imageset` now checks `len(indices) != len(wavelengths_all)` and,
when they differ, selects `wavelengths = [wavelengths_all[i] for i in indices]`.

### Root cause 3: `ExperimentListDict.decode` scan override (dxtbx)

When decoding an XFEL `.expt` file with `check_format=True`, the format class sets up the
imageset scan with per-frame wavelengths. The decode loop then unconditionally called
`imageset.set_scan(scan)` with the JSON-decoded scan (which may lack the full wavelength
array), discarding the format-provided data.

**Fix:** `experiment_list.py` decode now skips `set_scan` if the imageset already has a scan
with a `"wavelength"` property.

### Root cause 4: `CombineWithReference` scan clobber (dials)

`reference_from_experiment.scan=0` caused every experiment to receive experiment 0's
single-frame scan (containing only experiment 0's wavelength). All other per-frame
wavelengths were lost.

**Fix:** `CombineWithReference.__call__` in `dials/util/combine_experiments.py`: when both
`ref_scan.is_still()` and `experiment.scan.is_still()`, the experiment keeps its own scan.
Rotation scan behaviour is unchanged. The `compact_stills_scans=True` option is now passed
to all four `as_file` call sites in `_save_only_experiments` and
`_save_experiments_and_reflections`, so the N per-frame scans are consolidated into one
compact JSON object on write (safe no-op for rotation data).

### Correct invocation

```bash
dials.combine_experiments *_integrated.expt *_integrated.refl \
    reference_from_experiment.beam=0 \
    reference_from_experiment.detector=0 \
    reference_from_experiment.scan=0
```

`reference_from_experiment.scan=0` is accepted without error but is effectively a no-op for
stills — each experiment retains its own scan. The compact output gives ONE scan JSON object
with all per-frame wavelengths, matching the structure written by `stills_process`.
Do **not** use this option expecting it to copy experiment 0's scan verbatim for stills.

# Development Context: `stills_process_imagesequence` branch

**Branch:** `stills_process_imagesequence` in both dials and dxtbx  
**Status as of 2026-05-21:** Architecture refactored per upstream dxtbx maintainer review (Waterman, McDonagh). `XFELImageSequence` removed; per-frame wavelengths moved to scan properties table. `XFELBeam` retained as lightweight guard/marker type (no wavelength array). Build passes; all smoke tests pass.  
**2026-05-21 (scan consolidation):** `ExperimentList.to_dict/as_json` consolidates stills scans automatically; stills composite output writes one consolidated scan JSON object instead of N, saving ~8 % of file size (−266 KB on a 764-experiment test case). The `compact_stills_scans` option was later removed (2026-05-22) and consolidation made unconditional.  
**2026-05-21 (bug fixes):** Per-frame scan rebuild in `do_import` extended to all `ImageSequence` slices (not only XFEL stills) — rotation CBF frames processed via `convert_sequences_to_stills` also carry absolute scan ranges that would confuse the integrator. `FormatXTC.understand()` fixed to `return bool(ds)` instead of `True`. `Format.get_imageset` goniometer assertion relaxed for still sequences. `test_pseudo_scan` fixed and strengthened with round-trip assertions.  
**2026-05-22 (image_viewer):** Reflection overlays (spots, centres of mass, predictions, hkl, integrated shoeboxes) now draw correctly on all frames of composite still output in `dials.image_viewer`. Two independent root causes fixed; per-frame beam corrected for XFEL. Commit `875d1a8cc`. See "image_viewer changes" section below.  
**2026-05-22 (combine_experiments):** `dials.combine_experiments` now works with stills-as-ImageSequence output. Four fixes across both repos; see "combine_experiments changes" section below.
**2026-05-22 (frame ordering):** `stills_process` and `dials.combine_experiments` composite output now writes experiments (and paired reflections) in ascending frame order. New helper `_sort_experiments_by_frame()` sorts by `(source path, frame index)` and remaps reflection-table `id` and `experiment_identifiers` to match. Called in both the multiprocessing path (`_combine_multiprocessing_outputs`, after `_rebuild_shared_imageset_output`) and the MPI composite-stride path (`finalize`, after the rebuild loop, for all four stages). Previously, experiments were written in worker/rank encounter order.  
**2026-05-22 (scan consolidation made unconditional):** Removed `compact_stills_scans` parameter from `ExperimentList.to_dict()` and `as_json()`. Consolidation now runs automatically whenever all experiments are single-frame stills (oscillation == 0); rotation data is unaffected by the guard. Eliminates the need for ~200+ downstream call sites across dials/cctbx_project to opt in. Paired cleanup in dials drops the 8 now-redundant `compact_stills_scans=True` kwargs from `stills_process` and `combine_experiments`. dxtbx commit `e2d34c52`; dials commit `dd8bfea5c`.  
**2026-05-22 (sparse imageset decode crash fix):** `cctbx.xfel.filter_experiments_by_rmsd` and any downstream reader crashed with `AssertionError` on composite still `.expt` files written by this branch. Root cause: the `eobj_scan` merger in `ExperimentListDict.decode()` builds a merged scan spanning `(min_fn, max_fn)` from the per-experiment scans (after `_expand_consolidated_scans` gives them absolute frame numbers). For sparse output this merged range >> N frames, failing `FormatMultiImage`'s assertion `scan.get_num_images() <= num_images`. Fixed in `_make_sequence` (compact scan passed to `get_imageset`) and `_imageset_from_imageset_data` (guard extended to not overwrite the compact scan with the bad merged scan). dxtbx commit `1467eafe`.

**2026-05-22 (sparse imageset / image_viewer frame filter):** `_rebuild_shared_imageset_output` now builds the shared `ImageSequence` with only the integrated frame indices (`sorted(set(frame_indices))`) and a compact scan `(1, N)`, instead of a contiguous range `[min_f, max_f]`. This means `image_viewer` only shows frames that were successfully integrated. Three code paths had to be updated to support non-contiguous `single_file_indices` for still scan sequences: the C++ `ImageSequence` constructor (`imageset.h`), `FormatMultiImage.get_imageset()`, and `ExperimentListDict._make_sequence()` (bypasses `ImageSetFactory.make_sequence()` when `single_file_indices` is in the JSON). dxtbx commit `eea48ad7`. The `image_viewer` frame mapping for `viewing_still_scans` was updated to use `_still_scan_frame_list(imageset)[idx]` (a sorted list of absolute frame numbers from per-experiment scans) instead of `idx + offset`; falls back to offset for old-style contiguous `.expt` files. dials commit `a7b272a4f`.

**2026-05-22 (dials.refine):** `dials.refine` now runs on composite stills `.expt` output. `StillsReflectionManager` crashed on `XFELBeam.get_s0()`. New `Experiment.get_monochromatic_beam()` accessor combines the shared `XFELBeam` with the per-frame scan `"wavelength"` into a monochromatic `Beam`; refinement resolves it once at `_copy_experiments_for_refining` so no refinement-internal code changed, and `ExperimentList.to_dict()` re-consolidates the per-frame beams back to one `XFELBeam` on write. A second, independent bug also fixed: `_parameterise_crystals` treated the zero-oscillation still scan as a rotation experiment, raising a false "mixture of scan and still experiments" error. Non-refinement tools (e.g. `cctbx.xfel.detector_residuals`) that load composite stills resolve beams with a one-line loop at load. dxtbx commit `cd5d39e8`; dials commit `272a1b469`; cctbx_project commit `a34f7ddc`. See "dials.refine" section below.

**2026-05-24 (review hardening):** Branch review pass; P0 correctness + P1 blast-radius reductions. P0: (a) `Indexer.refined_experiments` is now nulled at the start of `index()`, so a frame that fails on a cached `Indexer` no longer silently inherits the previous frame's refined result; (b) `BeamFactory.make_xfel_beam` and `FormatXFEL.get_imageset` now pass the source beam's polarization, flux, transmission, probe and sample-to-source distance through to the per-frame monochromatic `Beam`, fixing a Lorentz-polarization regression (LP correction was using the default `polarization_fraction=0.5` instead of LCLS's ~0.9); (c) `ExperimentList.decode()` auto-resolves `Experiment.beam` from `XFELBeam` to a per-frame monochromatic `Beam` after `_expand_consolidated_scans`, so downstream tools (`cctbx.xfel.detector_residuals`, `cctbx.xfel.merge`, `dials.show`, etc.) no longer need the per-tool one-line fix — the imageset's `XFELBeam` and the scan's `"wavelength"` property are preserved unchanged. P1: route all open-coded stills detection through `Experiment.is_still()`; fix `spot_finding/factory.py` contract miss; align `combine_experiments`' shared imageset to the sparse layout used by `stills_process`; factor the JSON consolidate/expand helpers shared by `expt_set_to_sequence` / `expt_sequence_to_set` into `dials.util.stills_imageset_convert`; collapse `_still_scan_frame_offset` onto `_still_scan_frame_list[0]`; rewrite `_combine_multiprocessing_outputs` with `*.tmp` + `os.replace` (atomic) and logged cleanup failures; `XFELBeam.get_num_scan_points` returns 0 instead of throwing; `XFELBeam.from_dict` default `polarization_fraction` aligned to 0.5 (matching the C++ ctor); `FormatNXmxXFEL.understand()` wrapped in try/except; `Format.py` goniometer guard restored to `assert` idiom; pinned-`ruff format` pass to undo unrelated black churn. See "Review hardening" section below.

**2026-05-26 (consolidated scan oscillation unit-conversion fix):** `ExperimentList.to_dict()`'s stills-scan consolidation was writing `oscillation`/`oscillation_width` to JSON in **radians** (raw values from the C++ properties table), but every other writer/reader follows a **degrees** JSON convention (`to_dict<Scan>` in `dxtbx/src/dxtbx/model/boost_python/scan.cc:218–227` calls `get_oscillation_arr_in_deg()` explicitly; the reverse `extract_properties_table(..., convert_oscillation_to_rad=true)` at `scan.cc:60–88` does deg→rad on load). The mismatched units caused a "double-conversion" round-trip: `10° → 0.17453 rad` stored → `0.17453` written to JSON → read back as `0.17453°` → final `get_oscillation()` returns `0.17453`. The bug was latent because every prior test still scan used `oscillation=(0,0)` (identical in any unit) — caught by two new `dials.import` tests that use `scan.oscillation=10,0`. Fix in `dxtbx/src/dxtbx/model/__init__.py` consolidation loop: special-case the two angle keys to write `s.get_oscillation(deg=True)[0]` / `[1]` per scan instead of the raw `get_properties()` entry. New regression tests added in `dxtbx/tests/model/test_experiment_list.py` (`test_stills_consolidated_scan_oscillation_roundtrip` and `test_stills_consolidated_rotation_not_consolidated`) — closes open-questions item #7. No read-side change required.

**2026-05-26 (test fixtures aligned to ImageSequence contract):** Two test fixtures that pre-dated this branch built bare `Experiment` objects with no imageset. The branch's auto-routed sort path in `dials.combine_experiments` (`_sort_experiments_and_reflections`, triggered automatically when `expts[0].is_still()`) calls `iset.paths()[0]` on the sort key and crashed with `AttributeError: 'NoneType' object has no attribute 'paths'`. Per established policy on this branch — *don't weaken upstream production code to accommodate malformed test data* — both fixtures were updated, not the sort path. (a) `tests/algorithms/profile_model/ellipsoid/conftest.py` (commit `ae871b468`, 2026-05-21): `test_experiment` fixture rebuilt to construct a real `ImageSequence` with `Format.Reader`, beam/detector/goniometer/scan, then build `Experiment` from the imageset's models. (b) `tests/command_line/test_combine_experiments.py::test_min_max_reflections_per_experiment` (commit `963d588ea`, 2026-05-26): new `_attach_still_imagesets()` helper rewrites the `multi_stills_combined.json` dials_data fixture before invoking the subprocess, attaching one shared zero-oscillation `ImageSequence` (dummy reader) covering N frames with each experiment carrying its own single-frame per-experiment scan. The shared-object trick keeps `_consolidate_stills_imagesets` a no-op so it never tries to open `dummy.cbf`. The rest of the test (subprocess invocation, expected experiment counts, parametrize matrix) is unchanged. See "Testing: fixture alignment" section in `stills_process_imagesequence.md` for full pattern.

**2026-05-24 (split_experiments → combine_experiments round-trip):** `dials.split_experiments` was writing each single-experiment `.expt` with a monochromatic `Beam` rather than the canonical `XFELBeam`; round-tripping the splits back through `dials.combine_experiments reference_from_experiment.beam=0 detector=0 scan=0` then failed the beam tolerance check (per-frame wavelengths visibly differ). Two fixes: (a) `ExperimentList.to_dict()` beam re-consolidation no longer requires `len(beam_models) > 1` and is no longer nested under the scan-consolidation block — it fires whenever every experiment is a stills with a `"wavelength"` scan property, covering both the single-experiment-split case and the single-shared-beam-after-combine case. The `make_xfel_beam` call also now passes polarization/flux/transmission/probe/sample_to_source_distance through (matched the existing P0.2 load-side passthrough — without this, the write side silently dropped the LCLS ~0.999 polarization fraction); (b) `CombineWithReference.__call__` relaxes `wavelength_tolerance` to infinity when both `ref_scan` and `experiment.scan` are stills — after the decode auto-resolve every experiment's `Beam.wavelength` is its own per-frame value, so a strict tolerance was guaranteed to fail; geometry (direction, polarization) is still compared. Verified: 764-experiment composite split + combine reproduces the original consolidated file (beam, scan, imageset, detector match exactly; crystal differs only at ~1e-14 round-off). See "split_experiments round-trip" section below.

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
- `_combine_multiprocessing_outputs()` (~line 615): Merges per-worker composites, calls `_rebuild_shared_imageset_output()`, writes the combined output via `*.tmp` + `os.replace` (atomic) so a crash mid-write cannot corrupt worker-0's filename; logs cleanup `OSError` on workers 1..N-1 rather than swallowing.
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

Every place this changes is a potential regression site if any code path was missed. Touched: indexer, integrator, spotfinder factory, reflection predictor, refinement crystal parameterisation (`configure.py`), stills_process itself.

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

### Solution: automatic stills scan consolidation

**Write side — `ExperimentList.to_dict/as_json` in `src/dxtbx/model/__init__.py`**

Whenever every scan in the experiment list is a single-frame still
(`image_range[0] == image_range[1]`, oscillation width == 0), consolidation
runs automatically (no flag required):

- The N scan dicts are replaced by **one** consolidated object:
  ```json
  {"__stills_consolidated": true, "batch_offset": 0,
   "frame_numbers": [575, 821, ...],
   "properties": {"wavelength": [1.3087, 1.3078, ...]},
   "valid_image_ranges": {}}
  ```
- All-zero property arrays (epochs, exposure_time, oscillation,
  oscillation_width) are omitted; they are reconstructed as zeros on load.
- Angle-valued property arrays (`oscillation`, `oscillation_width`) are
  written in **degrees**, matching the standard `to_dict<Scan>` JSON
  convention. The C++ properties table stores them in radians, so the
  consolidation loop explicitly converts via `s.get_oscillation(deg=True)`
  rather than reading `get_properties()` raw — the read side's
  `extract_properties_table(..., convert_oscillation_to_rad=true)` then
  round-trips them back to radians cleanly.
- Each experiment gains `"scan_point": i` and `"scan": 0`.
- Rotation data (oscillation width > 0) is never consolidated — the guard
  condition fails and the standard per-scan format is written unchanged.

**Read side — `ExperimentListDict._expand_consolidated_scans()` in `src/dxtbx/model/experiment_list.py`**

Called in `__init__` before `_extract_models` runs. Expands the consolidated
object back to N per-frame scan dicts (no-op for files without
`__stills_consolidated`). The rest of `decode()` sees the standard format
unchanged.

### Numbers (764-experiment test case)

| | Size |
|---|---|
| Original | 3,269 KB |
| Compact | 3,002 KB |
| Savings | **−267 KB (8.1 %)** |

### Backward compatibility

- New reader + old file: **transparent** (`_expand_consolidated_scans` exits immediately).
- Old reader + new file: **fails** (`ScanFactory.from_dict` gets a dict with no `image_range`). The dxtbx reader (`_expand_consolidated_scans`) must land before any code that writes stills experiments reaches users.

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
| 2 | Is scan-based stills detection the right contract everywhere? | Resolved (2026-05-24) — all open-coded variants now route through `Experiment.is_still()` (the C++ accessor that already encodes the scan-is-still / scan-is-None / goniometer-is-None contract) |
| 3 | Tests for `XFELBeam` / scan `"wavelength"` property / `ImageSequence.get_beam(i)` | Deferred — tests will be added once the approach is agreed upon and stable |
| 4 | Downstream consumers (cctbx.xfel, merging) against new output | Mitigated (2026-05-24) — `ExperimentList.decode()` now auto-resolves `Experiment.beam` from `XFELBeam` → per-frame `Beam`, so most tools work transparently; still verify against cctbx.xfel.merge and the XFEL GUI |
| 5 | `combine_all_ranks` output semantics | Not settled |
| 6 | Performance benchmarking (file size + runtime) | Not done yet |
| 7 | Dedicated `to_dict()` stills scan consolidation test | Resolved (2026-05-26) — `test_stills_consolidated_scan_oscillation_roundtrip` (non-zero oscillation round-trips in degrees) and `test_stills_consolidated_rotation_not_consolidated` (rotation scans never consolidated) added in `dxtbx/tests/model/test_experiment_list.py`. Added while fixing the radians/degrees unit bug the consolidation write side had been carrying — would have caught it. |

---

## How it would ship (if the prototype validates)

1. **dxtbx PR** — `XFELBeam` C++ model (lightweight guard type, with `get_num_scan_points` returning 0 for transparent generic-code handling) + `ImageSequence.get_beam(index)` Python override + `FormatXFEL` mixin (with full polarization passthrough) + `FormatNXmxXFEL`/`FormatXTCXFEL` + `BeamFactory.make_xfel_beam()` (full polarization/flux/transmission/probe arguments) + `CachedWavelengthBeamFactory.get_wavelengths()` + automatic stills scan consolidation in `ExperimentList.to_dict/as_json` + `ExperimentListDict._expand_consolidated_scans` + `ExperimentList.decode()` auto-resolve of `Experiment.beam` from XFELBeam → per-frame monochromatic Beam (preserves imageset XFELBeam + scan wavelengths, so the resolve is transparent to downstream consumers). Must land first.
2. **dials PR 1 (correctness)** — `do_import()` pivot + scan-based stills detection everywhere + `_rebuild_shared_imageset_output()` + `_combine_multiprocessing_outputs()`. Stills scan consolidation is now automatic (no call-site flag needed).
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
| `src/dxtbx/model/beam.py` | `BeamFactory.make_xfel_beam(direction, divergence, sigma_divergence, polarization_normal, polarization_fraction, flux, transmission, probe, sample_to_source_distance)`, routing in `from_dict()` |
| `src/dxtbx/model/__init__.py` | `get_monochromatic_beam(wavelength)` injection on `XFELBeam`; `Experiment.get_monochromatic_beam()` injection; `ExperimentList.to_dict/as_json` automatic stills scan consolidation + beam re-consolidation |
| `src/dxtbx/imageset.h` | `ImageSequence` C++ constructor — sequential-indices assertion conditioned on `!scan->is_still()` |
| `src/dxtbx/imageset.py` | `ImageSequence.get_beam(index)` override (XFELBeam dispatch) |
| `src/dxtbx/format/FormatXFEL.py` | `FormatXFEL` mixin — builds ImageSequence with XFELBeam + scan wavelength property; wavelength subset fix for composite output; passes source-beam polarization/flux/transmission/probe/sample_to_source_distance through to `BeamFactory.make_xfel_beam` |
| `src/dxtbx/format/FormatMultiImage.py` | `get_imageset` model-fill guards; sequential-indices assertion conditioned on `not scan.is_still()` for sequences |
| `src/dxtbx/format/FormatNXmx.py` | `FormatNXmxXFEL` class |
| `src/dxtbx/format/FormatXTC.py` | `FormatXTCXFEL` class |
| `src/dxtbx/nexus/__init__.py` | `CachedWavelengthBeamFactory.get_wavelengths()` |
| `src/dxtbx/model/experiment_list.py` | `_expand_consolidated_scans()` — expands compact scan format on load; `decode()` preserves imageset scan with wavelength property; `decode()` auto-resolves `Experiment.beam` from `XFELBeam` → per-frame monochromatic `Beam` (imageset XFELBeam slot and scan wavelength property left intact); `_make_sequence()` bypasses `make_sequence()` when `single_file_indices` stored in JSON (supports sparse round-trip) |

### dials
| File | What's there |
|------|-------------|
| `src/dials/command_line/stills_process.py` | Everything: do_import, helpers, caching, output rebuild |
| `src/dials/algorithms/spot_finding/finder.py` | Static/dynamic mask split, update_imageset, SpotFinder cache |
| `src/dials/algorithms/spot_finding/factory.py` | `is_stills` param for SpotFinderFactory |
| `src/dials/algorithms/indexing/indexer.py` | Stills detection refactor, deleted ImageSet conversion; `index()` resets `self.refined_experiments` at entry (so a frame failure on the cached Indexer cannot inherit the previous frame's result); call sites routed through `Experiment.is_still()` |
| `src/dials/util/stills_imageset_convert.py` | **New.** Shared `.expt` JSON helpers: `source_file`, `frame_index`, `build_consolidated_scan`, `expand_consolidated_scan`, `per_frame_scan_map`.  Used by `expt_set_to_sequence` / `expt_sequence_to_set` to eliminate the divergent third copy of the consolidate/expand logic |
| `src/dials/algorithms/indexing/stills_indexer.py` | hardcoded_phil cache, deepcopy removal, experiment loop fix |
| `src/dials/algorithms/integration/integrator.py` | Memory calc stills detection fix |
| `src/dials/algorithms/spot_prediction/reflection_predictor.py` | Stills dispatch fix |
| `src/dials/algorithms/refinement/refiner.py` | `_copy_experiments_for_refining()` resolves per-frame `XFELBeam` via `Experiment.get_monochromatic_beam()` — single chokepoint for refinement |
| `src/dials/algorithms/refinement/parameterisation/configure.py` | `_parameterise_crystals()` — still scan no longer counts as a rotation experiment |

### dials/combine_experiments
| File | What's there |
|------|-------------|
| `src/dials/util/combine_experiments.py` | `CombineWithReference.__call__` scan logic — stills keep per-frame scan |
| `src/dials/command_line/combine_experiments.py` | `_save_only_experiments`, `_save_experiments_and_reflections` — stills scan consolidation automatic on write; `_sort_experiments_and_reflections` rewired to sort stills automatically; `_consolidate_stills_imagesets` now builds the shared `ImageSequence` with the sparse layout (`sorted(set(frame_indices))` + compact `Scan((1, N))`) to match `stills_process._rebuild_shared_imageset_output`, so image_viewer's sparse-mapping path applies uniformly |

### dials/format conversion scripts
| File | What's there |
|------|-------------|
| `src/dials/command_line/expt_set_to_sequence.py` | Converts old per-frame ImageSet `.expt` → new ImageSequence format. Groups experiments by source file, deduplicates detector, builds XFELBeam, writes consolidated scan. Sorts experiments by `(src, fi)` ascending; remaps `.refl` `id` column and `experiment_identifiers` to match. Default output: `sequence.expt` / `sequence.refl`. JSON-level helpers (consolidated-scan dict builder, `source_file`, `frame_index`) imported from `dials.util.stills_imageset_convert`. |
| `src/dials/command_line/expt_sequence_to_set.py` | Converts new ImageSequence `.expt` → old per-frame ImageSet format. Expands consolidated scan, creates one ImageSet + monochromatic Beam + Detector copy per experiment. Experiment order is preserved so no `.refl` remapping is needed (refl is copied to output filename unchanged). Default output: `imageset.expt` / `imageset.refl`. Consolidated-scan reader (`expand_consolidated_scan`, `per_frame_scan_map`) imported from `dials.util.stills_imageset_convert`. |

**Usage:**

```bash
# Old → new (with refl remapping)
dials.expt_set_to_sequence integrated.expt integrated.refl \
    output.experiments_filename=sequence.expt \
    output.reflections_filename=sequence.refl

# New → old
dials.expt_sequence_to_set integrated.expt integrated.refl \
    output.experiments_filename=imageset.expt \
    output.reflections_filename=imageset.refl
```

### dials/util/image_viewer
| File | What's there |
|------|-------------|
| `src/dials/util/image_viewer/spotfinder_frame.py` | `_identifiers_for_frame` (order-independent dict); `_still_scan_frame_list` (sorted absolute frame list for sparse imageset mapping); `_still_scan_frame_offset` is now a one-line wrapper returning `frame_list[0]` (single `id()`-keyed cache); `get_spotfinder_data` frame-offset/list fix; fixed per-overlay colours for `viewing_still_scans` |
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
Rotation scan behaviour is unchanged. The N per-frame scans are automatically consolidated
into one compact JSON object on write by `ExperimentList.to_dict()` (safe no-op for
rotation data).

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

---

## Sparse imageset / image_viewer frame filter (2026-05-22)

**Commits:** dxtbx `eea48ad7`, dials `dd8bfea5c` (stills_process), `a7b272a4f` (image_viewer)

`dials.image_viewer` was showing all frames from `min_f` to `max_f` in the composite
output, including frames that were not successfully integrated. Root cause:
`_rebuild_shared_imageset_output` built the shared `ImageSequence` with
`indices = range(min_f, max_f + 1)` — a contiguous slice covering every frame in the
range, regardless of integration success.

### Fix: sparse imageset indices

`_rebuild_shared_imageset_output` now uses:
```python
sorted_frames = sorted(set(frame_indices))  # only integrated frames
full_seq = _make_stills_sequence(
    parent_iset.data(),
    flex.size_t(sorted_frames),
    shared_beam,
    shared_detector,
    (1, len(sorted_frames)),              # compact scan
)
```

The shared imageset has exactly N entries (one per unique integrated frame), so
`image_viewer` sees only N chooser positions — one per integrated frame.

### Three assertions blocking sparse `ImageSequence` round-trip

`single_file_indices` is already written to the JSON for single-file readers. Three
code paths required changes to accept non-contiguous indices:

**1. C++ `ImageSequence` constructor (`imageset.h` ~line 1278)**  
The sequential-indices loop is now conditioned on `!scan->is_still()`. Still scans
(oscillation width == 0) have independent shots; rotation scans still require sequential.

**2. `FormatMultiImage.get_imageset()` (`FormatMultiImage.py` ~line 286)**  
The Python-level sequential assertion is similarly conditioned on `not (scan is not None
and scan.is_still())`.

**3. `ExperimentListDict._make_sequence()` (`experiment_list.py` ~line 628)**  
When `single_file_indices` is stored in the JSON, `_make_sequence` now bypasses
`ImageSetFactory.make_sequence()` entirely and calls `format_class.get_imageset()`
directly. `make_sequence()` recomputes a contiguous range from the scan's `image_range`
and asserts `array_range == scan.get_array_range()`, both of which fail for a sparse
imageset with compact scan `(1, N)`.

**Backward compatibility:** old `.expt` files with contiguous `single_file_indices`
take the same new bypass path and produce an identical result. Template-based files
(no `single_file_indices` in JSON) still go through `make_sequence()` unchanged.

### image_viewer frame mapping fix

With a sparse imageset, chooser position `k` → the `k`-th integrated frame
(in ascending absolute frame order). The old `i_frame = idx + offset` mapping
assumed contiguous frames. New helper `_still_scan_frame_list(imageset)` returns
the sorted list of `expt.scan.get_array_range()[0]` values for experiments sharing
the imageset (cached per imageset by `id()`). `get_spotfinder_data` uses:

```python
frame_list = self._still_scan_frame_list(imageset)
if frame_list and len(frame_list) == len(imageset.indices()):
    i_frame = frame_list[i_frame]   # sparse: direct lookup
else:
    i_frame += self._still_scan_frame_offset(imageset)  # fallback for old files
```

---

## dials.refine: per-frame beam for composite stills (2026-05-22)

**Commits:** dxtbx `cd5d39e8`, dials `272a1b469`, cctbx_project `a34f7ddc`

`dials.refine` crashed on composite stills `.expt` output from this branch:

```
RuntimeError: dxtbx Error: XFELBeam has no fixed s0
```

Each `Experiment.beam` loaded from the file is the shared `XFELBeam` (a guard type
with no wavelength). `StillsReflectionManager.__init__` builds
`self._s0vecs = [matrix.col(e.beam.get_s0()) for e in self._experiments]`
(`reflection_manager.py:545`), and `XFELBeam.get_s0()` throws by design.

### Design constraint

`Experiment.beam` must stay the shared `XFELBeam` on load: the file holds one beam
object, and `dials.combine_experiments ... reference_from_experiment.beam=0` must
remain a true no-op (every experiment genuinely shares the one beam). The per-frame
s0 is therefore derived by combining the beam direction with the frame wavelength
from the scan — the wavelength is not put back on the beam.

### Fix 1: upstream accessor + write-side re-consolidation (dxtbx)

`Experiment.get_monochromatic_beam()` — injected in `model/__init__.py` alongside
`load_models`. When `self.beam` is an `XFELBeam` and `self.scan` carries a
`"wavelength"` property, returns a monochromatic `Beam` built from the `XFELBeam`
direction + `self.scan.get_property("wavelength")[0]`. Returns `self.beam`
unchanged for any non-`XFELBeam`, so callers may use it unconditionally.

`ExperimentList.to_dict()` — once experiments carry per-frame monochromatic beams
(after refinement), the N beams would serialize separately and bloat the file.
Inside the existing stills scan-consolidation block, when the consolidated scan
carries a `"wavelength"` array and every experiment beam is a plain `Beam` sharing
direction / divergence / polarization, the N beams are collapsed to one
reconstructed `XFELBeam` (`BeamFactory.make_xfel_beam`) and every experiment is
pointed at it. The per-frame wavelengths are already written by the scan
consolidation, so nothing is lost. If any beam was genuinely refined to a different
geometry the guard fails and all N beams are written — no data loss.

### Fix 2: resolve once at the refinement chokepoint (dials)

`_copy_experiments_for_refining()` in `refinement/refiner.py` is the single point
through which experiments enter the refiner (and leave it, via
`Refiner.get_experiments()`). After the existing model-copy loop it now sets
`new_exp.beam = new_exp.get_monochromatic_beam()`. Every refinement-internal
consumer — `StillsReflectionManager`, `StillsReflectionPredictor`,
`BeamParameterisation`, the prediction parameterisations, residual helpers — then
sees an ordinary monochromatic `Beam` and required no change. Refined output flows
back through the same helper and is re-consolidated to one `XFELBeam` by
`to_dict()`. A no-op for non-XFEL data.

This deliberately avoids editing the many individual `expt.beam.get_s0()` call
sites across refinement (reflection manager, predictor factories, beam
parameterisation, prediction parameters, residual helpers) — each was confirmed to
crash when tried piecemeal. `decode()` is left unchanged: `Experiment.beam` loaded
from a `.expt` file is still the shared `XFELBeam`, so the on-disk file keeps its
single `XFELBeam` and `combine_experiments reference_from_experiment.beam=0` stays
a true no-op. Only the refinement *copy* carries monochromatic beams.

### Fix 3: a still scan is not a rotation experiment (dials)

Once the s0 crash cleared, refinement failed with "A crystal model appears in a
mixture of scan and still experiments". `_parameterise_crystals()` in
`refinement/parameterisation/configure.py` classified an experiment as a *scan*
experiment whenever `experiment.scan is not None` — but on this branch every still
carries a zero-oscillation scan. The guard now counts a scan as a rotation
experiment only when `not scan.is_still()`. This is another instance of the
scan-based stills-detection contract change; the site had simply not been exercised
until refinement ran against this branch's output.

### Fix 4: non-refinement tools resolve their beams at load (cctbx_project)

The refinement chokepoint only covers `dials.refine`. Other tools that load a
composite stills `.expt` and call `expt.beam.get_s0()` directly hit the same crash.
`cctbx.xfel.detector_residuals` does so both via the shared `StillsReflectionPredictor`
factory and in four direct sites (beam centre, two-theta, per-reflection energy).
The fix is one line right after the experiments are loaded:

```python
for expt in experiments:
    expt.beam = expt.get_monochromatic_beam()
```

A no-op for ordinary beams. It is safe for any tool: if the tool writes a `.expt`,
`ExperimentList.to_dict()` re-consolidates the per-frame beams back to one
`XFELBeam` on the way out. `detector_residuals.run()` applies this immediately
after `flatten_experiments`.

**Update (2026-05-24):** this per-tool one-liner is superseded for `.expt`-loaded
experiments by the `ExperimentList.decode()` auto-resolve added in the review
hardening pass — `decode()` now runs the same `expt.beam = expt.get_monochromatic_beam()`
loop centrally, so any consumer that gets its experiments from
`ExperimentListFactory.from_*` (i.e. by loading a `.expt` file) sees a regular
`Beam` on `Experiment.beam` without modification.  The `cctbx_project` commit
remains valid (still needed for in-memory pipelines that build experiments without
going through `decode()`), but no further cctbx.xfel tools need patching to load
composite stills.

### Verified

726-experiment filtered XFEL `.expt`:
- `dials.refine refine_level0.phil` completes; output `refined_level0.expt` has
  exactly one `"__id__": "xfel"` beam, the scan stays `__stills_consolidated` with
  the 726-value per-frame `wavelength` array, and reload gives `Experiment.beam`
  back as an `XFELBeam` that resolves correctly via `get_monochromatic_beam()`.
- `cctbx.xfel.detector_residuals hierarchy=1` completes and prints full residual
  statistics by panel group.

### Known limitation

If a refinement phil refines beam parameters for XFEL stills, `BeamParameterisation`
operates on the resolved monochromatic beam (a fresh object), and `to_dict()` will
write the N differing beams rather than re-consolidating. Beam refinement for XFEL
stills needs separate design. `refine_level0.phil` fixes the beam (`beam.fix=*all`),
so this does not affect the verified case.

---

## Review hardening (2026-05-24)

**Commits:** dxtbx `64b7c062`, `16294afa`, `266edaf9`, `2498250f`; dials
`054bb1fb1`, `6ee2a8f8e`, `3b4083a9f`, `c02872d72`, `5e4d11e15`, `61dfdab69`.

A review pass on the branch checkpoint produced three P0 correctness fixes and
seven P1 hardening / blast-radius reductions.

### P0.1 Indexer fails-silently across cached frames (dials)

`Processor` caches a single `Indexer` instance across frames (warm-start for
`d_min`, geometry, etc.) and uses `update_indexer()` to reset per-frame inputs.
`self.refined_experiments`, however, was nulled only in `Indexer.__init__`.
Result: a frame whose indexing failed before the success block at the end of
`index()` could not trigger the `if self.refined_experiments is None: raise
DialsIndexRefineError(...)` guard — `self.refined_experiments` was still
holding the previous frame's result, and the failure was silently promoted.

**Fix:** null `self.refined_experiments` at the start of `index()`.  The other
warm-start state (`self.d_min`, detector/crystal geometry) is left intact —
the warm-start intent is for *inputs*, not for *outputs*.

### P0.2 XFEL polarization regression (dxtbx)

`BeamFactory.make_xfel_beam` accepted only `direction`, `divergence`,
`sigma_divergence` — silently dropping the source beam's `polarization_normal`,
`polarization_fraction`, `flux`, `transmission`, `probe`, and
`sample_to_source_distance`.  `XFELBeam.get_monochromatic_beam(wl)` then
constructed the per-frame `Beam` without those values either, so every
per-frame beam used for indexing and integration carried
`polarization_fraction=0.5` (the C++ default) rather than the instrument
value (~0.9 for LCLS).  The Lorentz-polarization correction in integration
uses this directly → systematically wrong intensities.

**Fix:** extend `make_xfel_beam` to accept the full polarization tuple
(delegating to the existing C++ `XFELBeam` full constructor); have
`FormatXFEL.get_imageset` read the values from the source beam and pass them
through; have `XFELBeam.get_monochromatic_beam` forward them when building
the per-frame `Beam`.

**Verified:** `make_xfel_beam(polarization_fraction=0.92)` → `to_dict`/`from_dict`
→ `get_monochromatic_beam(1.5)` returns a `Beam` with `pf=0.92, wavelength=1.5`.

### P0.3 Auto-resolve XFELBeam on Experiment.beam at load time (dxtbx)

Previously, every downstream tool that loaded a composite stills `.expt` and
called `expt.beam.get_s0()` needed a per-tool one-line fix to call
`expt.beam = expt.get_monochromatic_beam()` first (see Fix 4 of the dials.refine
section above).  The blast radius was unbounded across cctbx.xfel and any other
consumer.

**Fix:** `ExperimentList.decode()` now does the auto-resolve once at load time,
right after `_expand_consolidated_scans` has restored the per-frame scans.
For every `Experiment` whose `beam` is an `XFELBeam`, `expt.beam` is replaced
with `expt.get_monochromatic_beam()` — a regular `Beam` with the correct
per-frame wavelength.

**Critically:** the imageset's `XFELBeam` and the scan's `"wavelength"`
property are left untouched.  `imageset.get_beam(i)` still dispatches to
`XFELBeam.get_monochromatic_beam` via the Python injection, so per-frame
wavelength variations remain queryable from the imageset.  The write-side
re-consolidation in `to_dict()` is unaffected: per-frame monochromatic beams
post-refinement still get collapsed back to one `XFELBeam` on write.

After this lands, the per-tool fix in `cctbx.xfel.detector_residuals.py` is
redundant for `.expt`-loaded experiments (still needed for runtime-built
experiments that bypass `decode()`).  The refinement chokepoint in
`_copy_experiments_for_refining` is similarly redundant for the loaded case.

### P1 hardening

- **`Experiment.is_still()` routed everywhere.**  Six open-coded variants of
  `expt.scan is None or expt.scan.is_still()` (and several `isinstance(...,
  ImageSequence)` checks) are now routed through the existing
  `Experiment.is_still()`.  Touches `indexer.py` (3 sites), `spot_finding/factory.py`
  (fixed a contract miss: the `is_stills=False` fallback was still using the
  isinstance check, inverting `no_shoeboxes_2d` for stills entering via
  non-`stills_process` callers), `spot_prediction/reflection_predictor.py`,
  `refinement/parameterisation/configure.py`, and
  `command_line/combine_experiments.py`.
- **`combine_experiments` sparse imageset layout.**  `_consolidate_stills_imagesets`
  was building `ImageSequence` indices as a contiguous
  `range(min_fi, max_fi+1)` with `Scan((min_fi+1, max_fi+1))`.  `stills_process`'
  `_rebuild_shared_imageset_output` was later updated to a sparse layout
  (`sorted(set(frame_indices))` + compact `Scan((1, N))`), but
  `combine_experiments` was missed — so the two producers emitted different
  imageset structures and the image_viewer `len(frame_list) ==
  len(imageset.indices())` check silently fell back to offset mapping for
  combine_experiments output.  `_consolidate_stills_imagesets` now uses the
  sparse layout.
- **Shared utility for `.expt` JSON conversion scripts.**  `expt_set_to_sequence`
  and `expt_sequence_to_set` manipulate `.expt` JSON directly (not via
  `ExperimentList`) and had their own copies of the helpers for reading and
  writing the `__stills_consolidated` scan dict — a third implementation of the
  consolidation logic.  Factored into `dials.util.stills_imageset_convert`:
  `source_file`, `frame_index`, `build_consolidated_scan`,
  `expand_consolidated_scan`, `per_frame_scan_map`.  Both scripts now import
  from it.  Kept as standalone command-line entry points (users rely on them).
- **`image_viewer` frame-helper duplication.**  Removed
  `_still_scan_frame_offset` (a 24-line near-duplicate of
  `_still_scan_frame_list`); kept the latter and have the offset accessor read
  `frame_list[0]`.  Single `id()`-keyed cache instead of two.
- **`_combine_multiprocessing_outputs` atomic write.**  Used to write the merged
  result directly over worker-0's filename and `os.remove` the other workers'
  files while swallowing `OSError`; a crash mid-write could corrupt worker-0 and
  leave no recovery path.  Now writes to a `*.tmp` sibling and renames with
  `os.replace` (atomic); `OSError` on cleanup of workers 1..N-1 is logged
  rather than swallowed.
- **`XFELBeam::get_num_scan_points()` returns 0.**  Throwing here forced generic
  Experiment copy/compare/serialize code (which probes "is this scan-varying?"
  as a benign check) to need XFEL-awareness.  Returning 0 matches the real
  semantics — `XFELBeam` genuinely has no scan points.
- **`XFELBeam.from_dict` default `polarization_fraction` 0.5.**  The Python
  binding default was 0.999 while the C++ ctor and `make_xfel_beam` used 0.5;
  align to 0.5.  This is the case where the field is absent from the dict
  entirely — to_dict always writes it, so existing files round-trip unchanged.
- **`FormatNXmxXFEL.understand()` guarded.**  The HDF5 open + wavelength check
  is now in `try/except Exception: return False`, so an HDF5 error returns
  "not my format" instead of aborting format detection.
- **`Format.py` goniometer guard idiom.**  Restored to
  `assert goniometer is not None or scan.is_still(), ...` from the
  `if/raise AssertionError` form (matches the surrounding assert style).
- **Pinned `ruff format` pass.**  An earlier `black` pass on the branch used
  an older version than the repos' pinned pre-commit; rerunning the pinned
  ruff-format / ruff-check brings the diff back to canonical style and
  removes ~15 hunks of pure formatter churn that were polluting the diff.

### Files changed (summary)

| Repo | File | What's changed |
|------|------|----------------|
| dxtbx | `src/dxtbx/model/beam.py` | `make_xfel_beam` polarization passthrough |
| dxtbx | `src/dxtbx/model/__init__.py` | `XFELBeam.get_monochromatic_beam` polarization passthrough |
| dxtbx | `src/dxtbx/format/FormatXFEL.py` | source-beam polarization → XFELBeam |
| dxtbx | `src/dxtbx/model/beam.h` | `get_num_scan_points` returns 0 |
| dxtbx | `src/dxtbx/model/boost_python/beam.cc` | `from_dict` pf default 0.5 |
| dxtbx | `src/dxtbx/model/experiment_list.py` | `decode()` auto-resolve XFELBeam → Beam |
| dxtbx | `src/dxtbx/format/Format.py` | restore assert idiom |
| dxtbx | `src/dxtbx/format/FormatNXmx.py` | `FormatNXmxXFEL.understand` try/except |
| dials | `src/dials/algorithms/indexing/indexer.py` | reset `refined_experiments`; route through `is_still()` |
| dials | `src/dials/algorithms/spot_finding/factory.py` | fix contract miss via `is_still()` |
| dials | `src/dials/algorithms/spot_prediction/reflection_predictor.py` | route through `is_still()` |
| dials | `src/dials/algorithms/refinement/parameterisation/configure.py` | route through `is_still()` |
| dials | `src/dials/command_line/combine_experiments.py` | sparse layout + `is_still()` |
| dials | `src/dials/util/stills_imageset_convert.py` | **new** — shared JSON helpers |
| dials | `src/dials/command_line/expt_set_to_sequence.py` | use shared helpers |
| dials | `src/dials/command_line/expt_sequence_to_set.py` | use shared helpers |
| dials | `src/dials/util/image_viewer/spotfinder_frame.py` | collapse frame-offset helper |
| dials | `src/dials/command_line/stills_process.py` | atomic combine write |

---

## split_experiments round-trip (2026-05-24)

`dials.split_experiments combined.expt combined.refl` writes one
`ExperimentList([experiment]).as_json(...)` per experiment. Loading the
original 764-experiment composite via `decode()` gives each experiment a
monochromatic `Beam` (after the auto-resolve added in the review hardening
pass). When `as_json` then ran on a single-experiment list, the existing
beam re-consolidation guards in `ExperimentList.to_dict()` blocked the
collapse back to `XFELBeam`:

```python
if (
    "wavelength" in consolidated_props          # requires scan consolidation
    and len(beam_models) > 1                    # excludes 1-experiment writes
    and all(type(b) is Beam for b in beam_models)
):
```

Both guards excluded the 1-experiment case, and the outer block itself was
gated on `len(scan_models) > 1`. The split file therefore serialized
`__id__: monochromatic` with a wavelength on the beam, instead of
`__id__: xfel` with the wavelength on the scan.

Running `dials.combine_experiments split_*.expt split_*.refl
reference_from_experiment.beam=0 detector=0 scan=0` then hit a separate
failure: after decode's auto-resolve every loaded `Experiment.beam` carries
its own per-frame wavelength, so `BeamComparison(wavelength_tolerance=1e-6)`
fails for every experiment after experiment 0 (e.g. 1.308481 vs 1.307819).

### Fix 1: write-side beam re-consolidation (dxtbx)

`ExperimentList.to_dict()` in `src/dxtbx/model/__init__.py` now computes a
single `all_stills` flag (every scan single-frame still, `frame > 0`,
oscillation == 0). Scan consolidation still gates on `len(scan_models) > 1`
(no JSON change for the rotation or 1-experiment-rotation cases). Beam
re-consolidation moved out of the scan block and runs whenever:

- `all_stills` holds, **and**
- every scan carries a `"wavelength"` property (so the per-frame info is
  preserved without the beam-side wavelength), **and**
- every `beam_models` entry is `type(b) is Beam` and they share geometry
  (`sample_to_source_direction`, `divergence`, `sigma_divergence`,
  `polarization_normal`, `polarization_fraction`).

The `BeamFactory.make_xfel_beam` call also now passes
`polarization_normal`, `polarization_fraction`, `flux`, `transmission`,
`probe`, `sample_to_source_distance` through (it previously only passed
the geometric three, silently dropping the LCLS ~0.999 polarization at
re-consolidation — analogous to the P0.2 load-side passthrough fix).

### Fix 2: stills-aware beam comparison (dials)

`CombineWithReference.__call__` in `src/dials/util/combine_experiments.py`
relaxes `wavelength_tolerance` to `float("inf")` when both `self.ref_scan`
and `experiment.scan` are stills. For stills the per-frame wavelength
lives in the scan property, not on the beam — comparing beam wavelengths
is meaningless and the strict 1e-6 default tolerance fails on any
multi-frame XFEL stills input. Direction, polarization-normal and
polarization-fraction tolerances are unchanged.

### Verified

`/global/cfs/cdirs/m4734/users/dwmoreau/refactor/test_cyt_sequence/test`:

```bash
dials.split_experiments ../idx-0000_integrated.expt ../idx-0000_integrated.refl
dials.combine_experiments split_*.expt split_*.refl \
    reference_from_experiment.beam=0 \
    reference_from_experiment.detector=0 \
    reference_from_experiment.scan=0
```

`combined.expt` matches `idx-0000_integrated.expt` on every model relevant
to the round-trip: beam (one `__id__: xfel`, `polarization_fraction:
0.999`), scan (one `__stills_consolidated` object with 764 `frame_numbers`
and 764 `wavelength` entries — bit-identical to the original arrays), one
detector (identical), one imageset (same template, same 764
`single_file_indices`). Crystal vectors differ only at ~1e-14 round-off
(re-serialization noise). Pre-existing: `CombineWithReference.__call__`
does not propagate `experiment.profile` — orthogonal to this fix and not
addressed here.

### Files changed

| Repo | File | What's changed |
|------|------|----------------|
| dxtbx | `src/dxtbx/model/__init__.py` | `to_dict()` beam re-consolidation decoupled from scan consolidation; runs for 1-experiment and 1-shared-beam stills writes; passes full polarization tuple to `make_xfel_beam` |
| dials | `src/dials/util/combine_experiments.py` | `CombineWithReference.__call__` skips wavelength tolerance when both ref and current experiment are stills |

# Development Context: `stills_process_imagesequence` branch

**Branch:** `stills_process_imagesequence` in both dials and dxtbx (cctbx_project has a paired commit).
**Status:** Build passes; smoke tests pass; review-hardening landed; split→combine round-trip verified bit-identical on 764-experiment composite.
**Companions:** `stills_process_imagesequence.md` (deep-dive narrative incl. testing policy), `stills_process_imagesequence_overview.md` (architecture brief). This doc is the operational reference.

Run Python as `libtbx.python` (not bare `python3`).

| Repo | Path |
|------|------|
| dials | `/pscratch/sd/d/dwmoreau/cctbx_05_14_26/cctbx/modules/dials` |
| dxtbx | `/pscratch/sd/d/dwmoreau/cctbx_05_14_26/cctbx/modules/dxtbx` |
| cctbx_project | `/pscratch/sd/d/dwmoreau/cctbx_05_14_26/cctbx/modules/cctbx_project` |

---

## What this branch does

`dials.stills_process` normally converts every imported image into a plain `ImageSet` (one object per frame). This branch keeps data as an `ImageSequence` for the entire pipeline. Two payoffs: (1) experiments from one source file share a single imageset/beam/detector, so `.expt` files shrink; (2) the spotfinder and indexer can be cached across frames because imageset identity is now stable. The pivot is `do_import()` no longer calling `imageset_from_anyset()` — that conversion was destroying object identity and blocking both the JSON deduplication and the caches.

For XFEL data specifically, per-frame wavelengths are stored in `scan.get_property("wavelength")` (the dxtbx scan properties table). The shared beam is an `XFELBeam` — a lightweight guard/marker type that throws on `get_wavelength()` / `get_s0()`. `ImageSequence.get_beam(i)` dispatches to `XFELBeam.get_monochromatic_beam(wl)` to return a per-frame monochromatic `Beam`.

---

## Cross-cutting contract change

**Old:** stills detected by `isinstance(imageset, ImageSequence)` — False for `ImageSet`, True for `ImageSequence`.
**New:** stills detected by `Experiment.is_still()` (the C++ accessor that encodes `scan is None` OR `scan.is_still()` OR `goniometer is None`). Every prior open-coded variant has been routed through this single accessor.

**Touch sites:** `algorithms/indexing/indexer.py` (3 sites incl. `from_parameters` routing and `_xyzcal_mm_to_px` guard), `algorithms/spot_finding/factory.py` (fixed contract miss: `is_stills=False` fallback was using the old isinstance check and inverting `no_shoeboxes_2d` for stills entering via non-`stills_process` callers), `algorithms/spot_prediction/reflection_predictor.py`, `algorithms/refinement/parameterisation/configure.py` (`_parameterise_crystals` — a zero-oscillation still no longer counts as a rotation experiment), `algorithms/integration/integrator.py` (`_determine_max_memory_needed`), `command_line/combine_experiments.py`, `command_line/stills_process.py`.

---

## dxtbx surface

### `XFELBeam` (C++) — `src/dxtbx/model/beam.h`, `src/dxtbx/model/boost_python/beam.cc`

Inherits `Beam`. Lightweight guard/marker — **no** wavelength array (unlike `PolychromaticBeam`). XFEL pulses are individually monochromatic; per-shot energies live in scan properties, not on the beam.

- Guards throw: `get_wavelength()`, `set_wavelength()`, `get_s0()`, `set_s0()`, all scan-point methods.
- `get_num_scan_points()` returns **0** (not throws), so generic Experiment copy/compare/serialize code that probes "is this scan-varying?" works without XFEL awareness.
- Constructors: default; `(direction, divergence, sigma_divergence)`; full with `polarization_normal, polarization_fraction, flux, transmission, probe, sample_to_source_distance`.
- `operator==` and `is_similar_to` compare direction + divergence + polarization (no wavelength).
- Serializes with `"__id__": "xfel"`; no `wavelengths` key in dict.
- `from_dict` `polarization_fraction` default is **0.5** (matches the C++ ctor; the field is always written by `to_dict`, so existing files round-trip unchanged).

### `BeamFactory.make_xfel_beam` — `src/dxtbx/model/beam.py`

Accepts the **full** polarization tuple (`polarization_normal, polarization_fraction, flux, transmission, probe, sample_to_source_distance`). Both the load-side (`FormatXFEL.get_imageset`) and write-side (`ExperimentList.to_dict` beam re-consolidation) pass these through. Dropping the polarization values silently regressed Lorentz-polarization correction (LCLS ~0.999 fell back to the 0.5 C++ default).

### `XFELBeam.get_monochromatic_beam(wavelength)` — Python injection in `src/dxtbx/model/__init__.py`

Returns a monochromatic `Beam` with the given wavelength, carrying the `XFELBeam`'s direction/divergence/polarization through. Called by `ImageSequence.get_beam(i)` and by `Experiment.get_monochromatic_beam()`.

### `Experiment.get_monochromatic_beam()` — Python injection in `src/dxtbx/model/__init__.py`

When `self.beam` is an `XFELBeam` and `self.scan` carries a `"wavelength"` property, returns a monochromatic `Beam` built from the XFELBeam direction + `self.scan.get_property("wavelength")[0]`. Returns `self.beam` unchanged for any non-`XFELBeam`, so callers may use it unconditionally.

### `ImageSequence.get_beam(index)` Python override — `src/dxtbx/imageset.py`

Replaces the removed `XFELImageSequence` class. When the imageset beam is `XFELBeam` and the scan has a `"wavelength"` property: `get_beam(i)` → `xfel_beam.get_monochromatic_beam(wl[i])`. Index=None or non-XFELBeam: returns the shared beam unchanged.

### Format classes

- **`FormatXFEL` mixin** — `src/dxtbx/format/FormatXFEL.py`. `get_imageset()` builds an `ImageSequence` with `XFELBeam` (no wavelength) + a zero-oscillation `Scan` carrying the `"wavelength"` property. Passes source-beam `polarization_normal/_fraction, flux, transmission, probe, sample_to_source_distance` through to `make_xfel_beam`. **Subset fix:** when the raw imageset covers only a subset of source-file frames (composite stills `.expt` reload), `wavelengths = [wavelengths_all[i] for i in raw_iset.indices()]`.
- **`FormatNXmxXFEL`** — `src/dxtbx/format/FormatNXmx.py`. `get_wavelengths()` delegates to `CachedWavelengthBeamFactory.get_wavelengths()`. `understand()` wrapped in `try/except Exception: return False`.
- **`FormatXTCXFEL`** — `src/dxtbx/format/FormatXTC.py`. `get_wavelengths()` calls `self.get_beam(i).get_wavelength()` per frame. `FormatXTC.understand()` returns `bool(ds)` (previously returned `True` whenever psana was importable — matched arbitrary files).
- **`CachedWavelengthBeamFactory.get_wavelengths()`** — `src/dxtbx/nexus/__init__.py`. Returns `list[float]` Ångström (renamed from `make_xfel_beam`).

### Scan consolidation — write side

`ExperimentList.to_dict()` / `as_json()` in `src/dxtbx/model/__init__.py` automatically consolidate stills scans whenever every scan is a single-frame still (oscillation width == 0). N scan dicts become one object:

```json
{"__stills_consolidated": true, "batch_offset": 0,
 "frame_numbers": [575, 821, ...],
 "properties": {"wavelength": [1.3087, ...]},
 "valid_image_ranges": {}}
```

- All-zero arrays (epochs, exposure_time, oscillation, oscillation_width) are omitted; reconstructed as zeros on load.
- **Angle-valued arrays (`oscillation`, `oscillation_width`) are written in degrees**, matching `to_dict<Scan>` convention. The C++ properties table stores them in radians, so the consolidation loop converts via `s.get_oscillation(deg=True)`. Read side's `extract_properties_table(..., convert_oscillation_to_rad=true)` round-trips them back. Latent bug: every prior test still scan used `oscillation=(0,0)` so radians-vs-degrees was identical until a 10° test case exposed it. Regression tests at `dxtbx/tests/model/test_experiment_list.py:339,370`.
- Each experiment gains `"scan_point": i` and `"scan": 0`.
- Rotation data (oscillation width > 0) is never consolidated.

### Beam re-consolidation — write side (same `to_dict`)

Independent of scan consolidation: fires whenever (a) every experiment is a single-frame still with `frame > 0` and zero oscillation, (b) every scan carries a `"wavelength"` property, (c) every beam is a plain `Beam` sharing geometry (`sample_to_source_direction`, `divergence`, `sigma_divergence`, `polarization_normal`, `polarization_fraction`). Collapses the N monochromatic beams produced by load-side decode or by refinement back to one `XFELBeam` via `make_xfel_beam` with the full polarization tuple. Decoupled from `len(scan_models) > 1`, so it fires on `split_experiments` 1-experiment writes too.

### Scan consolidation — read side

`ExperimentListDict._expand_consolidated_scans()` in `src/dxtbx/model/experiment_list.py`. Called in `__init__` before `_extract_models`. Expands the consolidated object back to N per-frame scan dicts; no-op for files without `__stills_consolidated`. Rest of `decode()` sees the standard format.

**Backward compatibility:** new reader + old file = transparent. Old reader + new file = fails (no `image_range` on the consolidated dict). The dxtbx reader must land before any writer reaches users.

### Decode auto-resolve of `Experiment.beam`

`ExperimentList.decode()` in `src/dxtbx/model/experiment_list.py` (~lines 589–593): for every `Experiment` whose `beam` is an `XFELBeam`, replace `expt.beam` with `expt.get_monochromatic_beam()`. Runs once at load, after `_expand_consolidated_scans`. **Critically:** the imageset's `XFELBeam` slot and the scan's `"wavelength"` property are left untouched — `imageset.get_beam(i)` still dispatches via the Python override, and `to_dict()` re-consolidation puts the single `XFELBeam` back on write.

This is the single point that lets every downstream consumer (`cctbx.xfel.detector_residuals`, `cctbx.xfel.merge`, `dials.show`, `dials.refine`, etc.) handle composite stills `.expt` files without a per-tool fix. The previous per-tool one-liner pattern is now redundant for `.expt`-loaded experiments.

### Sparse imageset round-trip

When `single_file_indices` is stored in the JSON, `ExperimentListDict._make_sequence()` bypasses `ImageSetFactory.make_sequence()` and calls `format_class.get_imageset()` directly — `make_sequence()` recomputes a contiguous range from the scan's `image_range` and asserts `array_range == scan.get_array_range()`, both of which fail for a sparse imageset with compact scan `(1, N)`. Backward-compatible: old files with contiguous `single_file_indices` take the same bypass path and produce an identical result.

### C++ / format guards for stills sequences

- `ImageSequence` C++ constructor in `src/dxtbx/imageset.h`: sequential-indices assertion conditioned on `!scan->is_still()`. Still scans (oscillation width == 0) have independent shots.
- `FormatMultiImage.get_imageset()`: Python-level sequential assertion conditioned on `not (scan is not None and scan.is_still())`. Model-fill guards include `and format_instance is not None` (matches `Format.get_imageset`), needed because `check_format=False` callers (like `combine_experiments`) have `format_instance=None`.
- `Format.get_imageset` goniometer assertion: `assert goniometer is not None or scan.is_still(), ...` — relaxed for still sequences so output `.expt` files (where `expt.goniometer=None` for stills) can be read with `check_format=False`.

---

## dials surface

### `command_line/stills_process.py`

**Helpers (module-level):**
- `_make_stills_sequence()`: builds an `ImageSequence` with zero-oscillation `Scan` and no goniometer. Shared by `do_import()`, its error fallback, and `_rebuild_shared_imageset_output()`.
- `_frame_index(expt)`: recovers 0-based frame index from either a per-frame imageset (`imageset.indices()`) or a full-sequence imageset (`scan.get_array_range()`).
- `extend_with_bookkeeping()`: remaps reflection `id` and `experiment_identifiers` when merging experiment lists. Used by both MPI and multiprocessing paths.
- `_sort_experiments_by_frame(expts, refls)`: sorts experiments by `(imageset.paths()[0], _frame_index(expt))` and remaps reflection-table `id` + `experiment_identifiers` dict to match. Called after `_rebuild_shared_imageset_output` in both output paths.

**`do_import()` pivot:** for multi-image formats, calls `format_class.get_imageset(as_sequence=True)` to get a shared-identity `ImageSequence`. If that raises `RuntimeError` (format's `get_scan()` returns None), wraps the `ImageSet` in a synthetic `ImageSequence` via `_make_stills_sequence()`. Per-frame slices are `imageset[i:i+1]`, then rebuilt with a fresh `Scan((1,1))` to normalize `array_range` to `(0,1)`. **Every `ImageSequence` slice is rebuilt regardless of whether the original scan was a still or a rotation** — rotation CBF frames processed via `convert_sequences_to_stills` carry absolute scan ranges (e.g., frame 2 → `array_range=(1,2)`) that would confuse the integrator just as much as XFEL absolute ranges. `load_models()` re-reads the synthetic scan back onto the experiment, so scan/goniometer must be re-nulled after every `load_models()` call (in both `do_import` and `do_work`).

**Output rebuilding — `_rebuild_shared_imageset_output()`:** groups experiments by source file; re-opens source to recover full `ImageSetData`; builds the shared `ImageSequence` with **only the integrated frame indices** (`sorted(set(frame_indices))`) and a compact scan `(1, N)`. Creates one `Scan((fi+1, fi+1))` per unique indexed frame (multi-lattice frames share). Guards against rotation data. This is what shrinks the file *and* makes `image_viewer` show only successfully integrated frames.

**Output rebuilding — `_combine_multiprocessing_outputs()`:** merges per-worker composites, calls `_rebuild_shared_imageset_output`, then `_sort_experiments_by_frame`. Writes via `*.tmp` + `os.replace` (atomic) so a crash mid-write cannot corrupt worker-0's filename; logs cleanup `OSError` on workers 1..N-1 rather than swallowing.

**`combine_all_ranks` phil:** triggers merge at end (requires `composite_output=True`). MPI path sets `composite_stride = comm.size` so rank 0 gets all frames.

**Indexer caching on `Processor`:** `self.spot_finder_factory` built once via `get_spot_finder_factory()` and reused per frame. `self.idxr`, `self.idxr_method_list`, `self.idxr_known_crystal_models` cached; `update_indexer()` resets per-frame state instead of reconstructing. The `copy.deepcopy(all_params)` previously around the indexer call is replaced with direct mutate-and-restore of `refinement.parameterisation.scan_varying` (forced False) and `basis_vector_combinations.max_refine` (pinned to 5, then restored). **`max_refine=5` is a behavioral change**, not just a refactor — it was previously derived from params.

### `algorithms/spot_finding/finder.py`

`ExtractPixelsFromImage`: added `is_stills` param; `self.mask` → `self.image_mask` (static component only); `__call__()` merges static mask with dynamic `imageset.get_mask(index)` per call. New `update_imageset()` swaps in next frame's data without rebuilding — static mask survives, dynamic mask is re-read. `ExtractSpots` and `SpotFinder._find_spots()` cache `ExtractPixelsFromImage` in `self.function`; call `update_imageset()` on subsequent frames. **Spotfinder XFEL detection uses `scan.has_property("wavelength")`**, not `isinstance(beam, XFELBeam)` — because per-frame imagesets carry monochromatic beams built by `do_import`, not XFELBeam.

### `algorithms/spot_finding/factory.py`

`SpotFinderFactory.from_parameters()` accepts explicit `is_stills`. When True, sets `no_shoeboxes_2d = True` directly rather than checking `isinstance(imageset, ImageSequence)`.

### `algorithms/indexing/indexer.py`

Stills detection via `Experiment.is_still()`. Deleted: the block that used to convert `ImageSequence → ImageSet` and null scan/goniometer for stills.

**P0 correctness:** `Indexer.refined_experiments` is nulled at the **start** of `index()`, not only in `__init__`. A frame whose indexing fails before the success block at the end of `index()` no longer silently inherits the previous frame's refined result on the cached `Indexer`. Warm-start state (`d_min`, detector/crystal geometry) is intentionally left intact — warm-start is for *inputs*, not *outputs*.

### `algorithms/indexing/stills_indexer.py`

- `self.hardcoded_phil` moved to `__init__()` (was re-parsed per call to `choose_best_orientation_matrix()`).
- `copy.deepcopy(self.all_params)` removed from `choose_best_orientation_matrix()` — params are not mutated here; mutation+restore happens in `Processor.index()`.
- `experiment_list_for_crystal()`: loop iterates experiments (not imagesets), so beam comes from `expt.beam` (wavelength-adjusted for XFEL after decode auto-resolve) not `imageset.get_beam()` (base shared beam).

### `algorithms/integration/integrator.py`

`_determine_max_memory_needed()`: prefers experiment's scan over imageset's scan for stills detection.

### `algorithms/spot_prediction/reflection_predictor.py`

Stills detection: `experiment.scan is not None and not experiment.scan.is_still()`. Removed `from dxtbx.imageset import ImageSequence` import.

### `algorithms/refinement/refiner.py`

`_copy_experiments_for_refining()` resolves per-frame `XFELBeam` via `Experiment.get_monochromatic_beam()` — `new_exp.beam = new_exp.get_monochromatic_beam()`. Single chokepoint for refinement: every internal consumer (`StillsReflectionManager`, `StillsReflectionPredictor`, `BeamParameterisation`, prediction parameterisations, residual helpers) sees an ordinary monochromatic `Beam`. Refined output flows back through the same helper and is re-consolidated to one `XFELBeam` by `to_dict()`.

Redundant for `.expt`-loaded experiments after the `decode()` auto-resolve landed, but still needed for runtime-built experiments and for keeping the on-disk `XFELBeam` intact on a true no-op `combine_experiments reference_from_experiment.beam=0`.

### `algorithms/refinement/parameterisation/configure.py`

`_parameterise_crystals()`: a zero-oscillation still scan no longer counts as a rotation experiment. Guard counts a scan as rotation only when `not scan.is_still()`.

### `util/combine_experiments.py` and `command_line/combine_experiments.py`

`CombineWithReference.__call__`:
- When both `ref_scan.is_still()` and `experiment.scan.is_still()`, the experiment keeps its own scan (`reference_from_experiment.scan=0` is effectively a no-op for stills — N per-frame scans get auto-consolidated on write).
- Relaxes `wavelength_tolerance` to `float("inf")` when both ref and current are stills. For stills the per-frame wavelength lives on the scan property, not the beam — comparing beam wavelengths is meaningless and the strict 1e-6 default fails on any multi-frame XFEL input. Direction, polarization-normal, polarization-fraction tolerances are unchanged.

`command_line/combine_experiments._consolidate_stills_imagesets`: builds the shared `ImageSequence` with the **sparse layout** (`sorted(set(frame_indices))` + compact `Scan((1, N))`) matching `stills_process._rebuild_shared_imageset_output`. `_sort_experiments_and_reflections` triggers automatically when `expts[0].is_still()`. Stills scan consolidation on write is automatic.

### `util/image_viewer/spotfinder_frame.py` + `slip_viewer/frame.py`

- `_identifiers_for_frame`: order-independent `{abs_frame: [identifiers]}` dict built by scanning all experiments. The previous binary-search implementation assumed ascending frame order; `_combine_multiprocessing_outputs` interleaves worker composites so the assumption was false — every frame except the first returned empty identifiers.
- `_still_scan_frame_list(imageset)`: sorted list of `expt.scan.get_array_range()[0]` for experiments sharing the imageset (cached per imageset by `id()`). `_still_scan_frame_offset` is a one-line wrapper returning `frame_list[0]`.
- `get_spotfinder_data` frame mapping: `i_frame = frame_list[idx]` when sparse, fallback `idx + offset` for old contiguous `.expt`.
- `viewing_still_scans`: integrated shoeboxes use fixed purple (`#984ea3`), predictions fixed green (`#4daf4a`). `prediction_colours` has 90 entries; composite-table id can reach 764+ on real cases.
- `slip_viewer/frame.py`: `chooser_wrapper.get_beam()` → `image_set.get_beam(self.index)`. Same pattern for `get_beam_center_px` and resolution rings.

### `command_line/expt_set_to_sequence.py` / `expt_sequence_to_set.py`

Round-trip conversion between old per-frame `ImageSet .expt` and new `ImageSequence` format. Both import JSON helpers from `dials.util.stills_imageset_convert`: `source_file`, `frame_index`, `build_consolidated_scan`, `expand_consolidated_scan`, `per_frame_scan_map` (eliminated a third copy of the consolidate/expand logic). Set→sequence sorts experiments by `(src, fi)` ascending and remaps `.refl` `id` + `experiment_identifiers`. Sequence→set preserves experiment order so `.refl` is copied through unchanged.

Default outputs: `sequence.expt`/`sequence.refl` and `imageset.expt`/`imageset.refl`.

---

## Sharp edges / gotchas

1. **C++ `ImageSequence` requires a non-null scan.** Two workarounds — don't conflate:
   - XFEL formats: `FormatXFEL.get_imageset()` builds zero-oscillation `Scan` with `"wavelength"` property. No error to catch.
   - Non-XFEL multi-image formats: `get_imageset(as_sequence=True)` raises `RuntimeError`, `do_import` catches it and wraps via `_make_stills_sequence()`.

2. **Per-frame scan range is absolute** for any `ImageSequence` slice — XFEL frame 3913 has `array_range=(3912,3913)`, rotation CBF frame 2 has `array_range=(1,2)`. The integrator expects 0-based z. Every `ImageSequence` slice is rebuilt with `Scan((1,1))` in `do_import` regardless of original scan type.

3. **`Format.get_imageset` goniometer assertion is relaxed for still sequences.** Required so composite `.expt` files (where `expt.goniometer=None`) can be read back with `check_format=False`.

4. **`load_models()` re-reads the synthetic still scan** back onto the experiment. Scan/goniometer must be re-nulled after every `load_models()` call.

5. **`ImageSequence.get_beam()` with no index returns `XFELBeam`** for XFEL data — callers needing a monochromatic beam must call `get_beam(i)`.

6. **Spotfinder XFEL detection uses `scan.has_property("wavelength")`**, not `isinstance(beam, XFELBeam)` — because per-frame imagesets carry monochromatic beams built by `do_import`.

7. **`easy_mp` tried to pickle C++ threshold/indexer objects** via the result queue. `do_work` returns `None` when `finalize=True` since output is already on disk.

8. **`max_refine=5` pin** in the indexer cache path is a behavioral change for stills.

9. **Consolidated scan oscillation is written in degrees**, not the raw radians from the C++ properties table. The consolidation loop must use `s.get_oscillation(deg=True)[0|1]`, not `get_properties()` directly. Tests at `test_experiment_list.py:339,370` catch the radians regression.

10. **Decode-side auto-resolve does not touch the imageset.** `Experiment.beam` becomes a monochromatic `Beam`; `imageset.get_beam(0)` still returns the `XFELBeam`. Don't conflate the two when reasoning about round-trip.

11. **Beam re-consolidation guard requires shared geometry.** Refining beam parameters per frame will (correctly) fail the guard and write N beams — see pending work.

---

## Pending work

| # | Item | Status |
|---|------|--------|
| 1 | `XFELBeam` permanence as a dxtbx model | Open — social/group decision; code is in place at `beam.h:887` with full constructor + `get_num_scan_points=0` |
| 3 | Dedicated tests for `XFELBeam` / `get_monochromatic_beam` / `ImageSequence.get_beam(i)` | **Verified missing** — zero hits across `dxtbx/tests/` for these symbols. Consolidation round-trip tests at `test_experiment_list.py:339,370` exercise the class indirectly but don't substitute for guard-throws / factory-passthrough / dict-roundtrip / polarization-preservation unit tests |
| 5 | `combine_all_ranks` output semantics | Open — param at `stills_process.py:125`; behavior is "merge at end with `composite_stride = comm.size`"; design intent not formally settled |
| 6 | Performance benchmarking (file size + runtime) | Not done — no benchmark scripts in either repo. File-size: −267 KB / 8.1 % on 764-experiment test case is the only number recorded. No spotfinder/indexer caching speedup numbers; no end-to-end runtime comparison |
| — | Beam refinement for XFEL stills | Open by design — `to_dict()` beam re-consolidation guard requires shared geometry; refining beam params per frame fails the guard (correctly) and writes N differing beams. `refine_level0.phil` fixes the beam (`beam.fix=*all`) and avoids this. Needs separate design |
| — | `image_viewer` overlay precompute | Not implemented — predictions path is still O(N_expt × N_refl) per frame flip (~3.7×10⁸ on 764-expt × 484,603-refl test case). A precomputed `{abs_frame: flex.size_t(row_indices)}` dict built at load would drop both `bbox/ctr-mass/shoebox` and `predictions` paths to O(rows-on-frame) per flip (~630 rows on the test case) |
| — | `CombineWithReference` does not propagate `experiment.profile` | Pre-existing on `main`; orthogonal to this branch; surfaced by split→combine round-trip verification but not fixed here |
| — | `cctbx.xfel.detector_residuals.py:453` redundant-but-kept | The `expt.beam = expt.get_monochromatic_beam()` loop is now redundant for `.expt`-loaded experiments after the `decode()` auto-resolve. Kept because runtime-built (non-decode) experiments still need it. Add a one-line comment noting intent so a future cleanup pass doesn't remove it |

**Resolved (per prior log):** Q2 (route via `Experiment.is_still()` — done 2026-05-24); Q4 (downstream consumers — mitigated 2026-05-24 by decode auto-resolve, modulo a still-outstanding verification against `cctbx.xfel.merge` and the XFEL GUI); Q7 (consolidation oscillation test — done 2026-05-26).

---

## Quick file reference

### dxtbx

| File | Role |
|------|------|
| `src/dxtbx/model/beam.h` | `XFELBeam` C++ class (guard, no wavelength array); `get_num_scan_points` returns 0 |
| `src/dxtbx/model/boost_python/beam.cc` | Python bindings + pickle + `to_dict`/`from_dict`; `from_dict` `polarization_fraction` default 0.5 |
| `src/dxtbx/model/beam.py` | `BeamFactory.make_xfel_beam(direction, divergence, sigma_divergence, polarization_normal, polarization_fraction, flux, transmission, probe, sample_to_source_distance)`; routing in `from_dict()` |
| `src/dxtbx/model/__init__.py` | `XFELBeam.get_monochromatic_beam(wavelength)` injection; `Experiment.get_monochromatic_beam()` injection; `ExperimentList.to_dict/as_json` scan consolidation + beam re-consolidation |
| `src/dxtbx/model/experiment_list.py` | `_expand_consolidated_scans()` on load; `decode()` auto-resolves `Experiment.beam` from `XFELBeam` → per-frame `Beam`; `_make_sequence()` bypasses `make_sequence()` when `single_file_indices` is in JSON |
| `src/dxtbx/imageset.h` | `ImageSequence` ctor — sequential-indices assertion conditioned on `!scan->is_still()` |
| `src/dxtbx/imageset.py` | `ImageSequence.get_beam(index)` override (XFELBeam dispatch) |
| `src/dxtbx/format/FormatXFEL.py` | Mixin — builds `ImageSequence` with `XFELBeam` + scan wavelength property; subset fix; full polarization passthrough |
| `src/dxtbx/format/FormatMultiImage.py` | Model-fill guards (`format_instance is not None`); sequential-indices assertion conditioned on `not scan.is_still()` |
| `src/dxtbx/format/FormatNXmx.py` | `FormatNXmxXFEL` class; `understand()` guarded with try/except |
| `src/dxtbx/format/FormatXTC.py` | `FormatXTCXFEL` class; `FormatXTC.understand()` returns `bool(ds)` |
| `src/dxtbx/format/Format.py` | Goniometer guard restored to `assert ... or scan.is_still()` idiom |
| `src/dxtbx/nexus/__init__.py` | `CachedWavelengthBeamFactory.get_wavelengths()` returns `list[float]` |
| `dxtbx/tests/model/test_experiment_list.py` | `test_stills_consolidated_scan_oscillation_roundtrip`, `test_stills_consolidated_rotation_not_consolidated` |

### dials

| File | Role |
|------|------|
| `src/dials/command_line/stills_process.py` | `do_import`, helpers, indexer/spotfinder caching, output rebuild, atomic `*.tmp` + `os.replace` combine, `combine_all_ranks` |
| `src/dials/util/stills_imageset_convert.py` | **New.** Shared JSON helpers: `source_file`, `frame_index`, `build_consolidated_scan`, `expand_consolidated_scan`, `per_frame_scan_map` |
| `src/dials/algorithms/spot_finding/finder.py` | Static/dynamic mask split; `update_imageset`; `SpotFinder` cache |
| `src/dials/algorithms/spot_finding/factory.py` | `is_stills` param (fixed contract miss via `Experiment.is_still()` routing) |
| `src/dials/algorithms/indexing/indexer.py` | Scan-based stills detection; deleted `ImageSet` conversion; `refined_experiments` reset at `index()` entry |
| `src/dials/algorithms/indexing/stills_indexer.py` | `hardcoded_phil` moved to `__init__`; deepcopy removed; experiment-loop fix |
| `src/dials/algorithms/integration/integrator.py` | Memory calc stills detection |
| `src/dials/algorithms/spot_prediction/reflection_predictor.py` | Stills dispatch |
| `src/dials/algorithms/refinement/refiner.py` | `_copy_experiments_for_refining` resolves XFELBeam via `get_monochromatic_beam()` |
| `src/dials/algorithms/refinement/parameterisation/configure.py` | `_parameterise_crystals` — still scan is not a rotation experiment |
| `src/dials/util/combine_experiments.py` | `CombineWithReference.__call__` — stills keep own scan; wavelength tolerance relaxed for stills |
| `src/dials/command_line/combine_experiments.py` | `_consolidate_stills_imagesets` sparse layout; auto-sort for stills; `Experiment.is_still()` routing |
| `src/dials/command_line/expt_set_to_sequence.py`, `expt_sequence_to_set.py` | Round-trip converters; use shared `stills_imageset_convert` helpers |
| `src/dials/util/image_viewer/spotfinder_frame.py` | `_identifiers_for_frame` (order-independent dict); `_still_scan_frame_list` (sparse mapping); fixed per-overlay colours |
| `src/dials/util/image_viewer/slip_viewer/frame.py` | Per-frame beam via `image_set.get_beam(self.index)` |

### cctbx_project

| File | Role |
|------|------|
| `xfel/command_line/detector_residuals.py` | `expt.beam = expt.get_monochromatic_beam()` loop after `flatten_experiments` — now redundant for `.expt`-loaded experiments after dxtbx `decode()` auto-resolve, still needed for runtime-built experiments |

---

## Verified end-to-end

- 764-experiment composite stills `.expt` reload via `ExperimentListFactory` → every `expt.beam` resolves to a monochromatic `Beam` with the per-frame wavelength.
- `dials.refine refine_level0.phil` on the same file → completes; `refined_level0.expt` has exactly one `"__id__": "xfel"` beam; consolidated scan keeps the 726-value per-frame `wavelength` array.
- `cctbx.xfel.detector_residuals hierarchy=1` on the same file → completes; prints full per-panel-group residual statistics.
- `dials.split_experiments combined.expt combined.refl` → `dials.combine_experiments split_*.expt split_*.refl reference_from_experiment.beam=0 detector=0 scan=0` → output matches original on beam (one `__id__: xfel`, `pf=0.999`), scan (one `__stills_consolidated` with bit-identical 764 `frame_numbers` and `wavelength` arrays), detector, imageset (same template + `single_file_indices`). Crystal differs only at ~1e-14 round-off.
- `dials.image_viewer` on composite output → reflection overlays (spots, COM, predictions, hkl, integrated shoeboxes) draw on every integrated frame; chooser shows only successfully integrated frames (sparse-imageset path).
- `dials.combine_experiments` over per-rank composites → produces a single composite matching the `stills_process` direct-output structure (sparse imageset layout, consolidated scan, single XFELBeam).

---

## Chronological appendix

Compact timeline. Use `git log` per commit hash for full diff.

- **2026-05-21** Architecture refactored per upstream maintainer review (Waterman, McDonagh). `XFELImageSequence` removed; per-frame wavelengths moved to scan properties table. `XFELBeam` retained as guard/marker. `ExperimentList.to_dict/as_json` introduces automatic stills scan consolidation (−267 KB on 764-experiment test). `compact_stills_scans` flag introduced.
- **2026-05-21** Per-frame scan rebuild in `do_import` extended to all `ImageSequence` slices (incl. rotation CBF). `FormatXTC.understand()` returns `bool(ds)`. `Format.get_imageset` goniometer assert relaxed for still sequences. `test_pseudo_scan` strengthened.
- **2026-05-21** `tests/algorithms/profile_model/ellipsoid/conftest.py` test fixture rebuilt to construct a real `ImageSequence` (commit `ae871b468`) — fixture-alignment policy: don't weaken production code for malformed test data.
- **2026-05-22** `dials.image_viewer` reflection overlays fixed for composite stills (commit `875d1a8cc`). Two root causes: `_identifiers_for_frame` binary-search on unsorted list; imageset scan loses start offset on JSON round-trip. Per-overlay colours fixed for `viewing_still_scans`.
- **2026-05-22** `dials.combine_experiments` made to work with stills-as-ImageSequence output (dxtbx `e3133a6e`, dials `9494e3a0b`). Four fixes: `FormatMultiImage.get_imageset` `format_instance is not None` guards; `FormatXFEL.get_imageset` wavelength subset; `ExperimentListDict.decode` scan-override skip; `CombineWithReference.__call__` keeps per-experiment scan for stills.
- **2026-05-22** Frame-order canonicalization in composite output: `_sort_experiments_by_frame` in `stills_process` and `combine_experiments`.
- **2026-05-22** Scan consolidation made unconditional — `compact_stills_scans` parameter removed (dxtbx `e2d34c52`; dials `dd8bfea5c`).
- **2026-05-22** Sparse imageset decode crash fix (dxtbx `1467eafe`) — `eobj_scan` merger no longer overwrites compact scan with bad merged range.
- **2026-05-22** Sparse imageset / image_viewer frame filter (dxtbx `eea48ad7`, dials `dd8bfea5c`, `a7b272a4f`). `_rebuild_shared_imageset_output` switches to `sorted(set(frame_indices))` + compact `Scan((1, N))`. Three code paths updated for non-contiguous `single_file_indices`. `_still_scan_frame_list` for sparse frame mapping in image_viewer.
- **2026-05-22** `dials.refine` works on composite stills (dxtbx `cd5d39e8`, dials `272a1b469`, cctbx_project `a34f7ddc`). `Experiment.get_monochromatic_beam()`; refinement chokepoint at `_copy_experiments_for_refining`; `to_dict()` re-consolidates per-frame beams. `_parameterise_crystals` no longer treats still scan as rotation.
- **2026-05-24** Review hardening pass (dxtbx `64b7c062`, `16294afa`, `266edaf9`, `2498250f`; dials `054bb1fb1`, `6ee2a8f8e`, `3b4083a9f`, `c02872d72`, `5e4d11e15`, `61dfdab69`). P0: `Indexer.refined_experiments` reset at `index()` entry; `make_xfel_beam` polarization passthrough end-to-end; `ExperimentList.decode()` auto-resolves `XFELBeam` → `Beam`. P1: `Experiment.is_still()` routed everywhere; sparse imageset layout in `combine_experiments`; shared `stills_imageset_convert` helpers; image_viewer frame-helper collapsed; atomic `*.tmp` + `os.replace` in `_combine_multiprocessing_outputs`; `XFELBeam.get_num_scan_points=0`; `from_dict` `pf` default 0.5; `FormatNXmxXFEL.understand` try/except; `Format.py` assert idiom; pinned `ruff format` pass.
- **2026-05-24** `split_experiments` → `combine_experiments` round-trip fixed. `to_dict()` beam re-consolidation decoupled from scan consolidation (fires on 1-experiment writes); full polarization passthrough to `make_xfel_beam`. `CombineWithReference.__call__` relaxes wavelength tolerance to infinity for stills.
- **2026-05-26** Consolidated scan oscillation unit fix — write side now uses `get_oscillation(deg=True)` instead of raw radians from `get_properties()`. Regression tests added at `test_experiment_list.py:339,370`. Closes Q7.
- **2026-05-26** `tests/command_line/test_combine_experiments.py::test_min_max_reflections_per_experiment` test fixture rebuilt to attach a real `ImageSequence` (commit `963d588ea`) — fixture-alignment policy applied again.

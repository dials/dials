# `stills_process_imagesequence`: architectural narrative

Branches: `stills_process_imagesequence` in **dials**, **dxtbx**, and **cctbx_project**.
A prototype, not intended for merge. The purpose of this document is to give senior
DIALS / dxtbx developers enough of the design to reason about it without reading the
code, and to anchor the group conversation that produces the production version.

## 1. Motivation

`dials.stills_process` has historically converted every imported image into a fresh
per-frame `ImageSet` inside `do_import()`. That conversion destroyed *object identity*:
imageset, beam, and detector were rebuilt as new Python objects for every frame.

Two consequences followed, both of which this branch removes:

- The JSON serializer deduplicates models by identity; with per-frame copies it could not
  deduplicate. Experiments from one source file each carried their own imageset, beam,
  and detector — N copies for N frames.
- Per-frame caches in the spotfinder and indexer key on imageset identity. With
  per-frame copies they could not be implemented cleanly.

The branch keeps the data as an `ImageSequence` for the entire pipeline. Everything
else in this document is either a direct consequence of removing the conversion, or a
consequence of the new contracts that removal forces on the rest of DIALS.

## 2. The pivot: `do_import()` keeps data as `ImageSequence`

`do_import()` calls `format_class.get_imageset(as_sequence=True)` instead of
`ImageSetFactory.imageset_from_anyset()`. Per-frame slices are `imageset[i:i+1]` and
are then rebuilt with a fresh `Scan((1, 1))` to normalize the absolute frame range to
`(0, 1)`. The slice is itself an `ImageSequence`, sharing the parent imageset's beam
and detector.

This is the root architectural change of the branch: object identity (the same Python
object visible across frames) is what JSON dedup and the per-frame caches both need.
No other piece of the work functions without it.

Every slice is rebuilt regardless of whether the original scan was a still or a
rotation. Rotation CBF frames processed via `convert_sequences_to_stills` also carry
absolute scan ranges (e.g., frame 2 → `array_range=(1, 2)`) that confuse the
integrator's 0-based z assumption just as much as XFEL absolute ranges. The
normalization is upstream of XFEL-specific concerns.

## 3. New dxtbx primitives for XFEL data

Three new primitives in dxtbx; each is small, and each carries an explicit decision.

### 3a. `XFELBeam` — lightweight guard / marker

C++ class inheriting `Beam` (in `model/beam.h` and `model/boost_python/beam.cc`).
Stateless with respect to wavelength: carries direction, divergence, polarization,
flux, transmission, probe, sample-to-source distance. `get_wavelength()`,
`set_wavelength()`, `get_s0()`, `set_s0()` and all scan-point methods **throw**.
Serializes with `"__id__": "xfel"`.

- **Why not `PolychromaticBeam`.** `PolychromaticBeam` models a source with an
  intrinsic spectral bandwidth (Laue / ToF) and stores a `wavelength_range`. XFEL is a
  sequence of monochromatic pulses with varying energies. The per-shot energies are
  *measurements*, not a property of the source. Conflating the two is incorrect physics.
- **Why a throwing guard rather than a default wavelength.** Code that calls
  `beam.get_s0()` on an XFEL beam has a real semantic gap — it has not chosen a frame.
  A silent fallback (mean, first, nominal) would mask bugs at exactly the boundaries
  where they need to surface.
- **Why `get_num_scan_points()` returns 0 instead of throwing.** Generic `Experiment`
  copy / compare / serialize code probes scan-varying state as a benign check.
  Returning 0 matches the real semantics — an `XFELBeam` genuinely has no scan points —
  and avoids forcing XFEL-awareness into generic code.

### 3b. Per-frame wavelengths in the scan properties table

*Design due to David Waterman and David McDonagh; implementation by McDonagh.*

Wavelengths are stored as `scan.set_property("wavelength", flex.double([...]))`. The
property is written by `FormatXFEL.get_imageset()` for XFEL data, and read by
`ImageSequence.get_beam(i)` (§3c), by `do_import` (which stamps a single-entry array
onto each per-frame scan so the spotfinder can detect XFEL data on the in-flight path),
and by `_rebuild_shared_imageset_output` (which stamps the full per-frame array onto
the output scan so the property survives JSON round-trip).

The design principle behind this is that **the beam should not need to carry
per-image information**. The beam represents the source; per-frame wavelength is a
per-image measurement and belongs on the scan alongside other per-image data such as
timestamps and exposure times. The same model handles the constant-wavelength case
without modification — the `"wavelength"` property is simply absent — and the
property gets consolidated serialization for free (§5c).

### 3c. `ImageSequence.get_beam(index)` override

Per-frame beam access is handled by a Python method override injected onto the
existing `ImageSequence` class (in `imageset.py`). When the imageset's beam is an
`XFELBeam` *and* the scan has a `"wavelength"` property, `get_beam(i)` reads
`scan.get_property("wavelength")[i]` and returns a regular monochromatic `Beam` built
from the `XFELBeam`'s direction / divergence / polarization tuple and the per-frame
wavelength. All other paths (no index, non-`XFELBeam`, non-XFEL imagesets) fall through
unchanged.

- **Why a method override rather than an `XFELImageSequence` subclass.** An earlier
  iteration of this work introduced an `XFELImageSequence` subclass. The dxtbx
  maintainers (Waterman, McDonagh) pushed back: a new ImageSet class proliferates type
  checks across dxtbx and DIALS. The override is invisible to non-XFEL callers, and
  the XFEL-specific dispatch is encapsulated inside dxtbx — it never leaks into DIALS.
- **No public per-frame beam API beyond `get_beam(i)`.** `imageset.get_beam(i)` is the
  existing contract. Code that already called it works without modification.

### 3d. `FormatXFEL` mixin and the two concrete formats

The mixin's `get_imageset()` lets the parent format build the raw `ImageSet`, reads
per-frame wavelengths via a subclass-provided `get_wavelengths()` returning
`list[float]` (Ångström), constructs the `XFELBeam` from the first frame's geometry
(passing the source beam's full polarization tuple through, see §7), and returns a
regular `ImageSequence` with that `XFELBeam` and a zero-oscillation `Scan` carrying
the `"wavelength"` property.

`FormatNXmxXFEL` delegates wavelength reading to
`CachedWavelengthBeamFactory.get_wavelengths()` (the 1-D `incident_wavelength` HDF5
dataset). Its `understand()` requires `wl.ndim > 0 and wl.size > 1`: scalar wavelengths
and 1-element arrays represent a single fixed wavelength for all frames; those files
belong to `FormatNXmx`, not the XFEL specialization. `FormatXTCXFEL` delegates
`understand()` to `FormatXTC` because XTC streams are always XFEL stills with
per-event wavelengths.

## 4. Cross-cutting contract: stills detection by scan, not imageset type

The old detection idiom — `isinstance(imageset, ImageSequence)` is False ⇒ still — no
longer separates stills from rotations, because both are now `ImageSequence`. The new
contract is:

> A still is `expt.scan is None or expt.scan.is_still()` (zero oscillation width).

Every check has been routed through `Experiment.is_still()`, the C++ accessor that
already encodes the full scan-is-still / scan-is-None / goniometer-is-None contract.
The sites touched include the indexer, the integrator's memory estimator, the
reflection predictor's static-vs-scan-varying dispatch, the refinement crystal
parameterisation, the spotfinder factory, and `combine_experiments`.

This is the widest blast radius in the branch. Centralizing through `is_still()` is
load-bearing: the 2026-05-24 review pass found that `spot_finding/factory.py`'s
`is_stills=False` fallback was still using the old isinstance check, silently
inverting `no_shoeboxes_2d` for stills entering via non-`stills_process` callers. The
question of whether scan-based detection is the right pipeline-wide contract is a
decision the group should ratify; this branch's position is that there should be a
single check, and it should be `Experiment.is_still()`.

A separate, narrower question is *XFEL* detection in the spotfinder. That uses
`scan.has_property("wavelength")` rather than `isinstance(beam, XFELBeam)`, because
the in-flight path carries monochromatic beams built by `do_import` — an
`XFELBeam` check would always be False there.

## 5. Output format: shared models, sparse imagesets, consolidated scans

Three layered transformations on the way out, each addressing a distinct cause of
file bloat.

### 5a. Shared models per source file

`_rebuild_shared_imageset_output()` groups experiments by source file and rebuilds them
so every experiment from one file points at one `ImageSequence` / beam / detector
object. Multi-lattice experiments from the same frame share one `Scan`. A no-op for
rotation data. This is the change that actually shrinks the file.

### 5b. Sparse imageset layout

The shared `ImageSequence` is built with indices `sorted(set(frame_indices))` — only
the frames that were successfully integrated — and a compact `Scan((1, N))`. The
alternative (a contiguous `range(min_f, max_f + 1)`) would show every frame in the
source-file range in `image_viewer`, including frames that never produced a result.

Allowing non-contiguous indices required three downstream changes in dxtbx: the C++
`ImageSequence` constructor's sequential-indices assertion is now conditioned on
`!scan->is_still()`; the Python-level assertion in `FormatMultiImage.get_imageset()`
is similarly conditioned; and `ExperimentListDict._make_sequence()` bypasses
`ImageSetFactory.make_sequence()` when `single_file_indices` is present in the JSON,
because `make_sequence()` recomputes a contiguous index range from the scan's
`image_range` and asserts equality with `array_range`, both of which fail for a sparse
imageset with compact scan.

### 5c. Consolidated stills scan, re-consolidated `XFELBeam`

When every scan in an `ExperimentList` is a single-frame still (oscillation width 0),
`ExperimentList.to_dict()` collapses the N scan dicts into one object:

```json
{"__stills_consolidated": true, "frame_numbers": [...],
 "properties": {"wavelength": [...]}, ...}
```

All-zero property arrays are omitted and reconstructed on load. Each experiment gains
`"scan_point": i`; all point to `"scan": 0`. On the read side,
`ExperimentListDict._expand_consolidated_scans()` runs in `__init__` before any model
extraction, so the rest of `decode()` sees the standard format. Measured reduction on
a 764-experiment file: about 8 % of total size.

The same `to_dict()` block re-consolidates per-frame monochromatic beams back into a
single `XFELBeam` when every experiment beam is a plain `Beam` sharing geometry
(direction, divergence, polarization). This preserves the canonical single-beam
representation through `refine` / `split` / `combine` round-trips. If beams have
diverged (e.g., genuinely refined to different geometries), all N beams are written —
no data is lost.

Consolidation runs automatically. An earlier opt-in flag (`compact_stills_scans=True`)
was removed: requiring 200+ downstream call sites to opt in is fragile, and the guard
(all scans single-frame, oscillation 0) makes the transformation a true no-op for
rotation data. Backward compatibility: a new reader handles old files transparently;
old readers fail on new files.

## 6. Performance caches in `Processor`

The spotfinder and indexer caches were developed and validated on a separate branch
(`stills_process_profiling`) and do not depend on `ImageSequence`. Nick Devenish
proposed exploring the `ImageSequence` conversion specifically so that those caches
*could be simplified* — stable imageset identity reduces the bookkeeping the cache
implementation needs. This branch carries the caches forward in the simpler form.

**Spotfinder.** `ExtractPixelsFromImage` splits its mask into a static component
(lookup mask + untrusted rectangles, frame-independent) and a dynamic component
(`imageset.get_mask(index)`). A new `update_imageset()` swaps the per-frame data
without rebuilding. `ExtractSpots` and `SpotFinder` cache the
`ExtractPixelsFromImage` instance and call `update_imageset()` on subsequent frames
instead of constructing a new one. `Processor` also caches the configured
`SpotFinderFactory` itself.

**Indexer.** `Processor` caches the `Indexer` instance. A new `update_indexer()`
resets per-frame state instead of reconstructing. `copy.deepcopy(all_params)` per
image has been replaced with direct mutate-and-restore of
`refinement.parameterisation.scan_varying` (forced False for stills). The
`max_refine` value is no longer touched at the call site: `Indexer.from_parameters`
already resolves the `libtbx.Auto` sentinel to 5 for stills on the first frame, and
the cached `Indexer` inherits the resolved value on subsequent frames — so behavior
matches `main`, where each frame's `copy.deepcopy` reached the same end-state through
the same resolver. An explicit user-supplied `max_refine` value is preserved.

A correctness item that earlier versions of this branch carried — `Indexer.refined_experiments`
being nulled only in `__init__`, allowing a frame whose indexing failed before the
success block of `index()` to silently inherit the previous frame's refined result
through the cached `Indexer` — has been fixed: `index()` now resets
`self.refined_experiments = None` at entry, so the
`if self.refined_experiments is None: raise DialsIndexRefineError(...)` guard fires
correctly on per-frame failures.

## 7. Downstream contract — what consumers must do

The XFEL portion of the contract — `Experiment.beam` may be an `XFELBeam` that throws
on `get_s0()`, with per-frame wavelengths living on the scan — has a non-trivial blast
radius. The branch addresses it with one central mitigation and a small set of
targeted patches in DIALS tools that bypass the central path.

**Central mitigation.** `ExperimentList.decode()` auto-resolves each `Experiment.beam`
from `XFELBeam` to a per-frame monochromatic `Beam` after `_expand_consolidated_scans`
has restored the per-frame scans. Any consumer that loads a composite `.expt` via
`ExperimentListFactory.from_*` sees a regular `Beam` on `Experiment.beam` and requires
no modification. Crucially, the imageset's `XFELBeam` and the scan's `"wavelength"`
property are left untouched, so `imageset.get_beam(i)` continues to dispatch through
the per-frame override, and the write-side re-consolidation in `to_dict()` is
unaffected.

**Targeted patches** were required where the auto-resolve does not apply — either
because the tool copies experiments instead of reloading them, or because the original
behavior was wrong in a way the auto-resolve does not catch:

- `dials.refine` builds an experiment *copy* for refinement
  (`_copy_experiments_for_refining`). The copy is now resolved via
  `Experiment.get_monochromatic_beam()` at that single chokepoint, so every
  refinement-internal consumer (`StillsReflectionManager`, `StillsReflectionPredictor`,
  `BeamParameterisation`, the prediction parameterisations) sees an ordinary
  monochromatic beam with no code change. Refined output flows back through the same
  helper and is re-consolidated to one `XFELBeam` on write. A separate fix in
  `_parameterise_crystals` stops classifying a zero-oscillation still scan as a
  rotation experiment.
- `dials.combine_experiments` needed four fixes spanning both repos: a
  `FormatMultiImage` guard against `format_instance=None` with `check_format=False`;
  a `FormatXFEL` wavelength-subset fix when loading composite output; an
  `ExperimentListDict.decode` change to not overwrite the format-provided scan with
  the JSON-decoded one; and a `CombineWithReference` change that keeps per-frame
  scans on stills and relaxes `wavelength_tolerance` to infinity for stills↔stills
  comparisons (after the auto-resolve, every loaded `Beam.wavelength` is its own
  per-frame value).
- `dials.image_viewer` overlays on composite stills required an order-independent
  `_identifiers_for_frame` (the previous binary search assumed sorted frame order,
  which neither `_rebuild_shared_imageset_output` nor the multiprocessing-combine
  path guarantee); a `_still_scan_frame_list` helper that recovers absolute frame
  indices from per-experiment scans (the shared imageset's scan loses its start
  offset on JSON round-trip); per-overlay color reassignment for composite-id
  reflection tables; and a per-frame beam via `image_set.get_beam(self.index)`.
- `dials.split_experiments` → `dials.combine_experiments` round-trip required
  decoupling `to_dict()` beam re-consolidation from the existing `len(beams) > 1`
  guard (so single-experiment and single-shared-beam stills writes also re-consolidate
  to `XFELBeam`), and passing the full polarization tuple through
  `BeamFactory.make_xfel_beam` at re-consolidation (a silent drop of the LCLS ~0.999
  polarization fraction was a P0 regression caught in the review pass).

**The pattern downstream tools must follow** is: if a tool builds experiments at
runtime *without* going through `decode()`, it must run
`expt.beam = expt.get_monochromatic_beam()` after assembly.
`Experiment.get_monochromatic_beam()` is a no-op for non-XFEL beams, so the call is
safe to make unconditionally. For tools that load experiments via `decode()`, no
change is required.

**Verified downstream tools.** The following `cctbx.xfel` consumers have been
exercised end-to-end against output produced by this branch and produce correct
results:

- `cctbx.xfel.filter_experiments_by_rmsd`
- `cctbx.xfel.detector_residuals`
- `cctbx.xfel.merge`

This is the main XFEL consumer toolchain, and the auto-resolve in `decode()` means
no per-tool patches were required after the review hardening pass.

## 8. Open architectural questions

- **Beam refinement for XFEL stills is not designed.** Currently `refine_level0.phil`
  fixes the beam (`beam.fix=*all`) and the write-side re-consolidation produces one
  `XFELBeam` as expected. A phil that refines beam parameters would produce N
  differing beams that the re-consolidation guard would refuse to collapse — the
  result would be a much larger file with no clear consumer story. Beam refinement
  for XFEL stills needs separate design before the group should consider this
  contract for production.
- **`combine_all_ranks` output semantics.** The new phil option merges per-rank (MPI)
  or per-worker (multiprocessing) composite files into a single combined file with
  shared models. The current semantics — including how merging interacts with
  per-source-file model identity — have not been settled with the group.

## 9. Pre-publication checklist

Before these branches are pushed:

1. Run the pinned `black` and `ruff` (matching the pre-commit configuration in each
   repo) across all three repos so the published diff is canonical style only, with
   no formatter churn that would distract reviewers from the architectural diff.

## 10. Testing: fixture alignment, not production weakening

The DIALS test suite contains fixtures that pre-date this branch and construct
`Experiment` objects in ways that no longer match the new contract — most often
by leaving `expt.imageset = None` and supplying only `{beam, detector, crystal}`.
The contract change in §4 (stills detection via scan, not imageset type) and the
auto-routed sort path added to `dials.combine_experiments` in §5 make
`expt.imageset.paths()[0]` an unconditional call on the still-detection path, so
fixtures whose experiments lack an imageset surface as `AttributeError: 'NoneType'
object has no attribute 'paths'` at first contact with the new code.

The policy on this branch is: **update the fixture, not the production code.**
Weakening upstream call sites with `if iset is None: ...` guards would let
genuinely malformed experiment data — produced by some future bug — propagate
silently through the pipeline; the loud crash at the boundary is the correct
behavior. The proper fix is to bring the fixture into alignment with the new
contract that every experiment carries an imageset.

**Two fixtures have been updated to date:**

- `tests/algorithms/profile_model/ellipsoid/conftest.py` — the `test_experiment`
  fixture previously constructed `Experiment(beam=..., crystal=..., scan=...,
  detector=...)` with no imageset. Rebuilt to construct a real `ImageSequence`
  with a `Format.Reader` (dummy filename, N copies for N frames), then build
  `Experiment` from the imageset's own models. Commit `ae871b468` (2026-05-21).

- `tests/command_line/test_combine_experiments.py::test_min_max_reflections_per_experiment`
  — the test loads `multi_stills_combined.json` from `dials_data`, which is an
  old-format file with no imagesets (the `"imageset"` section is empty;
  experiments reference only beam/detector/crystal). The fixture file itself is
  not editable (it ships with `dials_data`), so the rewrite happens at test
  setup. A new helper `_attach_still_imagesets(json_path, tmp_path)` loads the
  raw `.expt`, builds **one** shared zero-oscillation `ImageSequence` (dummy
  reader) covering all N frames, then for each experiment attaches that shared
  imageset and a single-frame per-experiment scan pointing at its frame
  position. The rewritten `.expt` is written to `tmp_path` and that path is
  what the `dials.combine_experiments` subprocess reads. Commit `963d588ea`
  (2026-05-26).

**The shared-imageset trick.** The natural reflex when each experiment "needs
an imageset" is to build a per-experiment `ImageSequence`. That fails on this
branch: `_consolidate_stills_imagesets` sees N distinct imageset Python objects
all referencing the same source path, fails the `id()` identity check, and tries
to re-open the source file via `get_format_class_for_file(path).get_imageset([path])`
to build a consolidated sequence. With a dummy filename that crashes. Making
**all** experiments share **one** `ImageSequence` Python object satisfies the
identity check (line 370 of `combine_experiments.py`), so consolidation is a
genuine no-op and the dummy file is never opened. This matches the in-flight
state on this branch: `stills_process` itself produces composite output where
all experiments share one shared `ImageSequence` per source file (§5a).

**What did not change:**

- The test bodies themselves — the subprocess invocation, the expected results
  matrix, the `parametrize` cases, the `assert` statements — are unchanged. The
  test still measures what it always measured (the min/max reflection filter on
  the same crystallographic input).
- No production code under `src/` was modified to accommodate the test. The
  fix lives entirely under `tests/`.
- The `dials_data` fixtures themselves are unmodified. The rewrite is a
  test-side preprocessing step that produces a transient `.expt` in `tmp_path`.

**The pattern for future test failures with this signature** (`AttributeError`
on `iset.paths()` or any similar `None.x` on a model accessor in a still path)
is the same: write a helper that takes the existing fixture, attaches a shared
`ImageSequence` with the right beam/detector/scan structure, and feeds the
rewritten file to whatever code path the test exercises. The test's intent and
assertions stay intact; only the input is brought into line with the contract.

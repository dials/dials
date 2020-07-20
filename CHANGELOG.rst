DIALS 3.0.4 (2020-07-20)
========================

- ``dials.scale``: Allow usage of ``mode=image_group`` with ``filtering.method=deltacchalf`` when
  only providing a single data set (#1334)
- ``dials.import``: When using a template and specifying an image_range, missing images outside of
  the range will not cause a failure (#1333)
- ``dials.stills_process``: Show better error message in specific spotfinding failure case (#1180)


DIALS 3.0.3 (2020-07-06)
========================

Features
--------

- Developer tool: On posix systems, sending SIGUSR2 to DIALS commands will now print a stack trace (#1277)

Bugfixes
--------
- HTML reports: Plot bin centres instead bin minimum for d_min line plots vs. resolution (#1323)
- ``dials.export``: Fix inconsistency in mtz export when given a non-reference (e.g. I2 or primitive) setting (#1279)
- ``dials.refine_bravais_settings``: Fix crash with large (>2gb) reflection tables and reduce memory use (#1274)
- ``dials.scale``: Fix bug in outlier rejection code causing misidentification of outliers (with outlier_rejection=standard).
- ``dials.scale``: Fix outlier rejection formula to avoid overconfidence in spuriously low values


DIALS 3.0.2 (2020-06-23)
========================

Bugfixes
--------

- Fix crash in scaling error model handling (#1243)


DIALS 3.0.1 (2020-06-11)
========================

Features
--------

- dials.reciprocal_lattice_viewer: Add an option to show lattice(s) in the crystal rather than laboratory frame. (#1259)
- Support for mtz project_name in export and scaling

Bugfixes
--------

- dials.reciprocal_lattice_viewer: fix multiple experiment view for integrated data (#1284)


DIALS 3.0 (2020-05-22)
======================

Features
--------

- Show more useful output when crashing in C++ code (#659)
- dials.image_viewer: for the unit cell tool, rename parameters for consistency and add a new show_hkl option to filter displayed powder rings to select only those of interest. (#1192)
- In dials.integrate: changed the background box size multiplier to be a parameter (sigma_b_multiplier) - setting to small values significantly reduces memory requirements. (#1195)
- dials.image_viewer: add an overlaying showing pixels marked as strong by the spot-finding operations. That is, the pixels picked out by the "threshold" image. (#1200)
- dials.scale report file was renamed from scaling.html to dials.scale.html
  dials.symmetry report file was renamed from dials-symmetry.html to dials.symmetry.html (#1202)
- dials.report output file was renamed from dials-report.html to dials.report.html (#1206)
- dials.image_viewer: faster navigation between different image types. (#1213)
- Crystal model now has a new recalculated_unit_cell attribute. This allows it to store
  a post-refined unit cell (e.g. from dials.two_theta_refine) in addition to that from
  traditional geometry refinement (which was used for prediction). Downstream programs
  such as dials.scale and dials.export will now use the recalculated unit cell 
  where appropriate. (#1214)
- New best_monoclinic_beta parameter for dials.refine_bravais_settings and dials.symmetry.
  Setting this to False will ensure that C2 is selected in preference to I2, where I2
  would lead to a less oblique cell (i.e. smaller beta angle). (#1226)
- New scaling model, model=dose_decay, implementing a shared exponential decay component for multicrystal experiments (#1183)


Bugfixes
--------

- Avoid empty "Unable to handle" messages on failed dials.import (#600)
- Functions from dials.export now raise exceptions on errors rather than exit. This improves their use elsewhere (such as in dials.scale). (#1205)
- Ensure dials.index chooses the C2 setting with the smallest beta angle (#1217)
- Fix propagation of best_unit_cell and application of resolution cutoffs in dials.scale and export_mtz.
  Add a new mtz.best_unit_cell parameter to dials.export (#1248)
- Make some of the DIALS tools furthest downstream (``dials.scale``, ``dials.symmetry``, ``dials.merge`` and ``dials.report``) more robust in the case of very few reflections. (#1263)


Misc
----

- #1221


DIALS 2.2 (2020-03-15)
======================

Features
--------

- dials.image_viewer: Add a choice between displaying the raw or the corrected image. (#634)
- Automatically choose between the `simple` and `glm` background determination
  algorithms, depending on whether the detector appears to be integrating or
  counting. (#706)
- Allow adjustment of font size for overlay text, such as Miller indices and
  resolution ring values. (#1074)
- Keep goniometer and scan objects in indexing of still data, if provided in input (#1076)
- Experimental: ``dials.image_viewer`` can be remotely controlled via a
  ZeroMQ endpoint with the ``zmq_endpoint`` PHIL parameter. Initially,
  the viewer can be commanded to load new images. This requires the
  (optional) ``pyzmq``package. (#1085)
- Programs now generate a unique identifier for each experiment created, and reflection tables are linked via the experiment_identifiers map (#1086)
- Introduce `dials.anvil_correction` to correct the absorption of the incident and diffracted X-ray beam by the diamond anvils in a pressure cell.
  Call `dials.anvil_correction` on the output of `dials.integrate` and then proceed to use post-integration tools as normal, just as though the sample had been measured in air. (#1090)
- Map of detector efficiency for photon counting detectors as a function of 
  detector position added to report, based on the qe value applied at the end 
  of integration. (#1108)
- Significantly reduce the amount of memory required to write .refl output files (#1115)
- Add maximum_trusted_value=N option to spot finding to temporarily allow override of trusted range, e.g. to find overloaded spots in spot finding. (#1157)
- array_family.flex interface has changed: background and centroid algorithms are
  set via public properties. Instead of flex.strategy use functools.partial with
  the same signature. as_miller_array() raises KeyError instead of Sorry.
  .extract_shoeboxes() lost its verbosity parameter, use log levels instead. (#1158)
- dials.stills_process now supports imagesets of length > 1 (e.g. grid scans) (#1174)


Bugfixes
--------

- Fixed prediction on images numbered zero, so integrating works correctly. (#1128)
- Fix an issue (#1097) whereby aggregating small numbers of reflections into resolution bins could sometimes result in empty bins and consequent errors. (#1130)
- Ensure that restraints are ignored for parameterisations that are anyway fixed (#1142)
- Fix dials.search_beam_centre to ensure that the correct detector models are
  output when multiple detector models are present in the input.
  Fix dials.search_beam_centre n_macro_cycles option (previously it was starting
  from the original geometry every macro cycle). (#1145)
- dials.find_spots_server no longer slows down 3x when using resolution filters (#1170)


Misc
----

- #932, #1034, #1050, #1077


DIALS 2.1 (2019-12-12)
======================

Features
--------

- We now fully support Python 3 environments.
- MessagePack is now the default reflection table file format. Temporarily, the
  environment variable ``DIALS_USE_PICKLE`` can be used to revert to the previous
  pickle-based format, however this will be removed in a future version. (#986)
- new option for dials.show 'show_shared_models=True' displays which beam, crystal, and detector models are used across experiments (#996)
- Import still image sequence as N experiments dereferencing into one image set
  rather than one experiment. (#1014)
- Add `reflection_table.get` method for defaulted column access (#1031)


Bugfixes
--------

- Don't use -2 to indicate masked pixels, except for DECTRIS detectors where this
  is to be expected. (#536)
- No longer show pixels that are above the trusted range upper bound as
  "saturated" on the "variance" image. (#846)
- Correctly account for scan-varying crystals while providing a scan range to
  dials.integrate (#962)
- Ensure that generated masks do not include pixels that are overloaded on a few
  images, but only pixels that are always outside the trusted range. (#978)
- Rewritten parameter auto-reduction code for dials.refine provides finer-grained
  fixing of individual parameters rather than whole parameterisations and
  correctly takes constrained parameters into account (#990)
- Fix output of predictions in dials.refine.
  A recently-introduced bug meant that the updated predictions weren't
  being copied to the output reflections file. (#991)
- Allow scan-varying refinement where either the crystal cell or
  orientation is fixed. (#999)
- Respect batch= option to dials.symmetry - can reduce time taken for finding
  the symmetry for large data sets. (#1000)
- Scan-varying refinement no longer fails when the scan is wider than the
  observed reflections (e.g. when the crystal has died). Instead, the scan
  is first trimmed to match the range of the diffraction. (#1025)
- If convert_sequences_to_stills then delete the goniometer and scan. (#1035)
- Correctly account for scan-varying crystals in dials.slice_sequence (#1040)
- Eliminate systematic absences before applying change of basis op to minimum 
  cell in dials.symmetry. (#1064)


Improved Documentation
----------------------

- Add "Extending DIALS" page to developer documentation (#893)


Deprecations and Removals
-------------------------

- The command dials.analyse_output was removed.
  Its replacement, dials.report, will give you more useful output. (#1009)


Misc
----

- #983, #1004


DIALS 2.0 (2019-10-23)
======================

Features
--------

- Support exporting multi-dataset and still experiments to XDS_ASCII (#637)
- Replace default spotfinder with improved dispersion algorithm (#758)
- ``dials.report`` now displays oscillation data with units and more significant figures (#896)
- A new program, ``dials.sequence_to_stills`` to create split a sequence into a
  separate still Experiment for every scan point in the sequence, splitting
  reflections as necessary. (#917)
- Moved ``dials.export format=best`` to ``dials.export_best`` as that one needed
  access to the format object, the rest do not, and having ``dials.export`` work
  in the general case seems like a better idea... (#921)
- Unified logging output for dials programs - logs are no longer split into .log
  and .debug.log. Use -v to get debug output. (#923)
- New command ``dials.resolutionizer`` (replaces ``xia2.resolutionizer``). Add support for ``expt``/``refl``
  in ``dials.resolutionizer``. (#933)
- Changed the selection of reflections used for determination of the reflection
  profile parameters in integration. Now uses reflections which were previously
  used in refinement rather than all reflections, resulting in a speed
  improvement for large data sets and a negligible difference in the quality
  of the integrated results. (#942)
- ``dials.image_viewer`` now allows the choice between
  ``dispersion_extended`` (new default) and ``dispersion`` (old default)
  thresholding algorithms for investigating the effect of different
  spot-finding parameters. (#948)
- ``dials.rs_mapper`` now respects masked regions of images (including
  the trusted range mask). (#955)


Bugfixes
--------

- Fix and reinstate normalisation option in ``dials.option`` (#919)


Misc
----

- #795, #862, #895, #915, #924

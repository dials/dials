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

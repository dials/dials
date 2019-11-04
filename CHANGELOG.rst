Dials 2.0 (2019-10-23)
======================

Features
--------

- Support exporting multi-dataset and still experiments to XDS_ASCII (#637)
- Replace default spotfinder with improved dispersion algorithm (#758)
- ``dials.report`` now displays oscillation data with units and more significant figures (#896)
- A new program, ``dials.sequence_to_stills`` to create split a sequence into a
  separate still Experiment for every scan point in the sequence, splitting
  reflections as necessary. (#917)
- Moved `dials.export format=best` to `dials.export_best` as that one needed 
  access to the format object, the rest do not, and having `dials.export` work
  in the general case seems like a better idea... (#921)
- Unified logging output for dials programs - logs are no longer split into .log
  and .debug.log. Use -v to get debug output. (#923)
- New command dials.resolutionizer (replaces xia2.resolutionizer). Add support for expt/refl
  in dials.resolutionizer. (#933)
- Changed the selection of reflections used for determination of the reflection 
  profile parameters in integration. Now uses reflections which were previously
  used in refinement rather than all reflections, resulting in a speed 
  improvement for large data sets and a negligable difference in the quality
  of the integrated results. (#942)
- ``dials.image_viewer`` now allows the choice between
  `dispersion_extended` (new default) and `dispersion` (old default)
  thresholding algorithms for investigating the effect of different
  spot-finding parameters. (#948)
- dials.rs_mapper now respects masked regions of images (including
  the trusted range mask). (#955)


Bugfixes
--------

- Fix and reinstate normalisation option in dials.option (#919)


Misc
----

- #795, #862, #895, #915, #924

+++++++++++++++++++++++++++++++++++++
Processing Sweeps with Missing Images
+++++++++++++++++++++++++++++++++++++

DIALS treats sweeps as a contiguous set of rotation images, and in some circumstances deviations from this will cause problems. A clear example of this is having one or more "bad" images in your data set - simply removing them *will not* be enough to allow processing to complete smoothly.

Importing
=========

Importing the images with ``template=blah_####.cbf`` will not work, as this checks for continuous images, and it will be necessary to use ``allow_multiple_sweeps=true``. After this, finding spots and indexing work as usual, but with multiple imagesets.

Refinement
==========

Refinement *without* ``scan_varying=true`` will work fine, but following that, scan varying refinement will fail with an error::

  Sorry: A single scan-varying crystal model cannot be refined when associated with more than one scan or goniometer

The issue here is that scan-varying refinement requires that each crystal being refined is associated with a single scan. In our case, we have a single crystal model from indexing, but this is associated with multiple scans. To proceed, we can split the experiment list into individual files. This breaks the sharing of models between experiments by making copies of each model for each file::

  dials.split_experiments indexed.refl indexed.expt

From this point, we could process each block as per the usual tutorial instructions (ideally in separate directories). However, this will refine the beam and the detector, which *should* be shared, separately for each process. A better way to proceed would be to recombine the experiments as follows::

  dials.combine_experiments experiments_*.expt reflections_*.refl \
    reference_from_experiment.goniometer=0 \
    reference_from_experiment.detector=0 \
    reference_from_experiment.beam=0

The ``reference_from_experiment`` options tells ``dials.combine_experiments`` to take the goniometer, detector and beam models only from the first experiment (although any would have done as these are merely copies of each other). The resulting ``combined.expt`` file has a separate (but identical) crystal model for each scan, alongside shared goniometer, detector and beam models. Scan-varying refinement and integration can now proceed as usual::

  dials.refine combined.expt combined.refl scan_varying=true
  dials.integrate refined.expt refined.refl nproc=4

Export
======

Currently, ``dials.export`` will not allow MTZ export of the multiple-experiment integration, but we can rely on ``dials.split_experiments`` again::

  dials.split_experiments integrated.expt integrated.refl
  i=0
  for e in $(ls experiments_*.expt)
  do
    r=reflections_$i.refl
    echo "dials.export $e $r mtz.hklout=integrated_$i.mtz"
    ((i++))
  done

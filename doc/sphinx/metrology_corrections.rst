Refining multi-tile detector metrology with DIALS
=================================================

Introduction
------------

At the end of the :doc:`advanced_tutorial`, we showed plots from
:program:`dials.analyse_output`, including images of the positional
residuals as a function of position on the detector face. There were
clearly some systematic effects, suggesting whole-tile shifts or rotations.
Generally, these are small (typically much less than the pixel size) and are
expected. DECTRIS provides calibration tables from factory metrology
experiments, which XDS can use. Currently we are not applying the correction
tables to reflection positions in DIALS. In some situations we may have to work
with a much less well calibrated or characterised detector than the Pilatus
instruments installed at Diamond. Nevertheless, the multi-panel detector models
in DIALS allow us to *discover* the appropriate corrections assuming we have
good diffraction data to refine against.

Warning!
--------

In this tutorial we shall use the standard tutorial data from the I04 beamline,
hosted at |thaumatin|.

.. |thaumatin| image:: https://zenodo.org/badge/doi/10.5281/zenodo.10271.png
               :target: http://dx.doi.org/10.5281/zenodo.10271


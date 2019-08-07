Processing in Detail (Thaumatin)
================================

.. highlight:: none

Introduction
------------

DIALS processing may be performed by either running the individual tools (spot
finding, indexing, refinement, integration, exporting to MTZ) or you can run
:samp:`xia2 pipeline=dials`, which makes informed choices for you at each stage. In
this tutorial we will run through each of the steps in turn, checking the output
as we go. We will also enforce the correct lattice symmetry.

Tutorial data
-------------

The following example uses a Thaumatin dataset collected using beamline I04
at Diamond Light Source, which is available for download from |thaumatin|.

.. |thaumatin| image:: https://zenodo.org/badge/doi/10.5281/zenodo.10271.svg
               :target: https://doi.org/10.5281/zenodo.10271

Import
^^^^^^

The first stage of step-by-step DIALS processing is to import the data - all
that happens here is that the image headers are read, and a file describing
their contents (:ref:`datablock.expt <datablock-json>`) is written.

::

  dials.import data/th_8_2_0*cbf


The output just describes what the software understands of the images it was
passed, in this case one sweep of data containing 540 images.

.. literalinclude:: logs/dials.import.log

Find Spots
^^^^^^^^^^

The first "real" task in any processing using DIALS is the spot finding.
Here we request multiple processors to speed this up (:samp:`nproc=4`).
It still takes a little while because we are finding spots on every image in the
dataset. This reflects the modular philosophy of the DIALS toolkit and will
enable us to do global refinement later on.

.. literalinclude:: logs/dials.find_spots.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. literalinclude:: logs/dials.find_spots.log

The default parameters for :doc:`dials.find_spots<../programs/dials_find_spots>`
usually do a good job for Pilatus images, such as these. However they may
not be optimal for data from other detector types, such as CCDs or image
plates. Issues with incorrectly set gain or sigma thresholds might lead to
far too many spots being extracted (for example). It is always worth
inspecting the images with :doc:`dials.image_viewer<../programs/dials_image_viewer>`,
especially if you are having issues with spot finding::

  dials.image_viewer datablock.expt

Viewing the various images from 'image' to 'threshold' gives an idea of how
the various parameters affect the spot finding algorithm. The final image,
'threshold' is the one on which spots are found, so ensuring this produces
peaks at real diffraction spot positions will give the best chance of success.

Having found strong spots it is worth checking the image viewer again::

  dials.image_viewer datablock.expt strong.refl

The :doc:`dials.image_viewer<../programs/dials_image_viewer>` tool is not as
fast as viewers such as ADXV, however it does integrate well with DIALS data
files. Information about the beam centre, spot centroids, reflection
shoeboxes and other data stored in the pickle files created by DIALS
programs can be overlaid on the diffraction images. You may need to adjust
the colour scheme and brightness to get the best out of it. A brightness of
20 with the 'invert' colour scheme works well with this data. Move forward a
few images to find a spot whose complete rocking curve is recorded. The
highest valued pixel in that three dimensional spot is marked with a pink
dot. The spot centre of mass is shown by a red cross. This is usually close to the
peak pixel, but slightly offset as the centroid algorithm allows calculation
of the spot centre at a better precision than the pixel size and image
angular 'width'. The strong pixels marked as being part of the peak are
highlighted with a green dot. The reflection shoebox you see with a blue
border is the smallest three dimensional box that can contain the continuous
peak region, that is, there is no background border region displayed here.

.. image:: /figures/found_spot.png

Another very powerful tool for investigating problems with strong spot positions
is :doc:`dials.reciprocal_lattice_viewer<../programs/dials_reciprocal_lattice_viewer>`.
This displays the strong spots in 3D, after mapping them from their detector
positions to reciprocal space. In a favourable case like this you should be
able to see the crystal's reciprocal lattice by eye in the strong spot
positions. Some practice may be needed in rotating the lattice to an
orientation that shows off the periodicity in reciprocal lattice positions::

  dials.reciprocal_lattice_viewer datablock.expt strong.refl

.. image:: /figures/reciprocal_lattice_strong.png

Indexing
^^^^^^^^

The next step will be indexing of the strong spots by
:doc:`dials.index<../programs/dials_index>`, which by default uses a
3D FFT algorithm, although the 1D FFT algorithm can be selected using the
parameter :samp:`indexing.method=fft1d`. We will pass in all the strong
spots found in the dataset - so no need to select subsets of images widely
separated in :math:`\phi`.

.. literalinclude:: logs/dials.index.cmd

If known, the space group and unit cell can be
provided at this stage using the :samp:`space_group` and :samp:`unit_cell`
parameters, otherwise indexing and refinement will be carried out in the
primitive lattice using space group P1.

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. literalinclude:: logs/dials.index.log

It is worth reading through this output to understand what the indexing
program has done. Note that this log is automatically captured in the file
:file:`dials.index.log`. There is also a somewhat more information written
into :file:`dials.index.debug.log`, but this is probably only helpful if
something has gone wrong and you are trying to track down why.

Inspecting the log shows that the indexing step is done at fairly low
resolution: ``Setting d_min: 3.89``. The resolution limit of data that can
be used in indexing is determined by the size of the 3D FFT grid and the
likely maximum cell dimension. Here we used the default :math:`256^3` grid
points: ``FFT gridding: (256,256,256)``. What follows then are macrocycles of
refinement at increasing resolution to bootstrap the indexing solution to as
many of the strong reflections as possible. In each case you can see that
only 8099 reflections are used in the refinement job. The diffraction
geometry is here described by only 16 parameters (6 for the detector, 1 beam
angle, 3 crystal 'misset' angles and 6 triclinic cell parameters). The
problem is thus hugely overdetermined. In order to save time, the refinement here
uses a subset of the input reflections, by default using 100 reflections for
every degree of the scan.

Continuing to look through the log, we see that the first macrocyle of
refinement makes a big improvement in the positional RMSDs. Second and
subsequent macrocycles are refined using the same number of reflections, but
after extending to higher resolution. The RMSDs at the start of each cycle
start off worse than at the end of the previous cycle, because the best fit
model for lower resolution data is being applied to higher resolution
reflections. As long as each macrocyle shows a reduction in RMSDs then
refinement is doing its job of extending the applicability of the model out to
a new resolution limit, until eventually the highest resolution strong
spots have been included. The final macrocycle includes data out to 1.26
Angstroms and produces a final model with
RMSDs of 0.044 mm in X, 0.036 mm in Y and 0.021 degrees in :math:`\phi`,
corresponding to 0.25 pixels in X, 0.21 pixels in Y and 0.14 image widths in
:math:`\phi`.

Despite the high quality of this data, we notice from the log that at each
macrocycle there were some outliers identified and removed from
refinement as resolution increases. Large outliers can dominate refinement
using a least squares target, so it is important to be able to remove these.
More about this is discussed below in :ref:`sec-refinement`.

After indexing it can be useful to inspect the reciprocal lattice again::

  dials.reciprocal_lattice_viewer indexed.expt indexed.refl

Now indexed/unindexed spots are differentiated by colour, and it is possible
to see which spots were marked by :doc:`dials.refine <../programs/dials_refine>`
as outliers. If you have a dataset with multiple lattices present, it may be
possible to spot them in the unindexed reflections.

If you want to specify the Bravais lattice for processing (i.e. include the
lattice constraints in the refinement) then you need to either request it at
this stage using :samp:`space_group=P4` as a command-line option to
:doc:`dials.index <../programs/dials_index>` or you can use
:doc:`dials.refine_bravais_settings<../programs/dials_refine_bravais_settings>`,
which will take the results of the P1 autoindexing and run refinement with
all of the possible Bravais settings applied - after which you may select
the preferred solution.

.. literalinclude:: logs/dials.refine_bravais_settings.cmd

gives a table containing scoring data and unit cell for
each Bravais setting. The scores include the metric fit (in degrees),
RMSDs (in mm), and the best and worse correlation coefficients for data
related by symmetry elements implied by the lowest symmetry space group from the
Bravais setting. This uses the raw spot intensity measurement from the
spot-finding procedure (uncorrected and unscaled) but provides a very
useful check to see if the data does appear to adhere to the proposed
symmetry operators.

.. literalinclude:: logs/dials.refine_bravais_settings.log

In this example we would continue processing (i.e. proceed to the refinement
step, perhaps) with :samp:`bravais_setting_9.expt`. Sometimes (that is, when
the change of basis operator is not equal to :samp:`a,b,c`) it is
necessary to reindex the :ref:`indexed.refl <reflection_pickle>` file output
by :doc:`dials.index<../programs/dials_index>`.
In this case as the change of basis operator to the chosen setting
is the identity operator (:samp:`a,b,c`) this step is not needed. We run it
anyway to demonstrate its use.

.. literalinclude:: logs/dials.reindex.cmd

This outputs the file :file:`reindexed.refl` which should be
used as input to downstream programs in place of :file:`indexed.refl`.

.. _sec-refinement:

Refinement
^^^^^^^^^^

The model is already refined during indexing, but we can also add explicit
refinement steps using :doc:`dials.refine <../programs/dials_refine>`
in here, to use all reflections in refinement rather than a subset and to
fit a scan-varying model of the crystal. There
are many options to refinement. As an
aside, to show all the options up to and including ``expert_level=1``
use this command::

  dials.refine -c -e 1

Equivalent command-line options exist for all the main DIALS programs. To
refine a static model including the tetragonal constraints we just do:

.. literalinclude:: logs/dials.refine.cmd

This used all reflections in refinement rather than a subset and provided a
small reduction in RMSDs. However, the refined model is still static over
the whole dataset. We may want to do an additional refinement job to fit a
more sophisticated model for the crystal, allowing small misset rotations to
occur over the course of the scan. There are usually even small changes to
the cell dimensions (typically resulting in a net increase in cell volume)
caused by exposure to radiation during data collection. To account for both
of these effects we can extend our parameterisation to obtain a smoothed
'scan-varying' model for both the crystal orientation and unit cell. This means
running a further refinement job starting from the output of the
previous job:

.. literalinclude:: logs/dials.sv_refine.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. literalinclude:: logs/dials.sv_refine.log

In this case we didn't alter the default choices that affect scan-varying
refinement, the most important of which is the number of intervals into which
the full scan is divided. This determines the number of samples that will be
used by the Gaussian smoother. More samples allows sharper changes to the model,
but overdoing this will lead to unphysical changes to the model that are just
fitting noise in the data. Figuring out the optimum number of points to use
is challenging. Here we are happy with the default interval width of 36 degrees
(this is a parameter at ``expert_level=1``).

See :ref:`html-report` for further information on how to view the smoothly
varying crystal cell parameters using :samp:`dials.report`.

Integration
^^^^^^^^^^^

After the refinement is done the next step is integration, which is performed
by the program :doc:`dials.integrate <../programs/dials_integrate>`. Mostly,
the default parameters are fine for Pilatus data, which will perform
XDS-like 3D profile fitting while using a generalized linear model in order
to fit a Poisson-distributed background model. We will also increase the
number of processors used to speed the job up.

.. literalinclude:: logs/dials.integrate.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. literalinclude:: logs/dials.integrate.log

Checking the log output, we see that after loading in the reference
reflections from :file:`refined.refl`, new predictions are made up to the
highest resolution at the corner of the detector. This is fine, but if we
wanted to we could have adjusted the resolution limits using parameters
:samp:`prediction.d_min` and :samp:`prediction.d_max`. The predictions are
made using the scan-varying crystal model recorded in
:file:`refined.expt`. This ensures that prediction is made using
the smoothly varying lattice and orientation that we determined in the
refinement step. As this scan-varying model was determined in advance of
integration, each of the integration jobs is independent and we can take
advantage of true parallelism during processing.

The profile model is calculated from the reflections in
:file:`refined.refl`. First reflections with a too small 'zeta'
factor are filtered out. This essentially removes reflections that are too
close to the spindle axis. In general these reflections require significant
Lorentz corrections and as a result have less trustworthy intensities anyway.
From the remaining reflection shoeboxes, the average beam divergence and
reflecting range is calculated, providing the two Gaussian width parameters
:math:`\sigma_D` and :math:`\sigma_M` used in the 3D profile model.

Following this, independent integration jobs are set up. These jobs
overlap, so reflections are assigned to one or more jobs. What follows are
blocks of information specific to each integration job.

After these jobs are finished, the reflections are 'post-processed', which
includes the application of the LP correction to the intensities. Then
summary tables are printed giving quality statistics first by frame, and
then by resolution bin.


.. _html-report:

HTML report
^^^^^^^^^^^

Much more information from the various steps of data processing can be found
within an HTML report generated using the program
:doc:`dials.report <../programs/dials_report>`.
This is run simply with:

.. literalinclude:: logs/dials.report.cmd

which produces the file :download:`dials-report.html <logs/dials-report.html>`.

This report includes plots showing the scan-varying crystal orientation and
unit cell parameters. The latter of these
is useful to check that changes to the cell during processing appear reasonable.
In this tutorial, we see an overall increase in all three cell parameters,
however the greatest change, in lengths *a* and *b*, is only about 0.02 Angstroms. If
significant cell volume increases had been observed that might be indicative of
radiation damage. However we can't yet conclude that there is *no* radiation
damage from the *lack* of considerable change observed. We can at least see from
this and the low final refined RMSDs that this is a very well-behaved dataset.

Some of the most useful plots are

* **Difference between observed and calculated centroids vs phi**,
  which shows how the average
  residuals in each of X, Y, and :math:`\phi` vary as a fuction of :math:`\phi`.
  If scan-varying refinement has been successful in capturing the real changes
  during the scan then we would expect these plots to be straight lines.

* **Centroid residuals in X and Y**, in which the X, Y residuals are shown
  directly. The key point here is to look for a globular shape centred at the origin.

* **Difference between observed and calculated centroids in X and Y**,
  which show the difference between predicted and observed reflection positions
  in either X or Y as functions of detector position. From these plots it is very
  easy to see whole tiles that are worse than their neighbours, and whether
  those tiles might be simply shifted or slightly rotated compared to the model
  detector.

* **Reflection and reference correlations binned in X/Y**.
  These are useful companions to the
  plots of centroid residual as a function of detector position above.
  Whereas the above plots show systematic errors in the positions and
  orientations of tiles of a multi-panel detector, these plots indicate what
  effect that (and any other position-specific systematic error) has on the
  integrated data quality. The first of these plots shows the correlation
  between reflections and their reference profiles for all reflections in the
  dataset. The second shows only the correlations between the strong reference
  reflections and their profiles (thus these are expected to be higher and do
  not extend to such high resolution).

* **Distribution of I/Sigma vs Z**. This reproduces the
  :math:`\frac{I}{\sigma_I}` information versus frame number given in the log
  file in a graphical form. Here we see that :math:`\frac{I}{\sigma_I}` is fairly
  flat over the whole dataset, which we might use as an indication that there
  were no bad frames, not much radiation damage occurred and that scale factors
  are likely to be fairly uniform.

Exporting as MTZ
^^^^^^^^^^^^^^^^

The final step of dials processing is to export the integrated results to mtz
format, suitable for input to downstream processing programs such as pointless_
and aimless_.

.. literalinclude:: logs/dials.export.cmd

And this is the output, showing the reflection file statistics.

.. literalinclude:: logs/dials.export.log

What to do Next
---------------

The following demonstrates how to take the output of dials processing and
continue with downstream analysis using the CCP4 programs pointless_, to sort the data and assign
the correct symmetry, followed by scaling with aimless_ and intensity analysis
using ctruncate_::

  pointless hklin integrated.mtz hklout sorted.mtz > pointless.log
  aimless hklin sorted.mtz hklout scaled.mtz > aimless.log << EOF
  resolution 1.3
  anomalous off
  EOF
  ctruncate -hklin scaled.mtz -hklout truncated.mtz \
  -colin '/*/*/[IMEAN,SIGIMEAN]' > ctruncate.log

to get merged data for downstream analysis. The output from this includes
the merging statistics which will give a better idea about data quality. It is
easiest to view these logfiles using the program :program:`logview`, e.g.::

  logview aimless.log

Often passing in a sensible resolution limit to aimless is helpful. Here we
assumed we ran first without a resolution limit to help decide where to cut
the data. This indicated slightly anisotropic diffraction, with diffraction along
the *c*\* direction a little better than *a*\* and *b*\* directions, which are
equivalent. Diffraction quality is good, however completeness falls off sharply,
especially in the *c*\* direction. Following this we chose to exclude all data
at a resolution higher than 1.3 Angstroms, to ensure about 80% completeness in
the outer shell. Here is the summary from aimless.log:

::

  Summary data for        Project: DIALS Crystal: XTAL Dataset: FROMDIALS

                                             Overall  InnerShell  OuterShell
  Low resolution limit                      150.00    150.00      1.32
  High resolution limit                       1.30      7.12      1.30

  Rmerge  (within I+/I-)                     0.061     0.024     0.416
  Rmerge  (all I+ and I-)                    0.069     0.026     0.488
  Rmeas (within I+/I-)                       0.075     0.030     0.575
  Rmeas (all I+ & I-)                        0.076     0.030     0.610
  Rpim (within I+/I-)                        0.043     0.017     0.395
  Rpim (all I+ & I-)                         0.033     0.014     0.358
  Rmerge in top intensity bin                0.029        -         -
  Total number of observations              308123      2257      5493
  Total number unique                        62352       499      2474
  Mean((I)/sd(I))                             10.7      27.1       1.4
  Mn(I) half-set correlation CC(1/2)         0.999     0.999     0.722
  Completeness                                98.2      99.8      80.1
  Multiplicity                                 4.9       4.5       2.2

  Anomalous completeness                      92.3     100.0      47.8
  Anomalous multiplicity                       2.4       3.0       1.5
  DelAnom correlation between half-sets     -0.002     0.279     0.065
  Mid-Slope of Anom Normal Probability       0.953       -         -


.. _pointless: http://www.ccp4.ac.uk/html/pointless.html
.. _aimless: http://www.ccp4.ac.uk/html/aimless.html
.. _ctruncate: http://www.ccp4.ac.uk/html/ctruncate.html

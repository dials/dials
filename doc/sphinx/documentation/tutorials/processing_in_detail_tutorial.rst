Processing in Detail
====================

Introduction
------------

DIALS processing may be performed by either running the individual tools (spot
finding, indexing, refinement, integration, exporting to MTZ) or you can run
:samp:`xia2 -dials`, which makes informed choices for you at each stage. In
this tutorial we will run through each of the steps in turn, checking the output
as we go. We will also enforce the correct lattice symmetry.

Tutorial data
-------------

The following example uses a Thaumatin dataset collected using beamline I04
at Diamond Light Source which is available for download from |thaumatin|.

.. |thaumatin| image:: https://zenodo.org/badge/doi/10.5281/zenodo.10271.png
               :target: http://dx.doi.org/10.5281/zenodo.10271

Import
^^^^^^

The first stage of step-by-step DIALS processing is to import the data - all
that happens here is that the image headers are read, and a file describing
their contents (:ref:`datablock.json <datablock-json>`) is written. It's worth noting that if
this file is changed subsequent processing can
use this.

::

  dials.import data/th_8_2_0*cbf

The output just describes what the software understands of the images it was
passed, in this case one sweep of data containing 540 images.

.. literalinclude:: logs/dials.import.log

Find Spots
^^^^^^^^^^

The first "real" task in any DIALS processing will be the spot finding.
Here we request multiple processors to speed up the spot-finding (:samp:`nproc=4`).
It takes a little while because we are finding spots on every image in the
dataset. This reflects the modular philosophy of the DIALS toolkit and will
enable us to do global refinement later on.

::

  dials.find_spots datablock.json nproc=4

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. literalinclude:: logs/dials.find_spots.log

The default parameters for :doc:`dials.find_spots<../programs/dials_find_spots>`
usually do a good job for Pilatus images, such as these. However they may
not be optimal for data from other detector types, such as CCDs or image
plates. Issues with incorrectly set gain or sigma thresholds might lead to
far too many spots being extracted (for example). It is always worth
inspecting the images with :program:`dials.image_viewer`, especially if you
are having issues with spot finding::

  dials.image_viewer datablock.json

Viewing the various images from 'image' to 'threshold' gives an idea of how
the various parameters affect the spot finding algorithm. The final image,
'threshold' is the one on which spots are found, so ensuring this produces
peaks at real diffraction spot positions will give the best chance of success.

Having found strong spots it is worth checking the image viewer again::

  dials.image_viewer datablock.json strong.pickle

The :program:`dials.image_viewer` tool is not as fast as tools such as ADXV,
however it does integrate well with DIALS data files. Information about the
beam centre, spot centroids, reflection shoeboxes and other data stored in
the pickle files created by DIALS programs can be overlaid on the
diffraction images. You may need to adjust the colour scheme and brightness
to get the best out of it. A brightness of 20 with the 'invert' colour
scheme works well with this data. Move forward a few images to find a spot
whose complete rocking curve is recorded. The highest valued pixel in that
three dimensional spot is marked with a pink dot. The spot centre of mass is
a red cross. This is usually close to the peak pixel, but slightly offset as
the centroid algorithm allows calculation of the spot centre at a better
precision than the pixel size and image angular 'width'. The strong pixels
marked as being part of the peak are highlighted with a green dot. The
reflection shoebox you see with a blue border is the smallest three
dimensional box that can contain the continuous peak region, that is, there
is no background border region displayed here.

.. image:: /figures/found_spot.png

Indexing
^^^^^^^^

The next step will be indexing of the strong spots, which by default uses a 3D FFT
algorithm, although the 1D FFT algorithm can be selected using the parameter
:samp:`indexing.method=fft1d`. We will pass in all the strong spots found in
the dataset - so no need to select subsets of images widely separated in
:math:`\phi`.

::

  dials.index datablock.json strong.pickle

If known, the space group and unit cell can be
provided at this stage using the :samp:`space_group` and :samp:`unit_cell`
parameters, otherwise indexing and refinement will be carried out in the
primitive lattice using space group P1.

.. literalinclude:: logs/dials.index.log

It is worth looking through this output to understand what the indexing program
has done. Note that this log
is automatically captured in the file :file:`dials.index.log`. There is also
a somewhat more information written into :file:`dials.index.debug.log`, but
this is probably only helpful if something has gone wrong and you are trying
to track down why.

Inspecting the log shows that the indexing step is done at fairly low
resolution: ``Setting d_min: 3.89``. The resolution limit of data that
can be used in indexing is determined by the size of the 3D FFT grid and the
likely maximum cell dimension. Here we
used :math:`256^3` grid points: ``FFT gridding: (256,256,256)``.
What follows are four macrocycles
of refinement at increasing resolution to bootstrap the indexing solution to as
many of the strong reflections as possible. In each case you can see that only
8099 reflections are used in the refinement job. The diffraction geometry is
here described by only 16 parameters (6 for the detector, 1 beam angle, 3
crystal 'misset' angles and 6 triclinic cell parameters). The problem is thus
hugely overdetermined. In order to save time, refinement uses a subset of the
input reflections, by default using 100 reflections for every degree of the scan.

Continuing to look through the log, we see that the first macrocyle of refinement makes
a big improvement in the positional RMSDs. The second macrocycle includes more reflections, after
extending to 3.2 Angstroms. The current model now shows slightly worse RMSDs
at the start, now that the higher resolution reflections are included, but refinement reduces
these again.
A similar situation is observed on the third and fourth macrocycles.
The RMSDs start higher again, now that more reflections are included, but refinement
is able to drive these down a little.
The final macrocycle includes data out to 1.3 Angstroms and refinement produces
a final model with
RMSDs of 0.043 mm in X, 0.033 mm in Y and 0.019 degrees in :math:`\phi`, corresponding
to 0.25 pixels in X, 0.19 pixels in Y and 0.13 image widths in :math:`\phi`.

Despite the high quality of this data, we notice from the ``Summary statistics``
tables that there were some outliers identified and removed from
refinement as resolution increases.
In the final macrocyle, prior to outlier rejection, we see the
distribution of positional residuals in the Y direction is tight around the
median, except for extreme values both positive and negative of more than 1 mm.
The angular residuals show a similar pattern with half the data having residuals
of less than about 0.14 degrees from the predicted positions, but the extreme
is as much as 1.4 degrees from the predicted diffraction angle. Large outliers
can dominate refinement using a least squares target, so it is important
to be able to remove these.

If you want to specify the Bravais lattice for processing (i.e. include the
lattice constraints in the refinement) then you need to either specify this
lattice at this stage as

::

  space_group=P4

as a command-line option to :doc:`dials.index <../programs/dials_index>`
or you can use
:doc:`dials.refine_bravais_settings <../programs/dials_refine_bravais_settings>`,
which will take the results of the P1 autoindexing and run refinement with all
of the possible Bravais settings applied - after which you may select the
preferred solution.

::

  dials.refine_bravais_settings experiments.json indexed.pickle

gives a table containing scoring data and unit cell for
each Bravais setting. The scores include the the metric fit (in degrees),
RMSDs (in mm), and the best and worse correlation coefficients for data
related by symmetry elements implied by the lowest symmetry space group from the
Bravais setting. This uses the raw spot intensity measurement from the
spot-finding procedure (uncorrected and unscaled) but provides a very
useful check to see if the data does appear to adhere to the proposed
symmetry operators.

.. literalinclude:: logs/dials.refine_bravais_settings.log

In this example we would continue processing (i.e. proceed to the refinement
step, perhaps) with :samp:`bravais_setting_9.json`. Sometimes it may be
necessary to reindex the :ref:`indexed.pickle <reflection_pickle>` file output by dials.index.
However, in this case as the change of basis operator to the chosen setting
is the identity operator (:samp:`a,b,c`) this step is not needed. We run it
anyway to demonstrate its use::

  dials.reindex indexed.pickle change_of_basis_op=a,b,c

This outputs the file :ref:`reindexed_reflections.pickle <reflection_pickle>` which should be
used as input to downstream programs in place of :ref:`indexed.pickle <reflection_pickle>`.


Refinement
^^^^^^^^^^

Although the model is already refined during indexing we can also add an
explicit refinement
step using :doc:`dials.refine <../programs/dials_refine>` in here. There
are many options to refinement. As an
aside, to show all the options up to and including ``expert_level = 1``
use this command::

  dials.refine -c -e 1

Equivalent command-line options exist for all the main DIALS programs.

The main reason
we may want to do an additional refinement job is to use a more sophisticated
model for the crystal,
allowing small misset rotations to occur over the course of the scan.
There are usually even small changes to the
cell dimensions (typically resulting in a net increase in cell volume) caused
by exposure to radiation during data collection. To account for both of these
effects we can extend our parameterisation to obtain a smoothed 'scan-varying'
model for both the crystal orientation and unit cell. To do this, we run a
further refinement job starting from the output of the previous job::

  dials.refine bravais_setting_9.json indexed.pickle scan_varying=true

The output for this job is

.. literalinclude:: logs/dials.refine.log

In this case we didn't alter the default choices that affect scan-varying
refinement, the most important of which is the number of intervals into which
the full scan is divided. This determines the number of samples that will be
used by the Gaussian smoother. More samples allows sharper changes to the model,
but overdoing this will lead to unphysical changes to the model that are just
fitting noise in the data. Figuring out the optimum number of points to use
is challenging. Here we are happy with the default interval width of 36 degrees
(this is a parameter at ``expert_level = 1``).

To view the smoothly varying crystal cell parameters use the following command::

  dials.plot_scan_varying_crystal refined_experiments.json

This program creates a directory :file:`scan-varying_crystal` containing
plots :file:`orientation.png` and :file:`unit_cell.png`. The latter of these
is useful to check that changes to the cell during processing appear reasonable.

.. image:: /figures/unit_cell.png

We see an overall increase in all three cell parameters, however the greatest
change, in lengths *a* and *b*, is only about 0.02 Angstroms. If
significant cell volume increases had been observed that might be indicative of
radiation damage. However we can't yet conclude that there is *no* radiation
damage from the *lack* of considerable change observed. We can at least see from
this and the low final refined RMSDs that this is a very well-behaved dataset
though.

Integration
^^^^^^^^^^^

After the refinement is done the next step is integration, which is performed
by the program :doc:`dials.integrate <../programs/dials_integrate>`. Mostly, the
default parameters are fine, which will perform XDS-like 3D profile fitting. However,
for datasets with very weak background, such as this, the default :samp:`nsigma`
background outlier rejection algorithm tends to underestimate the real background
value. This is because that method is only really appropriate for values from
a normal distribution, which is a poor approximation for a Poisson distibution
with a small mean, and significant skewness. For this reason we switch off
all outlier rejection from the background calculation.

From checking the output of :samp:`dials.integrate -c` we see that the full
parameter to do this is given by :samp:`integration.background.simple.outlier.algorithm=null`
but partial string matching can be used for command line parameters when the
partial match is unambiguous. This saves a lot of typing!

We will also increase the number of processors used to speed the job up.

::

  dials.integrate refined_experiments.json refined.pickle \
  background.algorithm=glm nproc=4

The log file is quite long.

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. literalinclude:: logs/dials.integrate.log

Checking this output we see that after loading in the reference reflections
from :file:`refined.pickle`,
new predictions are made up to the highest resolution at the corner of the
detector. This is fine, but if we wanted to we could have adjusted the
resolution limits using parameters :samp:`dmin` and :samp:`dmax`. The predictions
are made using the scan-varying crystal model recorded in
:file:`refined_experiments.json`. This ensures that prediction is made using
the smoothly varying lattice and orientation that we determined in the refinement
step. As this scan-varying model was determined in advance of integration, each
of the integration jobs is independent and we can take advantage of true
parallelism during processing.

The profile model is then calculated from the reflections in
:file:`refined.pickle`. First reflections with a too small 'zeta'
factor are filtered out. This essentially removes reflections that are too
close to the spindle axis. In general these reflections require significant
Lorentz corrections and as a result have less trustworthy intensities anyway.
From the remaining reflection shoeboxes, the average beam divergence and
reflecting range is calculated, providing the two Guassian width parameters
:math:`\sigma_D` and :math:`\sigma_M` used in the 3D profile model.

Following this, the independent integration jobs are set up. These jobs overlap,
so reflections are assigned to one or more jobs. What follows are blocks of
information specific to each integration job.

After these jobs are finished, the reflections are 'post-processed', which includes
the application of the LP correction to the intensities. Then summary tables
are printed giving quality statistics first by frame, and then by resolution bin.

Graphical analysis of the output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Much more information is available from the integration output in graphical form
using the command

::

  dials.analyse_output integrated.pickle

By default the plots will be written into a new directory :file:`analysis` with
subdirectories for different types of analysis::

  analysis
  ├── background
  ├── centroid
  ├── intensity
  ├── reference
  └── strong

Some of the most useful plots are

* :file:`background/background_model_mean_vs_xy.png`, which shows the mean
  background value as a function of detector position.

* :file:`centroid/centroid_mean_diff_vs_phi.png`, which shows how the average
  residuals in each of X, Y, and :math:`\phi` vary as a fuction of :math:`\phi`.
  If scan-varying refinement has been successful in capturing the real changes
  during the scan then we would expect these plots to be straight lines.

  .. image:: /figures/centroid_mean_diff_vs_phi.png

* :file:`centroid/centroid_xy_residuals.png`, on which the X, Y residuals are shown
  directly. The key point here is to look for a globular shape centred at 0.0.

  .. image:: /figures/centroid_xy_residuals.png

* :file:`centroid/centroid_diff_x.png` and :file:`centroid/centroid_diff_y.png`,
  which show the difference between predicted and observed reflection positions
  in either X or Y as functions of detector position. From these plots it is very
  easy to see whole tiles that are worse than their neighbours, and either whether
  those tiles might be simply shifted or slightly rotated compared to the model
  detector.

  .. image:: /figures/centroid_diff_x.png

  .. image:: /figures/centroid_diff_y.png

* :file:`reference/reflection_corr_vs_xy.png` and
  :file:`reference/reference_corr_vs_xy.png`. These are useful companions to the
  plots of centroid residual as a function of detector position displayed above.
  Whereas the earlier plots show systematic errors in the positions and
  orientations of tiles of a multi-panel detector, these plots indicate what
  effect that (and any other position-specific systematic error) has on the
  integrated data quality. The first of these plots shows the correlation
  between reflections and their reference profiles for all reflections in the
  dataset. The second shows only the correlations between the strong reference
  reflections and their profiles (thus these are expected to be higher and do
  not extend to such high resolution). The first plot is probably the most
  useful, and that is reproduced here.

  .. image:: /figures/reflection_corr_vs_xy.png

* :file:`intensity/ioversigma_vs_z.png`. This reproduces the
  :math:`\frac{I}{\sigma_I}` information versus frame number given in the log
  file in a graphical form. Here we see that :math:`\frac{I}{\sigma_I}` is fairly
  flat over the whole dataset, which we might use as an indication that there
  were no bad frames, not much radiation damage occurred and that scale factors
  are likely to be fairly uniform.

  .. image:: /figures/ioversigma_vs_z.png


Exporting as MTZ
^^^^^^^^^^^^^^^^

The final step of dials processing is to export the integrated results to mtz
format, suitable for input to downstream processing programs such as pointless_
and aimless_.

::

  dials.export integrated.pickle refined_experiments.json mtz.hklout=integrated.mtz

And this is the output, showing the reflection file statistics.

::

  The following parameters have been modified:

  mtz {
    hklout = "integrated.mtz"
  }
  input {
    experiments = refined_experiments.json
    reflections = integrated.pickle
  }

  Removing 23949 reflections with negative variance
  Removing 27432 profile reflections with negative variance
  Removing 2 reflections with I/Sig(I) < -5.0
  Removing 0 profile reflections with I/Sig(I) < -5.0
  Removing 4039 incomplete reflections
  Title: from dials.export
  Space group symbol from file: P4
  Space group number from file: 75
  Space group from matrices: P 4 (No. 75)
  Point group symbol from file: 4
  Number of batches: 540
  Number of crystals: 1
  Number of Miller indices: 323359
  Resolution range: 150.012 1.17004
  History:
  Crystal 1:
    Name: XTAL
    Project: DIALS
    Id: 1
    Unit cell: (57.7873, 57.7873, 150.012, 90, 90, 90)
    Number of datasets: 1
    Dataset 1:
      Name: FROMDIALS
      Id: 1
      Wavelength: 0.97625
      Number of columns: 15
      label        #valid  %valid    min     max type
      H            323359 100.00%   0.00   46.00 H: index h,k,l
      K            323359 100.00%   0.00   47.00 H: index h,k,l
      L            323359 100.00%   0.00  114.00 H: index h,k,l
      M_ISYM       323359 100.00%   1.00    8.00 Y: M/ISYM, packed partial/reject flag and symmetry number
      BATCH        323359 100.00%   2.00  539.00 B: BATCH number
      IPR          323359 100.00%  -2.22 3937.64 J: intensity
      SIGIPR       323359 100.00%   0.05   62.78 Q: standard deviation
      I            323359 100.00% -32.84 4228.04 J: intensity
      SIGI         323359 100.00%   0.11   65.18 Q: standard deviation
      FRACTIONCALC 323359 100.00%   1.00    1.00 R: real
      XDET         323359 100.00%   6.60 2456.31 R: real
      YDET         323359 100.00%   5.79 2520.56 R: real
      ROT          323359 100.00%  82.01  162.69 R: real
      LP           323359 100.00%   0.00    0.76 R: real
      DQE          323359 100.00%   0.71    0.86 R: real



What to do Next
---------------

The following demonstrates how to take the output of dials processing and
continue with downstream analysis using the CCP4 programs pointless_, to sort the data and assign
the correct symmetry, followed by scaling with aimless_ and intensity analysis
using ctruncate_::

  pointless hklin integrated.mtz hklout sorted.mtz > pointless.log
  aimless hklin sorted.mtz hklout scaled.mtz > aimless.log << eof
  resolution 1.3
  anomalous off
  eof
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

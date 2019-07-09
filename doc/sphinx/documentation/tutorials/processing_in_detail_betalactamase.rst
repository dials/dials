Processing in Detail
====================

.. highlight:: none

Introduction
------------

DIALS processing may be performed by either running the individual tools (spot
finding, indexing, refinement, integration, symmetry, scaling, exporting to MTZ)
or you can run :samp:`xia2 pipeline=dials`, which makes informed choices for you
at each stage. In this tutorial we will run through each of the steps in turn,
checking the output as we go. We will also enforce the correct lattice symmetry.


Tutorial data
-------------

The following example uses a Beta-Lactamase dataset collected using
beamline I04 at Diamond Light Source, and reprocessed especially for
these tutorials.

..  hint::
    If you are physically at Diamond on the CCP4 Workshop, then
    this data is already available in your training data area. After
    typing :samp:`module load ccp4-workshop` you'll be moved to a working
    folder, with the data already located in the :samp:`tutorial-data/summed`
    subdirectory.

The data is otherwise available for download from |lactamase|.
We'll only be using the first run of data in this tutorial,
:samp:`C2sum_1.tar`, extracted to a :samp:`tutorial-data/summed` subdirectory.

.. |lactamase|  image::  https://zenodo.org/badge/DOI/10.5281/zenodo.1014387.svg
                :target: https://doi.org/10.5281/zenodo.1014387

Import
^^^^^^

The first stage of step-by-step DIALS processing is to import the data - all
that happens here is that metadata are read for all the images, and a file
describing their contents (:ref:`datablock.expt <datablock-json>`) is written::

    dials.import tutorial-data/summed/C2sum_1*.cbf.gz

The output just describes what the software understands of the images it was
passed, in this case one sweep of data containing 720 images:

.. literalinclude:: logs_detail_betalactamase/dials.import.log

Now is a good point to take a first look at the data using the
:doc:`dials.image_viewer<../programs/dials_image_viewer>`, both to check that
the data is sensible and to anticipate any problems in processing::

  dials.image_viewer datablock.expt

You will be presented with the main image viewer screen:

.. image:: /figures/process_detail_betalactamase/image_viewer.jpg
   :width: 100%

Play with the brightness slider (①) a little until you can clearly see
the spots on the first image (something in the range 10-20 should make
the spots obvious). You can also change the colour scheme (sometimes
spots can be easier to identify in 'inverted' mode) , toggle
various information markers like beam center, and try different
configurations for the spot finding (②).

Find Spots
^^^^^^^^^^

The first "real" task in any processing using DIALS is the spot finding.
Since this is looking for spots on every image in the dataset, this process
can take some time, so we request multiple processors (:samp:`nproc=4`) to
speed this up:

.. literalinclude:: logs_detail_betalactamase/dials.find_spots.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. literalinclude:: logs_detail_betalactamase/dials.find_spots.log
        :linenos:

Once this has completed, a new :ref:`reflection file <reflection_pickle>`
'``strong.refl``' is written, containing a record of every spot found.

The :doc:`dials.image_viewer<../programs/dials_image_viewer>` tool is
not as fast as viewers such as ADXV, however it does integrate well with
DIALS data files. Having found strong spots open the image viewer again,
but giving it the newly found reflection list::

  dials.image_viewer datablock.expt strong.refl

Adjust the brightness so that you can see the spots, then zoom in so
that you can see the clustered individual pixels of a single spot.
Pixels determined to be part of a spot's peak are marked with green
dots. The blue outline shows the three-dimensional **shoebox** - the
extents over detector *x*, *y* and image number *z* of a all peak pixels
in a single spot. The single highest value pixel for any spot is marked
with a pink circle, and the centre of mass is marked with a red cross.

The spot centre-of-mass is usually close to the peak pixel, but slightly
offset as the algorithm allows calculation of the spot centre at a
better precision than the pixel size and image angular 'width'.

.. image:: /figures/process_detail_betalactamase/image_viewer_spot.png

The default parameters for spot finding usually do a good job for
Pilatus images, such as these. However they may not be optimal for data
from other detector types, such as CCDs or image plates. Issues with
incorrectly set gain might, for example, lead to background noise being
extracted as spots. You can use the image mode buttons (③) to preview
how the parameters affect the spot finding algorithm. The final image,
‘threshold’ is the one on which spots were found, so ensuring this produces
peaks at real diffraction spot positions will give the best chance of success.

Another very powerful tool for investigating problems with strong spot positions
is :doc:`dials.reciprocal_lattice_viewer<../programs/dials_reciprocal_lattice_viewer>`.
This displays the strong spots in 3D, after mapping them from their detector
positions to reciprocal space. In a favourable case you should be
able to see the crystal's reciprocal lattice by eye in the strong spot
positions. Some practice may be needed in rotating the lattice to an
orientation that shows off the periodicity in reciprocal lattice positions::

  dials.reciprocal_lattice_viewer datablock.expt strong.refl

.. image:: /figures/process_detail_betalactamase/reciprocal_lattice_strong.png

Although the reciprocal spacing is visible, in this data, there are clearly
some systematic distortions. These will be solved in the indexing.

Indexing
^^^^^^^^

The next step will be indexing of the strong spots by
:doc:`dials.index<../programs/dials_index>`, which by default uses a
3D FFT algorithm (although the 1D FFT algorithm can be selected, using the
parameter :samp:`indexing.method=fft1d`). We pass in all the strong
spots found in the dataset:

.. literalinclude:: logs_detail_betalactamase/dials.index.cmd

If known, the space group and unit cell can be provided at this stage
using the :samp:`space_group` and :samp:`unit_cell` parameters, and will
be used to constrain the lattice during refinement, but otherwise
indexing and refinement will be carried out in the primitive lattice
using space group P1.

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    ..  literalinclude:: logs_detail_betalactamase/dials.index.log
        :linenos:

If successful, ``dials.index`` writes two output data files - an
:ref:`indexed.expt <experiments_json>` containing the tuned
experimental model and determined parameters, and a ``indexed.refl``
reflection file, including index data from the best fit.

It is worth reading through this output to understand what the indexing
program has done. Note that this log is automatically captured in the file
:file:`dials.index.log`, along with a somewhat more detailed log written
into :file:`dials.index.debug.log` - but this second log is probably only
helpful if something has gone wrong and you are trying to track down why.

Inspecting the beginning of the log shows that the indexing step is done
at a resolution lower than the full dataset; 1.84 Å:

.. literalinclude:: logs_detail_betalactamase/dials.index.log
    :start-at: Found max_cell
    :lines: 1-3
    :lineno-match:
    :linenos:

The resolution limit of data that can be used in indexing is determined
by the size of the 3D FFT grid, and the likely maximum cell dimension.
Here we used the default 256³ grid points. These are used to make
an initial estimate for the unit cell parameters.

What then follows are 'macro-cycles' of refinement where the experimental model
is first tuned to get the best possible fit from the data, and then the
resolution limit is reduced to cover more data than the previous cycle.  16
parameters of the diffraction geometry are tuned - 6 for the detector, one for
beam angle, 3 crystal orientation angles and the 6 triclinic cell parameters.
At each stage only 36000 reflections are used in the refinement job. In order
to save time, a subset of the input reflections are used - by default using 100
reflections for every degree of the 360° scan.

We see that the first macrocycle of refinement makes a big improvement in
the positional RMSDs:

.. literalinclude:: logs_detail_betalactamase/dials.index.log
   :start-after: Refinement steps
   :end-before: RMSD no longer decreasing
   :lineno-match:
   :linenos:

Second and subsequent macrocycles are refined using the same number of
reflections, but after extending to higher resolution. The RMSDs at the
start of each cycle start off worse than at the end of the previous
cycle, because the best fit model for lower resolution data is being
applied to higher resolution reflections. As long as each macrocyle
shows a reduction in RMSDs then refinement is doing its job of extending
the applicability of the model out to a new resolution limit, until
eventually the highest resolution strong spots have been included. The
final macrocycle includes data out to 1.30 Å and produces a final model
with RMSDs of 0.050 mm in X, 0.049 mm in Y and 0.104° in φ,
corresponding to 0.29 pixels in X, 0.28 pixels in Y and 0.21 image
widths in φ.

Despite the high quality of this data, we notice from the log that at each
macrocycle there were some outliers identified and removed from
refinement as resolution increases. Large outliers can dominate refinement
using a least squares target, so it is important to be able to remove these.
More about this is discussed below in :ref:`detailbetal-sec-refinement`.
It's also worth checking the total number of reflections that were unable to
be assigned an index:

.. literalinclude:: logs_detail_betalactamase/dials.index.log.extract_unindexed
   :start-after: [START_EXTRACT]
   :end-before:  [END_EXTRACT]
   :lineno-match:
   :linenos:

because this can be an indication of poor data quality or a sign that more
care needs to be taken in selecting the strategy used by ``dials.index``.

After indexing it can be useful to inspect the reciprocal lattice again::

  dials.reciprocal_lattice_viewer indexed.expt indexed.refl

Now indexed/unindexed spots are differentiated by colour, and it is possible
to see which spots were marked by :doc:`dials.refine <../programs/dials_refine>`
as outliers. If you have a dataset with multiple lattices present, it may be
possible to spot them in the unindexed reflections.

In this case, we can see that the refinement has clearly resolved whatever
systematic error was causing distortions in the reciprocal space view, and the
determined reciprocal unit cell fits the data well:

.. image:: /figures/process_detail_betalactamase/reciprocal_lattice_indexed.png


Bravais Lattice Refinement
^^^^^^^^^^^^^^^^^^^^^^^^^^

Since we didn't know the Bravais lattice before indexing, we can now use
:doc:`dials.refine_bravais_settings<../programs/dials_refine_bravais_settings>`
to determine likely candidates. This takes the results of the P1
autoindexing and runs refinement with all of the possible Bravais
settings applied, allowing you to choose your preferred solution:

.. literalinclude:: logs_detail_betalactamase/dials.refine_bravais_settings.cmd

giving a table containing scoring data and unit cell for each Bravais
setting:

.. literalinclude:: logs_detail_betalactamase/dials.refine_bravais_settings.log
    :start-at: Chiral space groups
    :end-before: usr+sys


The scores include the metric fit (in degrees), RMSDs (in mm), and the
best and worse correlation coefficients for data related by symmetry
elements implied by the lowest symmetry space group from the Bravais
setting. This uses the raw spot intensity measurement from the spot-
finding procedure (uncorrected and unscaled) but provides a very useful
check to see if the data does appear to adhere to the proposed symmetry
operators.

A separate ``bravais_setting_N.expt`` experiments file is written for
each plausible lattice type, corresponding to the solution index. In this
example we choose to continue processing with
:samp:`bravais_setting_2.expt`, which is the highest symmetry suggested
result - the options 3, 4, 5 have higher symmetries, but at the cost of
a steep jump in RMSd's and worsening of fit.

In cases where the change of basis operator to the chosen setting is the
identity operator (:samp:`a,b,c`) we can proceed directly to further
refinement. However, we notice that the change of basis operator for our
chosen solution is :samp:`a+b,-a+b,c`, so it is necessary to reindex the
:ref:`indexed.refl <reflection_pickle>` file output by using
:doc:`dials.reindex<../programs/dials_reindex>`:

.. literalinclude:: logs_detail_betalactamase/dials.reindex.cmd

This outputs the file :file:`reindexed.refl` which we now
use as input to downstream programs, in place of the original
:file:`indexed.refl`.

.. _detailbetal-sec-refinement:

Refinement
^^^^^^^^^^

The model is already refined during indexing, but we can also add explicit
refinement steps using :doc:`dials.refine <../programs/dials_refine>`
in here, to use all reflections in refinement rather than a subset and to
fit a scan-varying model of the crystal. There are many options to
refinement - to show all the options up to and including ``expert_level=1``
use this command::

  dials.refine -c -e 1

and descriptions of each of the options can be included by adding ``-a1`` to
the command. All of the main DIALS tools have equivalent command-line options
to list available options.

To refine a static model including the monoclinic constraints
from ``dials.refine_bravais_settings`` run:

.. literalinclude:: logs_detail_betalactamase/dials.refine.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. literalinclude:: logs_detail_betalactamase/dials.refine.log
        :linenos:


This uses all reflections in refinement rather than a subset and provided a
small reduction in RMSDs, writing the results out to ``refined.expt``
and ``refined.refl``.

However, the refined model is still static over
the whole dataset. We may want to do an additional refinement job to fit a
more sophisticated model for the crystal, allowing small misset rotations to
occur over the course of the scan. There are usually even small changes to
the cell dimensions (typically resulting in a net increase in cell volume)
caused by exposure to radiation during data collection. To account for both
of these effects we can extend our parameterisation to obtain a smoothed
*scan-varying* model for both the crystal orientation and unit cell. This means
running a further refinement job starting from the output of the
previous job:

.. literalinclude:: logs_detail_betalactamase/dials.sv_refine.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. literalinclude:: logs_detail_betalactamase/dials.sv_refine.log
        :linenos:

which writes over the ``refined.expt`` and
``refined.refl`` from the previous refinement step. By default the
scan-varying refinement looks for smooth changes over an interval of 36°
intervals, to avoid fitting unphysical models to noise, though this
parameter can be tuned. We can use the :ref:`betalactamase-html-report`,
described shortly, to
view the results of fitting to smoothly varying crystal cell parameters:

.. image:: /figures/process_detail_betalactamase/scan_varying.png

In this tutorial, we see no overall increase in all three cell parameters. If
significant cell volume increases had been observed that might be indicative of
radiation damage. However we can't yet conclude that there is *no* radiation
damage from the *lack* of considerable change observed.


Integration
^^^^^^^^^^^

After the refinement is done the next step is integration, which is performed
by the program :doc:`dials.integrate <../programs/dials_integrate>`. Mostly,
the default parameters are fine for Pilatus data, which will perform
XDS-like 3D profile fitting while using a generalized linear model in order
to fit a Poisson-distributed background model. We will also increase the
number of processors used to speed the job up.

.. literalinclude:: logs_detail_betalactamase/dials.integrate.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. literalinclude:: logs_detail_betalactamase/dials.integrate.log
        :linenos:

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

Symmetry and Scaling
^^^^^^^^^^^^^^^^^^^^

At this point, we have the option to continue processing with dials or to
export the integrated data to do the symmetry and scaling steps with the
CCP4 programs pointless_ and aimless_. For instructions on how to export the
data for further processing, see `Exporting as MTZ`_. In this tutorial we shall
continue to process with dials.

Checking the symmetry
^^^^^^^^^^^^^^^^^^^^^

After integration we can return to our hypothesis of the space group of the
crystal. Although we made an assessment of that when we chose a Bravais lattice
after indexing, we now have better, background-subtracted, values for the
intensities, and for all reflections, not just the strong spots. So, it is
prudent to repeat the assessment to see if there is any indication that our
initial assessment should be revised.::

  dials.symmetry integrated.expt integrated.refl

The symmetry analysis scores all possible symmetry operations by looking at
the intensities of reflections that would be equivalent under that operation.
Then the symmetry operations are combined to score potential space groups::

  Scoring all possible sub-groups
  ---------------------------------------------------------------------------------------------
  Patterson group       Likelihood  NetZcc  Zcc+   Zcc-   CC     CC-    delta  Reindex operator
  ---------------------------------------------------------------------------------------------
  C 1 2/m 1        ***  0.913        9.76    9.76   0.00   0.98   0.00  0.0    -a,b,-c
  P -1                  0.087        0.09    9.81   9.72   0.98   0.97  0.0    -x-y,-x+y,-z
  ---------------------------------------------------------------------------------------------
  Best solution: C 1 2/m 1

Here we see clearly that the best solution is given by :samp:`C 1 2/m 1`, with
a high likelihood, in agreement with the result from
:samp:`dials.refine_bravais_settings`. As we remain confident with this choice,
we now continue to scaling.

Scaling
^^^^^^^

Before the data can be reduced for structure solution, the intensity values must be corrected for
experimental effects which occur prior to the reflection being measured on the
detector. These primarily include sample illumination/absorption effects
and radiation damage, which result in symmetry-equivalent reflections having
unequal measured intensities (i.e. a systematic effect in addition to any
variance due to counting statistics). Thus the purpose of scaling is to determine
a scale factor to apply to each reflection, such that the scaled intensities are
representative of the 'true' scattering intensity from the contents of the unit
cell.

During scaling, a scaling model is created, from which scale factors are calculated
for each reflection. By default, three components are used to create a physical model
for scaling, in a similar manner to that used in the program aimless_.
This model consists of a smoothly varying scale factor as a
function of rotation angle, a smoothly varying B-factor to
account for radiation damage as a function of rotation angle
and an absorption surface correction, dependent on the direction of the incoming
and scattered beam vector relative to the crystal. In this example, we shall
scale the dataset using the output of dials.symmetry with a resolution cutoff of
1.4 Angstrom::

  dials.scale symmetrized.expt symmetrized.refl d_min=1.4

As can be seen from the output text, 70 parameters are used to parameterise the
scaling model for this dataset. Outlier rejection is performed at several stages,
as outliers have a disproportionately large effect during scaling and can lead
to poor scaling results. During scaling, the distribution of the intensity
uncertainties are also analysed and an error model is optimised to transform the
intensity errors to an expected normal distribution. At the end of the output,
a table and summary of the merging statistics are presented, which give indications
of the quality of the scaled dataset::

             ----------Overall merging statistics (non-anomalous)----------

  Resolution: 35.33 - 1.40
  Observations: 276800
  Unique reflections: 41135
  Redundancy: 6.7
  Completeness: 94.11%
  Mean intensity: 82.6
  Mean I/sigma(I): 13.8
  R-merge: 0.065
  R-meas:  0.071
  R-pim:   0.027

The merging statistics, as well as additional output plots, are output into
a html report called :samp:`scaling.html`. This can be opened in your browser -
nativigate to the section "scaling model plots" and take a look.

What is immediately apparent is the periodic nature of the scale term, with peaks
and troughs 90° apart. This indicates that the illumated volume was changing
significantly during the experiment: a reflection would be measured as almost
twice as intense if it was measured at rotation angle of ~120° compared to at ~210°.
The absorption surface also shows a similar periodicity, as may be expected.
The form of the relative B-factor is less well defined, although is does show some
periodicity similar to the scale term. There is certainly no strong B-factor
reduction as a function of rotation angle, which would have suggested radiation
damage. The scaling can be repeated, omitting the :samp:`decay_term`::

  dials.scale symmetrized.expt symmetrized.refl d_min=1.4 decay_term=False

::

             ----------Overall merging statistics (non-anomalous)----------

  Resolution: 35.33 - 1.40
  Observations: 276792
  Unique reflections: 41135
  Redundancy: 6.7
  Completeness: 94.11%
  Mean intensity: 75.9
  Mean I/sigma(I): 14.3
  R-merge: 0.064
  R-meas:  0.069
  R-pim:   0.026


By inspecting the statistics in the output, we can see that removing the decay
term has slightly improved some of the R-factors and mean I/sigma(I). Therefore
it is probably best to exclude the decay correction for this dataset.

.. _betalactamase-html-report:

HTML report
^^^^^^^^^^^

Much more information from the various steps of data processing can be found
within an HTML report generated using the program
:doc:`dials.report <../programs/dials_report>`.
This is run simply with::

  dials.report scaled.expt scaled.refl

which produces the file :download:`dials-report.html <logs_detail_betalactamase/dials-report.html>`.

This report includes plots showing the scan-varying crystal orientation
and unit cell parameters. The latter of these is useful to check that
changes to the cell during processing appear reasonable. We can at least
see from this and the low final refined RMSDs that this is a very well-
behaved dataset.

Some of the most useful plots are

* **Difference between observed and calculated centroids vs phi**,
  which shows how the average
  residuals in each of X, Y, and φ vary as a fuction of φ.
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

The final step of dials processing is to either 1) export the integrated results to mtz
format, suitable for input to downstream processing programs such as pointless_
and aimless_.

.. literalinclude:: logs_detail_betalactamase/dials.export.cmd

2) export the scaled intensities for further downstream processing, making sure to
include the :samp:`intensity=scale` option::

  dials.export scaled.refl scaled.expt intensity=scale

Here is the output for exporting after integration, showing the reflection file statistics.

.. literalinclude:: logs_detail_betalactamase/dials.export.log
    :linenos:

Alternative processing with pointless and aimless
-------------------------------------------------

The following demonstrates how to take the output of dials processing after integration and
continue with downstream analysis using the CCP4 programs pointless_, to sort the data and assign
the correct symmetry, followed by scaling with aimless_ and intensity analysis
using ctruncate_::

  pointless hklin integrated.mtz hklout sorted.mtz > pointless.log
  aimless hklin sorted.mtz hklout scaled.mtz > aimless.log << EOF
  resolution 1.4
  EOF
  ctruncate -hklin scaled.mtz -hklout truncated.mtz \
  -colin '/*/*/[IMEAN,SIGIMEAN]' > ctruncate.log

to get merged data for downstream analysis. The output from this includes
the merging statistics which will give a better idea about data quality. It is
easiest to view these logfiles using the program :program:`logview`, e.g.::

  logview aimless.log

Often passing in a sensible resolution limit to aimless is helpful. Here
we assumed we ran first without a resolution limit to help decide where
to cut the data. Following this we chose to exclude all data at a
resolution higher than 1.4 Angstroms, to ensure about 90% completeness
in the outer shell. Here is the summary from aimless.log:

::

    Summary data for        Project: DIALS Crystal: XTAL Dataset: FROMDIALS

                                              Overall  InnerShell  OuterShell
    Low resolution limit                       69.19     69.19      1.42
    High resolution limit                       1.40      7.54      1.40

    Rmerge  (within I+/I-)                     0.056     0.028     0.598
    Rmerge  (all I+ and I-)                    0.066     0.039     0.670
    Rmeas (within I+/I-)                       0.067     0.033     0.733
    Rmeas (all I+ & I-)                        0.072     0.043     0.737
    Rpim (within I+/I-)                        0.036     0.017     0.417
    Rpim (all I+ & I-)                         0.027     0.017     0.302
    Rmerge in top intensity bin                0.029        -         -
    Total number of observations              276017      2016     11442
    Total number unique                        41113       300      1980
    Mean((I)/sd(I))                             14.4      47.7       2.1
    Mn(I) half-set correlation CC(1/2)         0.999     0.998     0.808
    Completeness                                94.3      99.3      90.5
    Multiplicity                                 6.7       6.7       5.8
    Mean(Chi^2)                                 0.92      0.77      0.83

    Anomalous completeness                      94.3     100.0      89.0
    Anomalous multiplicity                       3.4       3.7       2.9
    DelAnom correlation between half-sets      0.363     0.683     0.033
    Mid-Slope of Anom Normal Probability       1.140       -         -



.. _pointless: http://www.ccp4.ac.uk/html/pointless.html
.. _aimless: http://www.ccp4.ac.uk/html/aimless.html
.. _ctruncate: http://www.ccp4.ac.uk/html/ctruncate.html

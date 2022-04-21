.. raw:: html

  <a href="https://dials.github.io/dials-2.2/documentation/tutorials/processing_in_detail_betalactamase.html" class="new-documentation">
  This tutorial requires a DIALS 3 installation.<br/>
  Please click here to go to the tutorial for DIALS 2.2.
  </a>

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
describing their contents (:ref:`imported.expt <experiments_json>`) is written::

    dials.import tutorial-data/summed/C2sum_1*.cbf.gz

The output just describes what the software understands of the images it was
passed, in this case one sequence of data containing 720 images:

.. dials_tutorial_include:: betalactamase/dials.import.log

Now is a good point to take a first look at the data using the
:doc:`dials.image_viewer<../programs/dials_image_viewer>`, both to check that
the data is sensible and to anticipate any problems in processing::

  dials.image_viewer imported.expt

You will be presented with the main image viewer screen:

.. image:: https://dials.github.io/images/process_detail_betalactamase/image_viewer.png
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
can take some time, so DIALS will use multiple processors by default to
speed this up. Here we have limited it to 4, but feel free to omit this to
let DIALS make the choice:

.. dials_tutorial_include:: betalactamase/dials.find_spots.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. dials_tutorial_include:: betalactamase/dials.find_spots.log
        :linenos:

Once this has completed, a new :ref:`reflection file <reflection_pickle>`
'``strong.refl``' is written, containing a record of every spot found.

The :doc:`dials.image_viewer<../programs/dials_image_viewer>` tool is
not as fast as viewers such as ADXV, however it does integrate well with
DIALS data files. Having found strong spots open the image viewer again,
but giving it the newly found reflection list::

  dials.image_viewer imported.expt strong.refl

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

.. image:: https://dials.github.io/images/process_detail_betalactamase/image_viewer_spot.png

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

  dials.reciprocal_lattice_viewer imported.expt strong.refl

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

.. dials_tutorial_include:: betalactamase/dials.index.cmd

If known, the space group and unit cell can be provided at this stage
using the :samp:`space_group` and :samp:`unit_cell` parameters, and will
be used to constrain the lattice during refinement, but otherwise
indexing and refinement will be carried out in the primitive lattice
using space group P1.

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    ..  dials_tutorial_include:: betalactamase/dials.index.log
        :linenos:

If successful, ``dials.index`` writes two output data files - an
``indexed.expt`` containing the tuned
experimental model and determined parameters, and a ``indexed.refl``
reflection file, including index data from the best fit.

It is worth reading through this output to understand what the indexing
program has done. Note that this log is automatically captured in the file
:file:`dials.index.log`. A more verbose debug log can be generated by adding
the '-v' option to a dials command line program, but this is probably only
helpful if something has gone wrong and you are trying to track down why.

Inspecting the beginning of the log shows that the indexing step is done
at a resolution lower than the full dataset; 1.84 Å:

.. dials_tutorial_include:: betalactamase/dials.index.log
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

.. dials_tutorial_include:: betalactamase/dials.index.log
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

.. dials_tutorial_include:: betalactamase/dials.index.log.extract_unindexed
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

.. dials_tutorial_include:: betalactamase/dials.refine_bravais_settings.cmd

giving a table containing scoring data and unit cell for each Bravais
setting:

.. dials_tutorial_include:: betalactamase/dials.refine_bravais_settings.log
    :start-at: Chiral space groups

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

.. dials_tutorial_include:: betalactamase/dials.reindex.cmd

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

To refine over all reflections, and include the monoclinic constraints
from ``dials.refine_bravais_settings`` run:

.. dials_tutorial_include:: betalactamase/dials.refine.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. dials_tutorial_include:: betalactamase/dials.refine.log
        :linenos:

This provides a good reduction in RMSDs, indicating a better fit, and writes
the results out to ``refined.expt`` and ``refined.refl``.

Two passes of refinement are actually done here - an initial pass where unit
cell and crystal rotation is consistent over the length of the experiment, and
a second pass where these are allowed to vary. This *scan-varying* refinement
allows compensation for small missets in the rotation of the goniometer, and
compensation for changes to the unit cell dimensions; typically due to
radiation damage. By default, the refinement looks for smooth changes over
intervals of 30°, to avoid fitting unphysical models to noise, though this
interval can be configured.

We can use the :ref:`betalactamase-html-report`, described shortly, to view the
results of fitting to smoothly varying crystal cell parameters:

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

.. dials_tutorial_include:: betalactamase/dials.integrate.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. dials_tutorial_include:: betalactamase/dials.integrate.log
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


Symmetry analysis
^^^^^^^^^^^^^^^^^

After integration, further assessments of the crystal symmetry are possible.
Previously, we made an assessment of the lattice symmetry (i.e. the symmetry
of the diffraction spot positions), however now we have determined a set of
intensity values and can investigate the full symmetry of the diffraction
pattern (i.e. spot positions and intensities). The symmetry analysis consists
of two stages, determining the laue group symmetry and analysing absent
reflections to suggest the space group symmetry.

.. dials_tutorial_include:: betalactamase/dials.symmetry.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. dials_tutorial_include:: betalactamase/dials.symmetry.log
        :linenos:

The laue group symmetry is the 3D rotational symmetry of the diffraction
pattern plus inversion symmetry (due to Friedel's law that I(h,k,l) = I(-h,-k,-l)
when absorption is negligible). To determine the laue group symmetry, all
possible symmetry operations of the lattice are scored by comparing the
correlation of reflection intensities that would be equivalent under a given
operation. The scores for individual symmetry operations are then combined to
score the potential laue groups.

.. dials_tutorial_include:: betalactamase/dials.symmetry.log
    :start-at: Scoring all possible sub-groups
    :end-before: Analysing systematic absences

Here we see clearly that the best solution is given by C 1 2/m 1, with
a high likelihood. For macromolecules, their chirality means that mirror symmetry
is not allowed (the 'm' in C 1 2/m 1), therefore the determined symmetry
relevant for MX at this point is C2. For some laue groups, there are multiple
space groups possible due additional translational symmetries
(e.g P 2, P 2\ :sub:`1` for laue group P2/m), which requires an additional
analysis of systematic absences. However this is not the case for C 1 2/m 1,
therefore the final result of the analysis is the space group C2, in agreement
with the result from :samp:`dials.refine_bravais_settings`.

Scaling and Merging
^^^^^^^^^^^^^^^^^^^

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
for each reflection. Three physically motivated corrections are used to create an
scaling model, in a similar manner to that used in the program aimless_.
This model consists of a smoothly varying scale factor as a
function of rotation angle, a smoothly varying B-factor to
account for radiation damage as a function of rotation angle
and an absorption surface correction, dependent on the direction of the incoming
and scattered beam vector relative to the crystal.

.. dials_tutorial_include:: betalactamase/dials.scale.cmd

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. dials_tutorial_include:: betalactamase/dials.scale.log
        :linenos:

As can be seen from the output text, 70 parameters are used to parameterise the
scaling model for this dataset. Outlier rejection is performed at several stages,
as outliers have a disproportionately large effect during scaling and can lead
to poor scaling results. During scaling, the distribution of the intensity
uncertainties are also analysed and a correction is applied based on a prior
expectation of the intensity error distribution. At the end of the output,
a table and summary of the merging statistics are presented, which give indications
of the quality of the scaled dataset:

.. dials_tutorial_include:: betalactamase/dials.scale.log
    :start-at: ----------Merging statistics by resolution bin----------
    :end-before: Writing html report to dials.scale.html

Looking at the resolution-dependent merging statistics, we can see that the
completeness falls significantly beyond 1.4 Angstrom resolution.
If desired, a resolution cutoff can be applied and the
data rescaled (using the output of the previous scaling run as input to the
next run to load the existing state of the scaling model):

.. dials_tutorial_include:: betalactamase/dials.scale_cut.cmd

The merging statistics, as well as a number of scaling and merging plots, are
output into a html report called :samp:`dials.scale.html`.
This can be opened in your browser - nativigate to the section "scaling model plots" and take a look.
What is immediately apparent is the periodic nature of the scale term, with peaks
and troughs 90° apart. This indicates that the illuminated volume was changing
significantly during the experiment: a reflection would be measured as almost
twice as intense if it was measured at rotation angle of ~120° compared to at ~210°.
The absorption surface also shows a similar periodicity, as may be expected.
The relative B-factor shows low overall variation, suggesting little overall
radiation damage.

Once we are happy with the dataset quality, the final step of dials processing
is to merge the data and produce a merged mtz file, suitable for input to
downstream structure solution. To do this we can use the command::

  dials.merge scaled.expt scaled.refl

The log output reports intensity statistics, the symmetry equivalent reflections
are merged and a truncation procedure is performed, to give strictly positive
merged structure factors (Fs) in addition to merged intensities.

.. _betalactamase-html-report:

HTML report
^^^^^^^^^^^

Much more information from the various steps of data processing can be found
within an HTML report generated using the program
:doc:`dials.report <../programs/dials_report>`.
This is run simply with::

  dials.report scaled.expt scaled.refl

which produces the file :download:`dials.report.html <betalactamase-report.html>`.

This report includes plots showing the scan-varying crystal orientation
and unit cell parameters. The latter of these is useful to check that
changes to the cell during processing appear reasonable. We can at least
see from this and the low final refined RMSDs that this is a very well-
behaved dataset.

Some of the most useful plots are

* **Difference between observed and calculated centroids vs phi**,
  which shows how the average
  residuals in each of X, Y, and φ vary as a function of φ.
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

Exporting to unmerged MTZ
^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible that an unmerged mtz file is desired for further processing before
merging. To produce a scaled unmerged mtz file, one can use the ``dials.export``
command on the scaled datafiles::

  dials.export scaled.refl scaled.expt

It is also possible to export the integrated (unscaled) data in mtz
format using :samp:`dials.export`. If you have an installation of CCP4_, symmetry
analysis and scaling can then be continued with the ccp4 programs
pointless_, aimless_ and ctruncate_ to generate a merged mtz file::

  dials.export integrated.refl integrated.expt
  pointless hklin integrated.mtz hklout sorted.mtz > pointless.log
  aimless hklin sorted.mtz hklout scaled.mtz > aimless.log << EOF
  resolution 1.4
  anomalous off
  EOF
  ctruncate -hklin scaled.mtz -hklout truncated.mtz \
  -colin '/*/*/[IMEAN,SIGIMEAN]' > ctruncate.log

.. _CCP4: http://www.ccp4.ac.uk
.. _aimless: http://www.ccp4.ac.uk/html/aimless.html
.. _pointless: http://www.ccp4.ac.uk/html/pointless.html
.. _ctruncate: http://www.ccp4.ac.uk/html/ctruncate.html

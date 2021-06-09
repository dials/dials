Processing in Detail (Thaumatin / Eiger Edition)
================================================

.. highlight:: none

Introduction
------------

DIALS processing may be performed by either running the individual
tools (spot finding, indexing, refinement, integration, symmetry,
scaling, exporting to MTZ) or you can run `xia2`, which makes
informed choices for you at each stage. In this tutorial we will run
through each of the steps in turn, checking the output as we go. We will also enforce the correct lattice symmetry.


Tutorial data
-------------

The following example uses a small Thaumatin data set collected on
beamline i03 at Diamond Light Source, which is available from Zenodo
at https://doi.org/10.5281/zenodo.4916648

Files
^^^^^

DIALS creates two principle file types:

- experiment files called `something.expt`
- reflection files called `something.refl`

In most cases the filenames will correspond to the name of the DIALS
program which created them e.g. `indexed.refl` and `indexed.expt` from
`dials.index`. The only deviations from this are on import (see below)
where we are only reading experiment models and spot finding where we
find _strong_ reflections so write these to `strong.refl` - and we
create no models so (by default) there is no output experiment file. 

At any time you can _look_ at these files with `dials.show` which will
summarise the content of the files to the terminal. 

Parameters
^^^^^^^^^^

All DIALS programs accept parameters in the form of
`parameter=value` - in most cases this will be sufficient though some
less frequently used options may require "name space" clarification
e.g. `index_assignment.method=local`. All of the DIALS programs
support the option

::
   dials.program -c -e2

which will show you all possible configuration options - if you are
looking for an option this is the simplest way to search so e.g.

::
   dials.index -c -e2 | less

will allow you to scroll through the extensive list of options you can
adjust. In most cases the defaults are relatively sensible for
synchrotron data from a pixel array detector, as we are using in this
tutorial. 

Output
^^^^^^

In the majority of cases the `dials` programs write their output to
`dials.program.log` e.g. `dials.find_spots.log` etc. - everything
which is printed to the terminal is also saved in this file, so you
can review the processing later. In the case where you are reporting
an issue to the developers including these log files in the error
report (particularly for the step which failed) is very helpful. 

From most stages you can generate a more verbose _report_ of the
current state of processing with:

::
   dials.report step.expt step.refl

which will generate a detailed report as HTML describing the current
state of the processing. 
   
Import
^^^^^^

The starting point for any processing with DIALS is to _import_ the
data - here the metadata are read and a description of the data to be
processed saved to a file named `imported.expt`. This is "human
readable" in that the file is JSON format (roughly readable text with
brackets around to structure for computers). While you can edit this
file if you know what you are doing, usually this is not necessary. 

::
   dials.import xtal_1_5_master.h5

will read the metadata from this `master` file and write
`imported.expt` from this - equally in this case you could import from
the NeXus formatted file (which is functionally equivalent) with

::
   dials.import xtal_1_5.nxs

It is important to note that for well behaved data (i.e. anything
which is well collected from a well behaved sample) the commands below
will often be identical after importing.

At this point you can actually look at the images with the
`dials.image_viewer` tool - 

::
   dials.image_viewer imported.expt

in this tool there are many settings you can adjust, which could
depend on the source of the data and - most importantly - your
preferences. Personally the author finds for basic inspection of the
images the brightness is a bit high for pixel array data, and a value
of 10 may be better for viewing the diffraction pattern as a whole.

To get a sense of how the diffraction spots are spread, stacking
images can help - for example in this case setting the stack to 10
gives a good idea of the real separation between reflections. If the
data are not stacked the spot finding process can also be explored -
the controls at the bottom of the "Settings" window allow you to step
through these and can be very useful for getting a "computer's eye
view" of how the data look (particularly for establishing where the
diffraction is visible to.)

Find Spots
^^^^^^^^^^

The first "real" task in any processing using DIALS is the spot
finding. Since this is looking for spots on every image in the
dataset, this process can take some time so by default will use all of
the processors available in your machine - if you would like to
control this adjust with e.g. `nproc=4` - however the default is
usually sensible unless you are sharing the computer with many
others.

::
   dials.find_spots imported.expt

This is one of the two steps where every image in the data set is read
and processed and hence can be moderately time consuming. This
contains a reflection file `strong.refl` which contains both the
positions of the strong spots and also "images" of the spot pixels
which we will use later. You can view these spots on top of the images
with

::
   dials.image_viewer imported.expt strong.refl

to get a sense of what spots were found. You will see that the spots
are surrounded by little blue boxes - these are the _bounding boxes_ of
the reflections i.e. the outer extent of the connected regions of the
signal pixels. The signal pixels are highlighted with green blobs
giving a sense of what is and is not "strong."

The default parameters for spot finding usually do a good job for
Pilatus images, such as these. However they may not be optimal for data
from other detector types, such as CCDs or image plates. Issues with
incorrectly set gain might, for example, lead to background noise being
extracted as spots. You can use the image mode buttons to preview
how the parameters affect the spot finding algorithm. The final button
'threshold’ is the one on which spots were found, so ensuring this
produces peaks at real diffraction spot positions will give the best
chance of success. 

The second tool for visualisation of the found spots is the reciprocal
lattice viewer - which presents a view of the spot positions mapped to
reciprocal space.

::
   dials.reciprocal_lattice_viewer imported.expt strong.refl

No matter the sample orientation you should be able
to rotate the space to "look down" the lines of reflections. If you
cannot, or the lines are not straight, it is likely that there are
some errors in the experiment parameters e.g. detector distance or
beam centre. If these are not too large they will likely be corrected
in the subsequent analysis - however you may find it useful to run

::
   dials.search_beam_position imported.expt strong.refl

to determine an updated position for the beam centre - running the
reciprocal lattice viewer with the optimised experiment output:

::
   dials.reciprocal_lattice_viewer optimised.expt strong.refl

should show straight lines, provided everything has worked correctly. 


Indexing
^^^^^^^^

The next step will be indexing of the found spots with `dials.index` -
by default this uses a 3D FFT algorithm to identify periodicy in the
reciprocal space mapped spot positions, though there are other
algorithms available which can be better suited to e.g. narrow data
sets.

::
   dials.index imported.expt strong.refl

or

::
   dials.index optimised.expt strong.refl
   
are the ways to trigger the program, and the most common parameters to
set are the `space_group` and `unit_cell` if these are known in
advance. While this does index the data it will also perform some
refinement with a static crystal model, and indicate in the output the
fraction of reflections which have been indexed - ideally this should
be close to 100%. If it is significantly less than 100% it is possible
you have a second lattice - adding `max_lattices=2` (say) to the
command-line will indicate to the program that you would like to
consider attempting to separately index the unindexed reflections
after the first lattice has been identified. 

By default the triclinic lattice i.e. with `P1` no additional symmetry
is assumed - for the majority of data there are no differences in the
quality of the results from assigning the Bravais lattice at this
stage. 

If successful, `dials.index` writes the experiments and indexed
reflections to two new files `indexed.expt` and `indexed.refl` - if
these are loaded in the reciprocal lattice viewer you can see which
spots have been indexed and if you have multiple lattices switch them
"on and off" for comparison. 

The process that the indexing performs is quite complex -

- make a guess at the maximum unit cell from the pairwise separation
  of spots in reciprocal space
- transform spot positions to reciprocal space using the best
  available current model of the experimental geometry
- perform a Fourier transform of these positions or other algorithm to
  identify the _basis vectors_ of these positions e.g. the spacing
  between one position and the next
- determine a set of these basis vectors which best describes the
  reciprocal space positions
- transform this set of three basis vectors into a unit cell
  description, which is then manipulated according to some standard
  rules to give the best _triclinic_ unit cell to describe the
  reflections - if a unit cell and space group have been provided
  these will be enforced at this stage
- _assign indices_ to the reflections by "dividing through"
  the reciprocal space position by the unit cell parallelopiped (this
  is strictly the actual indexing step)
- take the indexed reflections and refine the unit cell parameters and
  model of the experimental geometry by comparing where the
  reflections should be and where they are found
- save the indexed reflections and experiment models to the output
  files

The indexing process takes place over a number of cycles, where low
resolution reflections are initially indexed and refined before
including more reflections at high resolution - this improves the
overall success of the procedure by allowing some refinement as a part
of the process. 
  
During this process an effort is made to eliminate "outlier"
reflections - these are reflections which do not strictly belong to
the crystal lattice but are accidentally close to a reciprocal space
position and hence can be indexed. Most often this is an issue with
small satellite lattices or ice / powder on the sample. Usually this
should not be a cause for concern. 


Bravais Lattice Determination (optional!)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have indexed the data you may optionally attempt to infer the
correct Bravais lattice and assign this to constrain the unit cell in
subsequent processing. If, for example, the unit cell from indexing
has all three angles close to 90 degrees and two unit cell lengths
with very similar values you could guess that the unit cell is
tetragonal. In `dials.refine_bravais_settings` we take away the
guesswork by transforming the unit cell to all possible Bravais
lattices which approximately match the triclinic unit cell, and then
performing some refinement - if the lattice constraints are correct
then imposing them should have little impact on the deviations between
the observed and calculated reflection positions (known as the R.M.S.
deviations). If a lattice constraint is incorrect it will manifest as
a significant increase in a deviation - however care must be taken as
it can be the case that the true _symmetry_ is lower than the shape of
the unit cell would indicate.

In the general case there is little harm in skipping this step,
however for information if you run

::
   dials.refine_bravais_settings indexed.expt indexed.refl

you will see a table of possible unit cell / Bravais lattice /
R.M.S. deviations printed in the output - in the case of this tutorial
data they will all match, as the true symmetry is tetragonal.

If you wish to use one of the output experiments from this process
e.g. `bravais_setting_9.expt` you will need to reindex the reflection
data from indexing to match this - we do not output every option of
reindexed data as these files can be large. In most cases it is
simpler to re-run `dials.index` setting the chosen space group. 

The reader is reminded here - in most cases it is absolutely fine to
proceed without worrying about the crystal symmetry at this stage `:-)` 


Refinement
^^^^^^^^^^

The model is already refined during indexing, but this is assuming
that a single crystal model is appropriate for every image in the data
set - in reality there are usually small changes in the unit cell and
crystal orientation throughout the experiment as the sample is
rotated. `dials.refine` will first re-run refinement with a fixed unit
cell and then perform scan-varying refinement. If you have indexed
multiple sweeps earlier in processing (not covered in this tutorial)
then the crystal models will be copied and split at this stage to
allow per-crystal-per-scan models to be refined. 

By and large one may run:

::
   dials.refine indexed.expt indexed.refl

without any options and the program will do something sensible - if
you compare the R.M.S. deviations from the end of indexing with the
end of refinement you should see a small improvement. If you look at
the output of `dials.report` at this stage you should see small
variations in the unit cell and sample orientation as the crystal is
rotated - if these do not appear small then it is likely that
something has happened during data collection e.g. severe radiation
damage. 
   

Integration
^^^^^^^^^^^

Once you have refined the model the next step is to integrate the
data - in effect this is using the refined model to calculate the
positions where all of the reflections in the data set will be found
and measure the background substracted intensities:

::
   dials.integrate refined.expt refined.refl

By default this will pass through the data twice, first looking at the
shapes of the predicted spots to form a reference profile model then
passing through a second time to use this profile model to integrate
the data, by being fit to the transformed pixel values. This is by far
the most computationally expensive step in the processing of the data!

by default all the processors in your computer are used, unless we
think this will exceed the memory available in the machine. At times,
however, if you have a large unit cell and / or a large data set you
may find that processing on a desktop workstation is more appropriate
than e.g. a laptop.

If you know in advance that the data do not diffract to anything close
to the edges of the detector you can assign a resolution limit at this
stage by adding `prediction.d_min=1.8` (say) to define a 1.8A
resolution limit - this should in general not be necessary. At the end
of integration two new files are created - `integrated.refl` and
`integrated.expt` - looking at these in the image viewer e.g.

::
   dials.image_viewer integrated.expt integrated.refl

can be very enlightening as you should see little red boxes around
every reflection - if you select "integrated only" you can see what
was and was not integrated. You may see a selection of reflections
close to the rotation axis are missed - these are not well modelled or
predicted in any program so typically excluded from processing. 


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

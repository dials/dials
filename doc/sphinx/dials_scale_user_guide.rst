User guide for scaling data with DIALS
======================================

This document aims to provide a guide to using :samp:`dials.scale`, at various levels
of depth. A new user is encouraged to read the Symmetry and Scaling sections of
the `Processing in detail
<https://dials.github.io/documentation/tutorials/processing_in_detail_betalactamase.html>`_
tutorial for a quick overview of scaling in DIALS.
For most users, it is likely to be sufficient to read only the
'Guide to common scaling options' below,
and return to the rest of the guide if further help is needed.

As a reminder, this is how to run routine data processing after integration to
obtain a merged MTZ file::

  dials.symmetry integrated.refl integrated.expt
  dials.scale symmetrized.refl symmetrized.expt
  dials.merge scaled.refl scaled.expt

The user is also advised to familiarise themselves with the standard program
output, which may contain useful information, and the html report generated
by scaling, which provides numerous plots relating to the merging statistics.

Guide to common scaling options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
These sections cover the most commonly used options (with example values)
for scaling routine macromolecular crystallography datasets.

**Cutting back data**
After inspecting the statistics in the :samp:`dials.scale.html` file, such as
the R-merge vs batch plot, it is often the case that not all of the data are
suitable for merging, perhaps due to radiation damage or nonisomorphism.
This can be the case within a  single sweep or across multiple sweeps in a
multi-sweep/multi-crystal experiment. These are example options to use:

- :samp:`d_min=2.0`  Applies a resolution cutoff at the given resolution (in Angstrom).
- :samp:`exclude_images="100:120"`  Removes a section of images for a single
  sweep dataset. Multiple commands like this can be used to exclude multiple ranges.
  In the case of multiple-sweeps, one must also provide the experiment ID that
  the exclusion should apply to, with the syntax :samp:`exclude_images="a:b:c"`
  where `a` is the experiment ID (a number starting at :samp:`0`), `b` is the initial
  image to exclude and `c` is the final image to exclude.
- :samp:`exclude_datasets="10 50 79"`  Removes whole datasets, based on
  the dataset number; useful for large multi-crystal datasets.

**Anomalous data**
During scaling, the option :samp:`anomalous=[True|False]` determines whether
anomalous pairs (I+/I-) are combined during scaling model minimisation and outlier
rejection. By default, :samp:`anomalous=False`, which is suitable for data with some anomalous signal, however for strongly anomalous data,
the anomalous signal strength may be enhanced when scaling with :samp:`anomalous=True`

**Controlling the absorption correction**
The default physical scaling model applies a relative absorption correction based
on the incoming and outgoing scattering vectors (this accounts for the relative
difference in absorption for different scattering paths through the crystal,
rather than absolute absorption of the beam by the crystal). This correction is
constrained, and the level of constraint and parameterisation can be changed
with the option :samp:`absorption_level=[low|medium|high]`. This aims to give
relative absorption corrections of around 1%, 5% and 25%, but will depend on
the dataset. To see the extent of the correction, check the 'scaling models'
section in the :samp:`dials.scale.html` file.

**Generating MTZ files**
For convenience, :samp:`dials.scale` can invoke the exporting and merging programs to
generate unmerged and merged MTZ files (you may want to use the individual
programs to have more extensive control over the program options):

- :samp:`merged_mtz=scaled.mtz`  Create a merged MTZ file, using the merging routines
  available in cctbx.
- :samp:`unmerged_mtz=unmerged.mtz`  Output the scaled data in unmerged MTZ format.

**Choosing which integrated intensity to use**
One choice that is made automatically during scaling is whether summation or
profile intensities seem give the best estimate of the integrated intensity
(or a combination of the two). To see the result of this combination, inspect the
table in the scaling log, which scores a set of Imid values on Rpim \& CC1/2.
To specify which intensity choice to use, there are a couple of options:

- :samp:`intensity_choice=[profile|sum|combine]`  Choose from profile, sum or combine (default is combine)
- :samp:`combine.Imid=700.0`  Specify the crossover value for profile-summation
  intensity combination.

**Adjusting the uncertainties/errors**
All scaling programs adjust the uncertainties (sigmas) of the integrated data, to
account for additional systematic errors not suffiently modelled during integration.
:samp:`dials.scale` adjusts the intensity errors by refining a two-component error model
(see the output log or :samp:`dials.scale.html` for the values). While this is
an important correction and should improve the data quality for typical
macromolecular crystallographic data, for poorer quality data the model parameters
may become overinflated.
If so, then this correction can be controlled with the parameters:

- :samp:`error_model=None`  Don't apply an error model.
- :samp:`error_model.basic.minimisation=None`  Don't refine the error model in this
  scaling run. Will keep the pre-existing error model parameters, or the default
  error model (:samp:`a=1.0, b=0.02`) on a first scaling run.

For the multi-sweep case, a single error model is applied to the combined dataset,
on the assumption that a similar systematic error is affecting all sweeps. This
approach may not be optimal for some datasets. As an alternative, a separate error
model can be refined on sweeps individually or as groups.

- :samp:`error_model.grouping=[individual|grouped|combined]`  If grouped is chosen,
  then the groups must be specified as below.
- :samp:`error_model_group='0 1' error_model_group='2 3'` e.g. groups the sweeps
  in pairs for error model refinement.

**Controlling partials**
By default, reflections with a partiality above 0.4 are included in the output
data files and merging statistics from dials.scale. This threshold can be changed
with the parameters:

- :samp:`partiality_threshold=0.95`  Disregard all measurements with partialities
  below this value.


Practicalities for large datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Depending on the computational resources available, scaling of large datasets
( > 1 million reflections) can become slow and memory intensive.
There are several options available for managing this.
The first option is separating the data in memory to allow blockwise calculations
and parallel processing, using the option :samp:`nproc=` (a value of 4 or 8 is
probably a reasonable choice).
One of the most computationally-intensive parts of the algorithm is the final
round of minimisation, which uses full-matrix methods. One can set
:samp:`full_matrix=False` to turn this off, however no errors for the scale
factors will be determined. A compromise is to set
:samp:`full_matrix_max_iterations=1` to do at least one iteration.
A third option is to reduce the number of reflections used by the scaling
algorithm during minimisation. If using :samp:`reflection_selection.method=auto`,
the number of reflections should be manageable even for very large datasets, but
this can always be controlled by the user. To get started, use the command
:samp:`dials.scale -ce2` to see the full set of available options in the section
:samp:`reflection_selection`. Try setting :samp:`reflection_selection.method=quasi_random`
alongside some of the :samp:`quasi_random` parameters.


Scaling against a reference dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
DIALS contains functionality for scaling against a reference dataset, also
referred to as targeted scaling.
This reference can either be a dataset scaled with dials.scale, or an mtz file
containing a scaled dataset. The scaled data (excluding the reference) will
be output in a single .refl/.expt file.

**Scaling against a dials reference dataset.**
In this example, reference.refl and reference.expt are from a dataset that has
already been scaled with dials.scale. To scale another dataset (datafiles
:samp:`integrated.refl integrated.expt`) against this reference, one should use the
following command::

  dials.scale only_target=True integrated.refl integrated.expt reference.refl reference.expt

This will scale the intensities of the dataset to agree as closely as possible
with the intensities of the reference dataset. The :samp:`only_target=True`
command is important, else all the data will be scaled together and output in
a joint output file.

**Scaling against a reference mtz file.**
In this case, it is assumed that the intensity and variance columns of the mtz
file have already been scaled. Reference scaling would be run with the following
command::

  dials.scale integrated.refl integrated.expt target_mtz=scaled.mtz

The reference scaling algorithm is the same regardless of the target datafile type.


Advanced use - Controlling the scaling models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are three available scaling models available in dials.scale, accessible
by the command line option :samp:`model = physical array KB *auto`.
The physical model is similar to the scaling model used in the program aimless_,
the array model is based on the approach taken in xscale_, while the KB model is
a simple two-component model suitable for still-image datasets or very small
rotation datasets (~ < 1 degree).

The auto option automatically chooses a default model and sensible parameterisation
based on the oscillation range of the experiment. This will choose the
physical model unless the oscillation range is < 1.0 degree, when the KB model
will be chosen. If the oscillation range is < 60 degrees, the absorption correction
of the physical model is disabled, as this may be poorly determined. The parameter
spacing as a function of rotation is also adjusted down from the defaults if the
oscillation range is below 90 degrees, to try to give a sensible automatic
parameterisation.

The physical model consists of up to three components; a smoothly varying
scale correction, a smoothly varying B-factor correction and an absorption surface
correction (all on by default). These are turned on/off with the command line options
:samp:`physical.scale_correction=True/False physical.decay_correction=True/False physical.absorption_correction=True/False`.
The smoothly varying terms have a parameter at regular intervals in rotation,
which can be specified with the :samp:`physical.scale_interval` and :samp:`physical.decay_interval`
options. The number of parameters in the absorption surface is determined by the
highest order of spherical harmonics function used, controlled by :samp:`physical.lmax`
(recommended to be no higher than 6, 4 by default). There is also a weak
:samp:`physical.decay_restraint` and strong :samp:`physical.surface_weight` to
restrain the parameters of the decay and absorption terms towards zero.
The physical model is suitable for most datasets, although the absorption correction
should be turned off for datasets with low reciprocal space coverage.

The KB model applies a single scale factor and single B-factor to the whole
dataset (B-factor can be turned off with :samp:`decay_term=False`). This is
only suitable for very thin wedge/single-image datasets. If the KB model is
used, it may be necessary to set :samp:`full_matrix=False`, as the full matrix
minimisation round can be unstable depending on the number of reflections per
dataset.

The array model consists of up to three components. The first (
:samp:`array.decay_correction`), consists of a smoothly varying correction
calculated over a 2D grid of parameters, as a function of rotation vs resolution
(d-value). The parameter interval in rotation is controlled by
:samp:`array.decay_interval`, while the number of resolution bins is
controlled by :samp:`array.n_resolution_bins`.
The second (:samp:`array.absorption_correction`) consists of a smoothly
varying correction calculated over a 3D grid of parameters, as a function of
rotation, x and y position of the measured reflection on the detector. The spacing
in rotation is the same as the decay correction, while the detector beginning is
controlled with :samp:`array.n_absorption_bins`.
Finally, an :samp:`array.modulation_correction` can be applied, which is a
smooth 2D correction as a function of x and y position, controlled with
:samp:`array.n_modulation_bins`, although this is off by default.
The array model is only suitable for wide-rotation datasets with a high
number of reflections and it should be tested whether the absorption
correction is suitable, as it may lead to overparameterisation.


Advanced use - Choosing reflections to use for minimisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To minimise the scaling model, a subset of reflections are used for efficiency.
Four methods are available with the following command:
:samp:`reflection_selection.method=auto quasi_random intensity_ranges use_all`.

By default, the auto method uses the quasi_random selection algorithm, with
automatically determined parameters based on the dataset properties. If the
dataset is small (<20k reflections), the :samp:`use_all` option is selected.

For each dataset, the quasi_random algorithm chooses reflection groups that
have a high connectedness across different areas of reciprocal space,
across all resolution shells. In multi-dataset scaling, a separate selection
is also made to find reflection groups that have a high connectedness across
the datasets (choosing from groups with an average I/sigma above a cutoff).
The parameters of the algorithm are therefore controllable with the following
options, if one explicitly chooses :samp:`reflection_selection.method=quasi_random`:
:samp:`quasi_random.min_per_area`, :samp:`quasi_random.n_resolution_bins`,
:samp:`quasi_random.multi_dataset.min_per_dataset` and
:samp:`quasi_random.multi_dataset.Isigma_cutoff`. The :samp:`auto` option sets these
parameters in order to give sufficient connectedness across reciprocal space/datasets
depending on the size of the dataset, number or parameters and number of datasets.

The :samp:`intensity_ranges` option chooses intensities between a range of
normalised intensities (:samp:`E2_range`), between a range of I/sigma (:samp:`Isigma_range`)
and between a resolution range (:samp:`d_range`). This will typically select
around 1/3 of all reflections.

The :samp:`use_all` method simply uses all suitable reflections for scaling model
minimisation, but may be prohibitively slow and memory-intensive for large datasets.


.. _aimless: http://www.ccp4.ac.uk/html/aimless.html
.. _xscale: http://xds.mpimf-heidelberg.mpg.de/html_doc/xscale_program.html

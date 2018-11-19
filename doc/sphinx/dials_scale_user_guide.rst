User guide for scaling data with DIALS
======================================

This document aims to provide an in-depth guide on how to use several
of the dials.scale command line options. It should be considered an 'expert'
level guide, and the reader is encouraged to first read the 
'scaling beta lactamase' tutorial for an overview of scaling in dials.
This guide includes:
- how to customise the scaling models and general tips
- how to exclude data after a first round of scaling
- some tips for how to help performance when scaling large datasets

Guide to the different scaling models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
There are three available scaling models available in dials.scale, accessible
by the command line option :samp:`model = physical array KB *auto`.
The physical model is similar to the scaling model used in the program aimless_,
the array model is based on the approach taken in xscale_, while the KB model is
a simple two-component model suitable for still-image datasets or very small
rotation datasets (~ < 3 degrees).

The auto option automatically chooses a default model and sensible parameterisation
based on the oscillation range of the experiment. model=auto will choose the
physical model unless the oscillation range is < 1.0 degree, when the KB model
will be chosen. The auto parameterisation rules are given at the bottom of this
section.

The physical model consists of up to three components - a smoothly varying
scale term, a smoothly varying B-factor term and an absorption surface
correction (all on by default). These are turned on/off with the command line options
:samp:`scale_term=True/False decay_term=True/False absorption_term=True/False`.
The smoothly varying terms have a parameter at regular intervals in rotation,
which can be specified with the :samp:`scale_interval` and :samp:`decay_interval`
options. The number of parameters in the absorption surface is determined by the
highest order of spherical harmonics function used, controlled by :samp:`lmax`
(recommended to be no higher than 6, 4 by default). There is also a weak
:samp:`decay_restraint` and strong :samp:`surface_weight` to restrain the
parameters of the decay and absorption terms towards 0.
The physical model is suitable for most datasets, although the absorption_term
should be turned off for datasets with low reciprocal space coverage.

The KB model applies a single scale factor and single B factor to the whole
dataset (B-factor can be turned off with :samp:`decay_term=False`). This is
only suitable for very thin wedge/single-image datasets. If the KB model is
used, it may be necessary to set :samp:`full_matrix=False`, as the full matrix
minimisation round can be unstable depending on the number of reflections per
dataset.

The array model consists of up to three components. The first (:samp:`decay_term`),
consists of a smoothly varying correction calculated over a 2D grid of
parameters, as a function of rotation vs resolution (d-value). The parameter
interval in rotation is controlled by :samp:`decay_interval`, while the number
of resolution bins is controlled by :samp:`n_resolution_bins`.
The second (:samp:`absorption_term`) consists of a smoothly varying correction
calculated over a 3D grid of parameters, as a function of rotation, x and y
position of the measured reflection on the detector, controlled with
:samp:`decay_interval` and :samp:`n_absorption_bins`.
Finally, a :samp:`modulation_term` can be applied, which is a smooth 2D correction as a
function of x and y position, controlled with :samp:`n_modulation_bins`,
although this is off by default. The array model is only suitable for
wide-rotation datasets with a high number of reflections and it should be tested
whether the absorption term is suitable, as it may lead to overparameterisation.

| Auto model rules:
| if oscillation range < 1.0 degrees - use KB model, else use physical model
| if oscillation range < 60.0 degrees, absorption_term = False
| scale and decay parameter intervals based on oscillation range:
| if 1.0 <= oscillation range < 10.0 degrees; intervals 2.0, 3.0
| if 10.0 <= oscillation range < 25.0 degrees; intervals 4.0, 5.0
| if 25.0 <= oscillation range < 90.0 degrees; intervals 8.0, 10.0
| if oscillation range >= 90.0 degrees; intervals 15.0, 20.0

Excluding data/image handling after initial scaling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
After a first round of scaling, it may be apparant that there are datasets,
or regions of datasets, that are in poor agreement with the rest of the
dataset, and it would be advantageous to remove this data and rescale (this is
particularly relevant for thin-wedge rotation datasets and still image datasets).
dials.scale provides two options for removing data, depending on whether
one wishes to exclude a whole dataset or only part of a dataset.

To exclude whole datasets, we can take advantage of the fact that unique
experiment identifiers are assigned to the datasets as labels - these are
currently assigned as strings of integers i.e. '0', '1', '2' etc. (these
can also be assigned manually with :samp:`dev.dials.assign_experiment_identifiers`)
The assignment of the identifiers can be seen in the scaling log / terminal
output, in one of the first lines of output::

  Dataset unique identifiers are ['0', '1', '2', '3']

To exclude datasets, one therefore uses the :samp:`exclude_datasets` option::

  dials.scale ...... exclude_datasets="0 2"

Alternatively, one can use the option :samp:`use_datasets`::

  dials.scale ...... use_datasets="1 3"

These datasets are removed at the start of the program before scaling occurs,
and will not be contained in the output :samp:`scaled.pickle` and
:samp:`scaled_experiments.json`.

To help with excluding parts of a dataset, image exclusion can be performed
using the command-line syntax :samp:`exclude_images="exp_id:start:stop"`. Here
exp_id is the experiment identifier (a string) indicating the dataset,
and start and stop are integers that define the image range to exclude (the
excluding region includes start and stop) i.e. to exclude images 101 to 200 from
experiment "0", one would use :samp:`exclude_images="0:101:200"`.

In the reflection_table, the reflections corresponding to these imags are
marked with the :samp:`user_excluded_for_scaling` flag, and the parameters of the
scaling models are adjusted to span the new image range. These data will not
be included in future scaling or data export, and further image exclusion
can be performed in subsequent scaling jobs.

Note that it is recommended to only exclude data at the beginning or end of a
sweep. One can use it to exclude data in the middle of a sweep, however care
must be taken that only a short image range is excluded. If the interior
excluded range is of the order of the scaling model parameter spacing, this can
cause the scaling model minimisation to fail. In this case it would be better to
split the experiment with :samp:`dials.slice_sweep` and then proceed with
excluding images at the edge of the new experiments.


Practicalities for large datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Depending on the computational resources available, scaling of large datasets
( > 1 million reflections) can become slow and memory intensive.
There are several options available for managing this.

The first option is separating the data in memory to allow blockwise calculations
and parallel processing, using the option :samp:`nproc=` (a value of 4 or 8 is probably a
reasonable choice).

One of the most intensive part of the algorithm is
full matrix minimisation, which is by default performed after a quicker LBFGS
minimisation round. One can set :samp:`full_matrix=False` to turn this off, however
no errors for the inverse scale factors will be determined. A compromise is
to set :samp:`full_matrix_max_iterations=1` to do at least one iteration.

A third option is to reduce the number of reflections used by the scaling
algorithm during minimisation. By default, a subset of reflections is chosen based on their
normalised intensities, with the default set chosen between E2 values of 0.8
and 5.0, which typically selects between 1/3 and 1/2 of the dataset. These limits
can be set with :samp:`E2_range=min,max`, or similary an :samp:`Isigma_range=min,max` or
:samp:`d_range=min,max` can be set to reduce the number of reflections
used to determine the scaling model. However, one should be
careful that the subset is representative of the whole dataset, and selecting
too few reflections will lead to overfitting of the subset and worse overall
merging statistics.

.. _aimless: http://www.ccp4.ac.uk/html/aimless.html
.. _xscale: http://xds.mpimf-heidelberg.mpg.de/html_doc/xscale_program.html
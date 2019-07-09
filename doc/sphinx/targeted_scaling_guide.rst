Scaling against a reference dataset with DIALS
==============================================

This document provides a guide on how to scale a dataset against a reference
dataset, referred to as the target dataset.
The target dataset can be either a dataset scaled with dials.scale, or an mtz
file containing a scaled dataset.
Any number of datasets can be scaled against the reference dataset at once,
giving the same scaling result as if each dataset were scaled as an independent
job. The only difference is that all scaled datasets would be output in one
scaled.refl and scaled.expt, which may be more or less convenient
for further processing.

Scaling against a dials reference dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this example, reference.refl and reference.expt are
from a dataset that has already been scaled with dials.scale. To scale another
dataset (datafiles integrated.refl, integrated.expt) against this
target/reference, one should use the following command::

  dials.scale only_target=True integrated.refl integrated.expt reference.refl reference.expt

This will scale the intensities of the dataset to agree as closely as possible
with the intensities of the reference dataset, and save the scaled dataset to
scaled.refl, scaled.expt (the reference files are unchanged).
The :samp:`only_target=True` command is important, else all the data will be
scaled together and output in a joint output file.

Scaling against a reference mtz file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this case, it is assumed that the intensity and variance columns of the mtz
file have already been scaled. Targeted scaling would be run with the following
command::

  dials.scale integrated.refl integrated.expt target_mtz=scaled.mtz

The targeted scaling algorithm is the same regardless of the target datafile type,
likewise the scaled dataset will be saved to scaled.refl and scaled.expt.


General considerations for suitable options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A common use case for scaling against a reference is to scale thin-wedge
datasets against a high quality full-sweep dataset. To give the best scaling, it
may be necessary to manually set the scaling model parameters: for more details
see the *In-depth guide to scaling options in DIALS*.
In the case of very thin wedge/stills datasets, or depending on the scientific question under investigation, it may be
suitable to set :samp:`model=KB`, to give a single global scale and relative B-factor
to each dataset. However, if significant intensity variation/decay is present in each
measurement, it may be best to use :samp:`model=physical`, setting :samp:`absorption_term=False`
and specifying values for :samp:`scale_interval` and :samp:`decay_interval`.

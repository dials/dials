Developer Tutorial
==================

DIALS has been designed to be extensible. If you've got a great idea for a new
integration algorithm, then, with a bit of work, you should be able to get it
running within the DIALS framework. The following tutorial applies for both
developers looking to put their algorithms within DIALS itself as well as for
those adding algorithms in their own installation.

Spot Finding
------------

Algorithms for computing the possible spot pixels used in spot finding should
implement the following interface. An example can be found under
"dials/extensions/kabsch_spotfinder_threshold_ext.py".

.. code-block:: python

  class SpotFinderThresholdIface(interface.Interface):

    def __init__(self, params, imageset):
      pass

    @interface.abstractmethod
    def compute_threshold(self, image, mask):
      pass

The algorithm is configured through phil parameters and the imageset.
The extension should also provide a :samp:`compute_threshold` method which
takes an image and a mask and should return a flex.bool array of the same
dimensions where :samp:`True` pixels are possible spot pixels and :samp:`False`
pixels are background.

Profile Modelling
-----------------

.. warning:: The profile modelling framework is not very mature and may be
   changed. I'll try and keep the documentation up-to-date.

Algorithms for computing the profile model should implement the following
interface. An example can be found under
"dials/extensions/gaussian_rs_profile_model_ext.py".

.. code-block:: python

  class ProfileModelCreatorIface(interface.Interface):

    @interface.abstractmethod
    def create(cls, params, experiments, reflections=None):
      pass

The profile algorithm should have a "create" method that is given phil
parameters, experiments and (optionally) reflections. In the absence of
reflections, the profile model should be constructable from phil parameters.
When reflections are present, the create method should be able to construct the
model from the input reflections and experiments.

The algorithm should return an instance of the ProfileModelList class (found in
:samp:`dials.algorithms.profile_model.model_list`). The profile model list 
should contain a list of objects which implement the interface given in
the :samp:`dials.algorithms.profile_model.interface` module and shown in brief
below.

.. code-block:: python

  class ProfileModelIface(object):

    @interface.abstractmethod
    def predict_reflections(self, experiment, **kwargs):
      pass

    @interface.abstractmethod
    def compute_bbox(self, experiment, reflections, **kwargs):
      pass

    @interface.abstractmethod
    def compute_partiality(self, experiment, reflections, **kwargs):
      pass

    @interface.abstractmethod
    def compute_mask(self, experiment, reflections, **kwargs):
      pass

    @interface.abstractmethod
    def dump(self):
      pass

The profile model should have methods for predicting the reflections, computing
the bounding box of reflections for a number of experiments, computing the
partiality of reflections and computing the foreground/background mask. Of these
the bounding box and mask methods are crucial for integration to work;
partiality is currently only used in reporting and can be a placeholder.

The extention should have the ability to dump the profile model to phil
parameters so that it can be input via a profile.phil file to, for example,
re-run integration with the same profile parameters.

Indexing
--------

FIXME

Refinement
----------

FIXME

Integration
-----------

Centroid algorithms
^^^^^^^^^^^^^^^^^^^

Algorithms for computing the reflection centroid should implement the following
interface. An example can be found under
"dials/extensions/simple_centroid_ext.py".

.. code-block:: python

  class CentroidIface(interface.Interface):

    def __init__(self, params, experiments):
      pass

    @interface.abstractmethod
    def compute_centroid(self, reflections):
      pass

The algorithm is configured through phil parameters and the list of experiments.
The extension should also provide a :samp:`compute_centroid` method which
takes a list of reflections with extracted shoebox data. The "shoebox" column of
the reflection table should contain a list of :samp:`dials.model.Shoebox` types.
The algorithm should fill the "xyzobs.px" column of the reflection table with
the observed centroid positions.

Background algorithms
^^^^^^^^^^^^^^^^^^^^^

Algorithms for computing the reflection background should implement the
following interface. An example can be found under
"dials/extensions/simple_background_ext.py".

.. code-block:: python

  class BackgroundIface(interface.Interface):

    def __init__(self, params, experiments):
      pass

    @interface.abstractmethod
    def compute_background(self, reflections):
      pass

The algorithm is configured through phil parameters and the list of experiments.
The extension should also provide a :samp:`compute_background` method which
takes a list of reflections with extracted shoebox data. The "shoebox" column of
the reflection table should contain a list of :samp:`dials.model.Shoebox` types.
The algorithm should fill the shoebox.background values and return the
reflection list.


.. Intensity algorithms
.. ^^^^^^^^^^^^^^^^^^^^

.. Algorithms for computing the reflection intensities should implement the
.. following interface. An example can be found under
.. "dials/extensions/summation_integration_ext.py".

.. .. code-block:: python

..   class IntensityIface(interface.Interface):

..     def __init__(self, params, experiments, profile_model):
..       pass

..     @interface.abstractmethod
..     def type(self, params, experiments):
..       pass

..     @interface.abstractmethod
..     def compute_intensity(self, reflections):
..       pass

.. The algorithm is configured through phil parameters, the list of experiments and
.. the list of profile models. The extension should also provide a
.. :samp:`@classmethod` named type which returns the type of integrator to use. The
.. supported return values for this function as shown below. Some algorithms may
.. choose to configure the appropriate type of integrator from the input phil
.. parameters and experiment list. Others may support only a single type of
.. integrator.

..  +----------+------------+-------------------------------+
..  | Value    | Experiment | Description                   |
..  +==========+============+===============================+
..  | 3d       | rotation   | 3D shoeboxes                  |
..  +----------+------------+-------------------------------+
..  | flat3d   | rotation   | 3D shoeboxes flattend         |
..  +----------+------------+-------------------------------+
..  | 2d       | rotation   | 2D partials                   |
..  +----------+------------+-------------------------------+
..  | single2d | rotation   | 2D partials on a single image |
..  +----------+------------+-------------------------------+
..  | stills   | stills     | 2D partials on a single image |
..  +----------+------------+-------------------------------+

.. Finally, the extension should provide a :samp:`compute_intensity` method which
.. takes a list of reflections with extracted shoebox data. The algorithm should
.. fill the "intensity.prf.value" and "intensity.prf.variance" columns in the
.. reflection table and return it.

Deploying algorithms
--------------------

Within the DIALS project
^^^^^^^^^^^^^^^^^^^^^^^^

The DIALS project has the following layout.

.. code-block:: none

  dials
  |
  |-- algorithms
  |   |
  |   |-- integration
  |       |
  |       |-- sum
  |           |
  |           |-- ...
  |
  |-- interfaces
  |   |
  |   |-- ...
  |
  |-- extensions
      |
      |-- summation_integration_ext.py
      |
      |-- ...

Each algorithm should have its implementation encapsulated within a package in
the appropriate place. For example, summation integration is implemented within
the "dials.algorithms.integration.sum" package. The extension class, which
inherits from the appropriate interface, and configures and calls this algorithm
should then be placed in the "dials/extensions" folder with an appropriate name.
For example, the summation integration extension is placed in the module
"dials.extensions.summation_integration_ext". Modules within the dials.extension
package will be automatically loaded when searching for algorithms and any class
within these modules that inherits from an interface will be registered for use
within the DIALS command line programs.

Within external projects
^^^^^^^^^^^^^^^^^^^^^^^^

If you have a project containing algorithms written for use within DIALS that is
built using the cctbx build system, it is easy to make DIALS aware of your new
algorithms.

A typical project layout is shown below.

.. code-block:: none

  my_project
  |
  |-- algorithms
  |   |
  |   |-- integration
  |       |
  |       |-- my_algorithm
  |           |
  |           |-- ...
  |
  |-- extensions
      |
      |-- my_algorithm_intensity_ext.py

If your project has this directory structure, with an intensity algorithm
implementation within the "my_algorithm" directory and the extension class
(inheriting from the IntensityIface class) in the "extension" directory you can
make DIALS aware of your algorithm by adding the following code to the
libtbx_refresh.py scripy in the top level of your project. This will add the
extensions directory in your project to the list of directories searched when
loading available algorithms.

.. code-block:: python

  # libtbx_refresh.py

  from __future__ import division

  try:
    from dials.framework import env
    import libtbx.load_env
    from os.path import join
    path = libtbx.env.dist_path("my_project")
    env.cache.add(join(path, "extensions"))
  except Exception:
    pass

Running "libtbx.refresh" or "make reconf" will update your build. You can check
that your algorithm has been found properly by using the "dials.show_extensions"
command-line program which should show a list of extensions implementing each
interface with your algorithm listed with the other available algorithms.


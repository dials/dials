+++++++++++++++
Extending DIALS
+++++++++++++++

.. contents::
   :depth: 2

Entry points
============

DIALS uses `entry points <https://packaging.python.org/specifications/entry-points/>`_
to define points in the code that can be extended by external developers. These entry
points are defined in the ``dials`` and ``dxtbx`` ``libtbx_refresh.py`` files:

.. literalinclude:: ../../../../libtbx_refresh.py
   :language: python
   :lines: 4-25

Developers implementing their own extensions can register their extensions either in
their package ``libtbx.refresh.py`` if writing a cctbx-style package, or in their
package ``setup.py`` if using the
`setuptools <https://setuptools.readthedocs.io/en/latest/>`_ framework.

The list of installed plugins can be obtained by running the ``dials.plugins`` command::

  $ dials.plugins
  dials.index.basis_vector_search_strategy  Basis vector search strategies
   fft1d (dials.algorithms.indexing.basis_vector_search.strategies via libtbx.dials 0.0.0)
   fft3d (dials.algorithms.indexing.basis_vector_search.strategies via libtbx.dials 0.0.0)
   real_space_grid_search (dials.algorithms.indexing.basis_vector_search.strategies via libtbx.dials 0.0.0)

  dxtbx.scaling_model_ext  scaling models
   KB (dials.algorithms.scaling.model.scaling_model_ext via libtbx.dials 0.0.0)
   array (dials.algorithms.scaling.model.scaling_model_ext via libtbx.dials 0.0.0)
   physical (dials.algorithms.scaling.model.scaling_model_ext via libtbx.dials 0.0.0)

  dxtbx.profile_model  profile models
   gaussian_rs (dials.extensions.gaussian_rs_profile_model_ext via libtbx.dials 0.0.0)

  dials.index.lattice_search_strategy  Lattice search strategies
   low_res_spot_match (dials.algorithms.indexing.lattice_search_strategies via libtbx.dials 0.0.0)


Adding new format classes
=========================

``dxtbx`` now discovers format classes during configuration time instead of
at runtime. New format classes can either be added into the dxtbx/format
directory, registered by other python packages using the
'dxtbx.format' entry point, or installed by the user via the
'dxtbx.install_format' command.

To add a new format class to be distributed with ``dials``, please submit a pull request
to the `dxtbx repository <https://github.com/cctbx/dxtbx>`_.

To register format classes stored in ~/.dxtbx you need to run
'dxtbx.install_format -u' whenever you add or remove format classes.

Writing a new format class
--------------------------

The `dxtbx` format class framework enables beamline staff and users to easily add
support for new detector types and beamlines. In essence all that is needed is to
implement a Python class which extends the Format class to add some specific details
about this detector and the associated beamline/experimental environment.

In particular there are two groups of things which need to be
implemented - a static method named `understand` which will take a
look at the image and return True if it understands it, and a number of
class methods which need to override the construction of the `dxtbx` models.

`understand` Static Method
--------------------------

This method is the key to how the whole framework operates - you write code
which looks at the image to decide whether it is right for this class. If it is
not you must return False - i.e. if you are making a custom class for a given
detector serial number and it is given an image from a different detector.

Ideally your implementation will inherit from a similar Format class
and just apply further customizations.  Your implementation will be
chosen to read the image if it is the most customized, i.e. it derives
from the longest chain of ancestors, all of which claim to understand
the image.

Class Methods
-------------

The class methods need to use the built in factories to construct descriptions
of the experimental apparatus from the image, namely the goniometer, detector,
beam and scan. In many cases the "simple" model will be the best which is
often trivial. In other cases it may be more complex but will hopefully
correspond to an already existing factory method.


As an example, let's pretend your beamline has a "reversed" rotation axis. We can create
a new format class that correctly understands images from your beamline and instantiates
a goniometer model with a reversed direction goniometer:

.. code-block:: python

  from __future__ import absolute_import, division, print_function

  from dxtbx.format.FormatCBFMiniPilatus import FormatCBFMiniPilatus

  class FormatCBFMiniPilatusMyBeamline(FormatCBFMiniPilatus):
      """A class for reading mini CBF format Pilatus images for MyBeamline."""

      @staticmethod
      def understand(image_file):
          """Check the detector serial number to check it is from MyBeamline."""

          header = FormatCBFMiniPilatus.get_cbf_header(image_file)
          for record in header.split("\n"):
              if (
                  "# Detector" in record
                  and "PILATUS" in record
                  and "S/N 42-4242" in header
              ):
                  return True
          return False

      def _goniometer(self):
          """Return a model for a simple single-axis reversed direction goniometer."""

          return self._goniometer_factory.single_axis_reverse()


We can then register this format class in the ``libtbx_refresh.py`` file of our local
``myproject`` ``cctbx`` package:

.. code-block:: python

  import libtbx.pkg_utils
  libtbx.pkg_utils.define_entry_points(
      {
          "dxtbx.format": [
              "FormatCBFMiniPilatusMyBeamline:FormatCBFMiniPilatus = myproject.my_format_module:FormatCBFMiniPilatusMyBeamline",
          ],
      }
  )

More generally, the format of an entry point for dxtbx.format is::

    "FormatMyClass:FormatBaseClass1,FormatBaseClass2 = myproject.myformat:FormatMyClass"

Format classes must be named 'Format*', and must inherit either from
other format classes or from the top-level format class, 'Format'.
Base classes must be given as their original name and must therefore not
contain '.'s.

To view the full hierarchy of registered format classes, run the command
``dxtbx.show_registry``::

  $ dxtbx.show_registry
  Showing hierarchy of classes in the dxtbx registry. The root classes are shown with depth of 1, and subclasses are shown indented and with a higher depth number.

  Depth  Class name
      0  Format
      1    FormatBruker
      2      FormatBrukerFixedChi
      2      FormatBrukerPhotonII
      ...


Extending dials.index
=====================

``dials.index`` defines two possible entry points,
``dials.index.basis_vector_search_strategy`` and
``dials.index.lattice_search_strategy``.


Basis vector search strategies
------------------------------
The ``dials.index.basis_vector_search_strategy`` entry point can be used to extend the
list of possible basis vector search strategies available in DIALS, by delegating the
search for a list of possible real space basis vectors to a strategy. DIALS currently
includes the `fft1d`, `fft3d` and `real_space_grid_search` strategies. A basis vector
search strategy should inherit from the class
:class:`dials.algorithms.indexing.basis_vector_search.strategies.Strategy` and provide
an implementation of the ``find_basis_vectors`` method.

.. code-block:: python

  from libtbx import phil
  from dials.algorithms.indexing.basis_vector_search.strategies import Strategy

  mystrategy_phil_str = """\
  magic_parameter = 42
      .help = "This is a magic parameter."
      .type = float
  """

  phil_scope = phil.parse(mystrategy_phil_str)

  class MyStrategy(Strategy):
      """Basis vector search using my magic algorithm."""

      def find_basis_vectors(self, reciprocal_lattice_vectors):
          """Find a list of likely basis vectors.

          Args:
              reciprocal_lattice_vectors (scitbx.array_family.flex.vec3_double):
                  The list of reciprocal lattice vectors to search for periodicity.

          Returns:
              A tuple containing the list of basis vectors and a flex.bool array
              identifying which reflections were used in indexing.

          """
          used_in_indexing = flex.bool(reciprocal_lattice_vectors.size(), True)
          # determine the list of candidate_basis_vectors
          ...
      return candidate_basis_vectors, used_in_indexing


We can now register this new basis vector search strategy in the ``libtbx_refresh.py``
file of our local ``myproject`` package:

.. code-block:: python

  import libtbx.pkg_utils
  libtbx.pkg_utils.define_entry_points(
      {
          "dials.index.basis_vector_search_strategy": [
              "mystrategy = myproject.mystrategy:MyStrategy",
          ],
      }
  )

Lattice search strategies
-------------------------

An alternative entry point into dials.index is
``dials.index.lattice_search_strategy``, where the entire crystal model search is
delegated to the strategy.

.. code-block:: python

  from libtbx import phil
  from dials.algorithms.indexing.lattice_search_strategies import Strategy

  mystrategy_phil_str = """\
  magic_parameter = 42
      .help = "This is a magic parameter."
      .type = float
  """

  class MyLatticeSearch(Strategy):
      """My magic lattice search strategy."""

      phil_scope = phil.parse(mystrategy_phil_str)

      def find_crystal_models(self, reflections, experiments):
          """Find a list of candidate crystal models.

          Args:
              reflections (dials.array_family.flex.reflection_table):
                  The found spots centroids and associated data

              experiments (dxtbx.model.experiment_list.ExperimentList):
                  The experimental geometry models

          Returns:
              A list of candidate crystal models.

          """
          # determine the list of candidate_crystal_models
          return candidate_crystal_models


As above, register this new lattice search strategy in the ``libtbx_refresh.py``
file of our local ``myproject`` package:

.. code-block:: python

  import libtbx.pkg_utils
  libtbx.pkg_utils.define_entry_points(
      {
          "dials.index.lattice_search_strategy": [
              "mylatticesearch = myproject.mylatticesearch:MyLatticeSearch",
          ],
      }
  )


Extending profile models
========================


Extending dials.scale
=====================

`dials.scale` can be extended by defining new scaling models using the
entry point ``dxtbx.scaling_model_ext``.


Defining a scaling model
------------------------
A new scaling model can be defined, which should inherit from the class
:class:`dials.algorithms.scaling.model.model.ScalingModelBase`. A new
scaling model must define the `from_dict`, `from_data` and
`configure_components` methods, and should also define an `__init__` method. The model
must also define `consecutive_refinement_order` to indicate which order the components
should be refined for the consecutive scaling mode.
The scaling model must be composed of multiplicative components, which must
inherit from
:class:`dials.algorithms.scaling.model.components.scale_components.ScaleComponentBase`.

.. code-block:: python

  from libtbx import phil
  from scitbx.array_family import flex
  from dials.algorithms.scaling.model.model import ScalingModelBase
  from mypath.components import SpecialComponent

  mymodel_phil_str = """\
  special_correction = True
      .help = "Option to toggle the special correction."
      .type = bool
  """

  class MyScalingModel(ScalingModelBase):
      """My scaling model."""

      id_ = "modelname"

      phil_scope = phil.parse(mymodel_phil_str)

      def __init__(self, parameters_dict, configdict, is_scaled=False):
          super(MyScalingModel, self).__init__(configdict, is_scaled)
          if "special" in configdict["corrections"]:
              self._components["special"] = SpecialComponent(
                  parameters_dict["special"]["parameters"],
                  parameters_dict["special"]["parameter_esds"],
              )

      @classmethod
      def from_dict(cls, obj):
          """Create a MyScalingModel from a dictionary."""
          configdict = obj["configuration_parameters"]
          is_scaled = obj["is_scaled"]
          if "special" in configdict["corrections"]:
              parameters = flex.double(obj["special"]["parameters"])
              if "est_standard_devs" in obj["special"]:
                  parameter_esds = flex.double(obj["special"]["est_standard_devs"])
          parameters_dict = {"special : {"parameters" : parameters, "parameter_esds" : parameter_esds}}
          return cls(parameters_dict, configdict, is_scaled)

      @classmethod
      def from_data(cls, params, experiment, reflection_table):
          """Create the MycalingModel from data."""
          configdict = {"corrections": []}
          parameters_dict = {}

          if params.modelname.special_correction:
              configdict["corrections"].append("special")
              parameters_dict["special"] = {
                  "parameters": flex.double([1.0, 1.0, 1.0]),
                  "parameter_esds": None,
              }
          configdict["important_number"] = len(reflection_table)

          return cls(parameters_dict, configdict)

      def configure_components(self, reflection_table, experiment, params):
          """Add the required reflection table data to the model components."""
          if "special" in self.components:
              self.components["special"].data = {"d": reflection_table["d"]}

      def consecutive_refinement_order(self):
          "A nested list of the refinement order".
          return [["special"]]

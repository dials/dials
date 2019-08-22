+++++++++++++++
Extending DIALS
+++++++++++++++

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

    "FormatMyClass:FormatBaseClass1,FormatBaseClass2 = myproject.my_format_module:my"

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
list of possible basis vector search strategies available in DIALS. DIALS currently
includes the `fft1d`, `fft3d` and `real_space_grid_search` strategies. A basis vector
search strategy should inherit from the class
:class:`dials.algorithms.indexing.basis_vector_search.strategies.Strategy` and provide
an implementation of the ``find_basis_vectors`` method.

.. code-block:: python

  from libtbx import phil
  from dials.algorithms.indexing.basis_vector_search_strategy.strategies import Strategy

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


Extending profile models
========================


Extending dials.scale
=====================

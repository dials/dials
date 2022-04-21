dials.cosym
===========

Introduction
------------

.. python_string:: dials.command_line.cosym.help_message


Basic parameters
----------------

.. phil:: dials.command_line.cosym.phil_scope
   :expert-level: 0
   :attributes-level: 0


Example
-------

Run ``dials.cosym``, providing the integrated experiment (``.expt``) and reflection (``.refl``)
files as input:

.. dials_tutorial_include:: multi_crystal/dials.cosym.cmd

The first step is to analyse the metric symmetry of the input unit cells, and perform
hierarchical unit cell clustering to identify any outlier datasets that aren't consistent
with the rest of the datasets. The largest common cluster is carried forward in the
analysis. You can modify the threshold that is used for determining outliers by setting
the ``unit_cell_clustering.threshold`` parameter.

.. dials_tutorial_include:: multi_crystal/dials.cosym.log
   :start-at: Hierarchical clustering of unit cells
   :end-at: Highest possible metric symmetry and unit cell using LePage

In this case, the unit cell analysis found 1 cluster of 4 datasets in :math:`P\,4\,2\,2`.
As a result, all datasets will be carried forward for symmetry analysis.

Each dataset is then normalised using maximum likelihood isotropic Wilson scaling,
reporting an estimated Wilson :math:`B` value and scale factor for each dataset:

.. dials_tutorial_include:: multi_crystal/dials.cosym.log
   :start-at: Normalising intensities for dataset
   :end-at: -2.67

A high resolution cutoff is then determined by analysis of CC½ and <I>/<σ(I)> as a
function of resolution:

.. dials_tutorial_include:: multi_crystal/dials.cosym.log
   :start-at: Estimation of resolution
   :end-at: reflections with d

Next, the program performs automatic determination of the number of dimensions for
analysis. This calculates the functional of equation 2 of
`Gildea, R. J. & Winter, G. (2018)`_ for each dimension from 2 up to the number of
symmetry operations in the lattice group. The analysis needs to use sufficient
dimensions to be able to separate any indexing ambiguities that may be present, but
using too many dimensions reduces the sensitivity of the procedure. In this case, it is
determined that 8 dimensions will be used for the analysis:

.. dials_tutorial_include:: multi_crystal/dials.cosym.log
   :start-at: Automatic determination of number of dimensions for analysis
   :end-at: Best number of dimensions

Once the analysis has been performed in the appropriate number of dimensions, the
results are analysed to score all possible symmetry elements, using algorithms similar
to those of `POINTLESS`_, using the fact that the angles between vectors represent
genuine systematic differences between datasets, as made clear by equation 5 of
`Diederichs, K. (2017)`_.

.. dials_tutorial_include:: multi_crystal/dials.cosym.log
   :start-at: Scoring individual symmetry elements
   :end-before: Scoring all possible sub-groups

Scores for the possible Laue groups are obtained by analysing the scores for the
symmetry elements that are present or absent from each group, and the groups are ranked
by their likelihood.

.. dials_tutorial_include:: multi_crystal/dials.cosym.log
   :start-at: Scoring all possible sub-groups
   :end-at: Laue group confidence

The program then concludes by reporting any reindexing operations that are necessary to
ensure consistent indexing between datasets. In this case, no indexing ambiguity is
present, so the reindexing operator is simply the identity operator for all datasets.

.. dials_tutorial_include:: multi_crystal/dials.cosym.log
   :start-at: Reindexing operators
   :end-before: Writing html report

The correctly reindexed experiments and reflections are then saved to file, along with
a `HTML report <https://dials.github.io/tutorial_data/master/multi_crystal/dials.cosym.html>`_:

.. dials_tutorial_include:: multi_crystal/dials.cosym.log
   :start-at: Writing html report
   :end-at: Saving reindexed reflections

The full log file can be viewed here:

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    .. dials_tutorial_include:: multi_crystal/dials.cosym.log


Full parameter definitions
--------------------------

.. phil:: dials.command_line.cosym.phil_scope
   :expert-level: 2
   :attributes-level: 2

.. _`Gildea, R. J. & Winter, G. (2018)`: https://doi.org/10.1107/S2059798318002978
.. _`Diederichs, K. (2017)`: https://doi.org/10.1107/S2059798317000699
.. _`POINTLESS`: https://doi.org/10.1107/S090744491003982X

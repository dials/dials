.. raw:: html

  <a href="https://dials.github.io/dials-2.2/documentation/tutorials/correcting_poor_initial_geometry_tutorial.html" class="new-documentation">
  This tutorial requires a DIALS 3 installation.<br/>
  Please click here to go to the tutorial for DIALS 2.2.
  </a>

DPF3 Part 1: Correcting poor initial geometry
=============================================

.. highlight:: none

Introduction
------------

The following example uses a dataset kindly provided by Wolfram Tempel, which
was collected at beamline 19-ID at the APS. This dataset is available for
download from |DPF3|.

.. |DPF3| image:: https://zenodo.org/badge/doi/10.5281/zenodo.45756.svg
          :target: https://doi.org/10.5281/zenodo.45756

This is a challenging dataset to process. There are a combination of problems,
including:

* A 'reversed' rotation axis
* Incorrect beam centre recorded in the image headers
* Split spots
* Multiple lattices
* Systematically weak spots that may correspond to pseudocentring

Despite these issues, the diffraction data is of reasonable quality and was
used to `solve the structure`_ after processing by XDS.

.. _solve the structure: http://www.rcsb.org/pdb/explore/explore.do?structureId=5I3L

In the first part of this tutorial we will look at how to use the DIALS toolkit
to address the poor initial model for the experimental geometry, which leads to
problems with indexing. The first listed problem, namely the inverted rotation
axis, is trivially dealt with. However the incorrect beam centre is
particularly pernicious in this case. Rather than resulting in an outright
failure to process, we instead obtain an incorrect indexing solution. If we
were being careless, this could have lead to the integration of a useless
data set.

This tutorial is a cautionary tale, the moral of which is that the user should
employ the diagnostic tools at their disposal and to think about the output of
the programs they run.

Import
------

The dataset consists of a tar archive of bz2-compressed images. DIALS can read
the compressed images directly, however we need to extract the archive first::

  tar xvf DPF3_247398.tar

At this point we have no reason not to trust the image headers. We shall just
go ahead and import the whole sequence as normal::

  dials.import x247398/t1.0*.img.bz2

This produces the file :file:`imported.expt`, containing an initial model for
the beamline geometry. You can inspect this model using :program:`dials.show`::

  dials.show imported.expt

Note how the goniometer rotation axis is given by ``{-1,0,0}`` rather than
``{1,0,0}``. This is because DIALS recognises that these images as being
from beamline 19-ID at the APS, which is known to have an inverted axis of
rotation compared with the more common direction. Settings such as inverse
:math:`\phi`, or vertical goniometers, can be the cause of problems with
processing data from currently unrecognised beamlines. As an aside, in such
a case we could force the rotation axis to be whatever we want like this::

  dials.import x247398/t1.0*.img.bz2 geometry.goniometer.axes=-1,0,0

Now that we have imported the data we should look at the images::

  dials.image_viewer imported.expt

Keen-eyed observers may already suspect that the beam centre is not correct,
however we shall continue through spot-finding as this is not affected by
the experimental geometry.

Find Spots
----------

Spot-finding in DIALS usually works well for Pilatus detectors, where
default assumptions about Poisson statistics of pixel counts, unity gain and
no point spread are accurate. These assumptions are not correct for CCD
detectors and this can be another source of problems with data processing.

To see the positions of strong pixels identified by the spot finding
algorithm, select the ``threshold`` button at the bottom of the image
viewer's ``Settings`` window. In this case, the default settings are not too
bad: the strong pixels clearly follow the diffraction pattern. So, we will
run the :program:`dials.find_spots` program with default settings::

  dials.find_spots imported.expt

After finding strong spots it is *always* worth viewing them using
:program:`dials.reciprocal_lattice_viewer`::

  dials.reciprocal_lattice_viewer imported.expt strong.refl

.. image:: https://dials.github.io/images/correcting_poor_initial_geometry_tutorial/dpf3_bad_found_spot.png

Presented with this view, we might already start to worry that something is
not quite right. Instead of neat columns of points corresponding to a
regular reciprocal lattice grid, the points are aligned in curved or even
spiral tracks. Extreme cases of this may indicate something grossly wrong,
like an inverted :math:`\phi` direction. In this instance the lattice is
still detectable, but distorted. We understand this as inaccurate mapping
from detector to reciprocal space. If the diffraction geometry model is
wrong, then :program:`dials.reciprocal_lattice_viewer` cannot calculate the
reciprocal lattice position for each centroid properly. This can cause
problems with indexing because that requires exactly the same step of
mapping centroid positions from detector to reciprocal space.

Notwithstanding these concerns, we press on into indexing.

Indexing
--------

::

  dials.index imported.expt strong.refl

It turns out that the reciprocal lattice positions were regular enough for
indexing to complete ('succeed' is the wrong word, as will become clear).
Remember that initial indexing uses fairly low resolution data only. At low
resolution the curved tracks of spots are straight enough to fit a lattice.
Macrocycles of refinement then extend the solution out to increasingly
high resolution. One might imagine this process as steps of unwarping the
distorted lattice from the centre outwards until a regular grid is formed.
Here's some output from the end of the indexing log::

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref  | RMSD_X  | RMSD_Y | RMSD_Z   |
  | id  |       | (px)    | (px)   | (images) |
  ---------------------------------------------
  | 0   | 20000 | 0.98416 | 1.6552 | 0.4345   |
  ---------------------------------------------

  Refined crystal models:
  model 1 (23317 reflections):
  Crystal:
      Unit cell: (118.74(3), 119.45(3), 126.41(3), 88.682(2), 89.257(3), 60.954(3))
      Space group: P 1

This is another point at which the experienced user may pause for thought.
Positional RMSDs of 0.98 and 1.66 pixels are rather bad. Good models for
synchrotron X-ray data
typically have values around 0.3 pixels or less. Split spots or other issues
with spot profiles may result in higher RMSDs for a solution that is still
correct, however we should always remain sceptical. Looking at the results
in :program:`dials.reciprocal_lattice_viewer` is instructive as ever::

  dials.reciprocal_lattice_viewer indexed.expt indexed.refl

.. image:: https://dials.github.io/images/correcting_poor_initial_geometry_tutorial/dpf3_bad_indexed.png

Refinement has done what it could to produce a regular lattice, but it is still
messy. We also see that the majority of the centroids remain unindexed, and
these are messier still.

.. image:: https://dials.github.io/images/correcting_poor_initial_geometry_tutorial/dpf3_bad_unindexed.png

At this point we should definitely heed the warnings and try to figure out
what happened and how to fix it. However, unfortunately a careless user could
go ahead and integrate with this model. Let's see what happens if we try
to refine compatible Bravais lattices::

  dials.refine_bravais_settings indexed.expt indexed.refl

::

  -------------------------------------------------------------------------------------------------------------------
  Solution Metric fit  rmsd    min/max cc #spots lattice                                 unit_cell  volume      cb_op
  -------------------------------------------------------------------------------------------------------------------
        12     1.7490 0.607   0.021/0.028  20000      hP 123.10 123.10 129.77  90.00  90.00 120.00 1702991    -a,b,-c
        11     1.7490 0.602  -0.043/0.057  20000      oC 123.82 215.23 130.75  90.00  90.00  90.00 3484342 b,-2*a+b,c
        10     1.7490 0.601   0.027/0.027  20000      mC 214.82 123.62 130.53  90.00  90.16  90.00 3466356  2*a-b,b,c
         9     1.3289 0.564  -0.043/0.091  20000      oC 120.83 212.43 128.48  90.00  90.00  90.00 3297608 a,-a+2*b,c
         8     1.3233 0.522  -0.043/0.040  20000      oC 127.04 215.20 132.57  90.00  90.00  90.00 3624346  a-b,a+b,c
         7     1.3289 0.485   0.091/0.091  20000      mC 119.74 210.39 127.26  90.00  89.00  90.00 3205385 a,-a+2*b,c
         6     1.3233 0.519 -0.043/-0.043  20000      mP 123.42 131.38 124.09  90.00 118.93  90.00 1761030     -a,c,b
         5     1.2564 0.437   0.033/0.033  20000      mC 123.60 210.10 129.34  90.00  90.97  90.00 3358310  a-b,a+b,c
  *      4     1.1535 0.353   0.057/0.057  20000      mC 118.64 205.80 125.24  90.00  88.65  90.00 3057089 b,-2*a+b,c
  *      3     1.0684 0.327 -0.031/-0.031  20000      mC 204.56 116.44 123.87  90.00  91.29  90.00 2949728  a-2*b,a,c
  *      2     0.6885 0.268   0.040/0.040  20000      mC 208.52 122.85 128.42  90.00  88.65  90.00 3288791 a+b,-a+b,c
  *      1     0.0000 0.195           -/-  19928      aP 118.97 119.67 126.65  88.68  89.25  60.96 1576060      a,b,c
  -------------------------------------------------------------------------------------------------------------------


It turns out that quite a few lattices can be forced to fit the putative
indexing solution, but again there are warnings everywhere that imply none
of these are right. First look at the ``Metric fit`` column. This value is
the `Le Page <https://doi.org/10.1107/S0021889882011959>`_ :math:`\delta`
value. For a correct indexing solution with a good dataset this should be a
small number, less than 0.1 say, such as in the
:doc:`processing_in_detail_betalactamase` tutorial. The ``rmsd`` column reports an
overall positional RMSD. Again, small numbers are better. Typically we would
look for a solution below a jump to higher values of RMSD. Here they are all
pretty bad, at around an order of magnitude larger than what we'd expect
from good data. Another clear indication that none of the symmetry operations
implied by the higher symmetry lattices is correct is given by the ``min/max
cc`` column. This reports the lowest and highest correlation coefficients
between the rough spot-finding intensities of subsets of reflections related
by symmetry elements of the ``lattice``. For a real solution without rather
extreme radiation damage or other scaling issues we would expect much larger
numbers than these, say >0.5 or so for both the ``min`` and ``max`` values.

Check indexing symmetry
-----------------------

The fact that none of the correlation coefficients is high is a hint that
although the spots we indexed may indeed be real, perhaps the indices are
shifted by some value. This would be equivalent to the beam centre latching
onto some very low resolution Bragg reflection rather than the direct beam
:math:`hkl = (0,0,0)`. DIALS offers a tool to check this. If we run::

  dials.check_indexing_symmetry indexed.expt indexed.refl grid=1

then all combinations of off-by-one offsets in :math:`h`, :math:`k` and :math:`l`
will be checked by testing correlation coefficients between sets of reflections
related by symmetry. Here the model crystal symmetry is :math:`P 1`, so we are
testing only the Friedel pairs. The results are printed as a table in the
output::

  Checking HKL origin:

  dH dK dL   Nref    CC
  -1 -1 -1   2996 0.171
  -1 -1  0   3151 0.241
  -1 -1  1   3147 0.256
  -1  0 -1   2924 0.159
  -1  0  0   3097 0.261
  -1  0  1   3232 0.266
  -1  1 -1   2729 0.134
  -1  1  0   2904 0.172
  -1  1  1   3139 0.136
   0  0  0   1573 -0.178
   1 -1 -1   2876 0.272
   1 -1  0   2992 0.331
   1 -1  1   3135 0.257
   1  0 -1   2851 0.254
   1  0  0   3005 0.265
   1  0  1   3156 0.339
   1  1 -1   2792 0.244
   1  1  0   3073 0.283
   1  1  1   3718 0.886

  Check symmetry operations on 23317 reflections:

                 Symop   Nref    CC
                 x,y,z  23317 0.996


In this case there is a much greater correlation coefficient for the shift
:math:`\delta h=1`, :math:`\delta k=1` and :math:`\delta l=1` than for all
others. In fact with nearly 90% correlation even in the unscaled, rough intensities
of the found spots, with no background subtraction, we can be very sure we
have found the right solution.

Although it is possible to apply the correction using :program:`dials.reindex`
like this::

  dials.reindex indexed.refl hkl_offset=1,1,1

it will be very difficult to take the result and continue to process the data.
There is a much better way to proceed.

Discover a better experimental model
------------------------------------

We have determined that there is a problem with indexing, which gives us a
mis-indexed solution. The typical culprit in such cases is a badly wrong
beam centre. DIALS provides the
:program:`dials.search_beam_position`, which can help out
here. This performs a search to improve the direct beam position using
the `methods <https://doi.org/10.1107%2FS0021889804005874>`_ originally
implemented in :program:`LABELIT`.

This sits in between the spot finding and the indexing operations, so that
we could have done::

  dials.search_beam_position strong.refl imported.expt n_macro_cycles=2

In particularly bad cases it may useful to perform this search iteratively.
Indeed that is what we have done here by requesting two macrocyles. The first
macrocycle was not sufficient to find the real beam centre, but it improved
the search enough that it could be found in the second round::

  Starting macro cycle 1
  Selecting subset of 10000 reflections for analysis
  Running DPS using 10000 reflections
  Found 6 solutions with max unit cell 167.93 Angstroms.
  Old beam centre: 159.98, 154.50 mm (1562.3, 1508.8 px)
  New beam centre: 159.76, 152.65 mm (1560.2, 1490.7 px)
  Shift: 0.22, 1.85 mm (2.1, 18.1 px)

  Starting macro cycle 2
  Selecting subset of 10000 reflections for analysis
  Running DPS using 10000 reflections
  Found 9 solutions with max unit cell 167.93 Angstroms.
  Old beam centre: 159.98, 154.50 mm (1562.3, 1508.8 px)
  New beam centre: 162.26, 153.39 mm (1584.6, 1498.0 px)
  Shift: -2.28, 1.11 mm (-22.3, 10.8 px)

Indexing with the corrected beam centre
---------------------------------------

::

  dials.index optimised.expt strong.refl

We now have a more convincing solution, which also indexes many more
reflections::

  RMSDs by experiment:
  ----------------------------------------------
  | Exp | Nref  | RMSD_X  | RMSD_Y  | RMSD_Z   |
  | id  |       | (px)    | (px)    | (images) |
  ----------------------------------------------
  | 0   | 20000 | 0.50948 | 0.56722 | 0.20791  |
  ----------------------------------------------

  Refined crystal models:
  model 1 (62669 reflections):
  Crystal:
      Unit cell: (56.259(2), 99.521(4), 121.212(5), 89.9765(8), 89.9914(11), 90.0028(11))
      Space group: P 1


The lattice looks orthorhombic, and indeed the top solution in the table
from :program:`dials.refine_bravais_settings` looks reasonable::

  dials.refine_bravais_settings indexed.expt indexed.refl

::

  --------------------------------------------------------------------------------------------------------------
  Solution Metric fit  rmsd  min/max cc #spots lattice                                 unit_cell volume    cb_op
  --------------------------------------------------------------------------------------------------------------
  *      5     0.0250 0.078 0.746/0.842  20000      oP  56.28  99.55 121.25  90.00  90.00  90.00 679306    a,b,c
  *      4     0.0237 0.078 0.746/0.746  20000      mP  56.29  99.57 121.27  90.00  90.00  90.00 679612    a,b,c
  *      3     0.0250 0.078 0.746/0.746  20000      mP  56.28 121.26  99.56  90.00  90.00  90.00 679516 -a,-c,-b
  *      2     0.0091 0.078 0.842/0.842  20000      mP  99.51  56.26 121.20  90.00  89.98  90.00 678570 -b,-a,-c
  *      1     0.0000 0.078         -/-  20000      aP  56.26  99.52 121.21  89.98  89.99  90.00 678646    a,b,c
  --------------------------------------------------------------------------------------------------------------

We may now go on to refine the solution and integrate, following the steps
outlined in the :doc:`processing_in_detail_betalactamase` tutorial. This is left
as an exercise for the reader. You can continue to solve
the structure in the primitive orthorhombic lattice, however model refinement
will present difficulties.

Could we have foreseen this difficulties as early as the indexing step in DIALS?
Can we circumvent them? These are the topics explored in the second part of this
tutorial at :doc:`centring_vs_pseudocentring`.

Conclusions
-----------

* Incorrect or wrongly-interpreted image headers are a fact of life. You will
  encounter these.
* When beam centre problems are suspected, try
  :program:`dials.search_beam_position`.
* :program:`dials.reciprocal_lattice_viewer` and
  :program:`dials.image_viewer` are excellent troubleshooting tools for all
  sorts of spot finding and indexing problems.
* Some issues manifest as outright failures in indexing, others are more
  insidious and may result in a misindexed solution.
* Look out for CCs to detect misindexed data, and remember
  :program:`dials.check_indexing_symmetry`.
* Always use the diagnostic tools!

Acknowledgements
^^^^^^^^^^^^^^^^

Thanks to Wolfram Tempel for making this dataset available and inspiring
the writing of this tutorial.

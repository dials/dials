##############################
MyD88\ :sup:`TIR` small wedges
##############################

.. highlight:: none

.. warning::

  This tutorial was prepared using DIALS version 3.5.4 for Linux, downloaded
  from :doc:`this site <../../../installation>`. Results may differ with other
  versions of the software.

Introduction
============

The structures of MyD88\ :sup:`TIR` domain higher-order assemblies were resolved
by MicroED and SFX, providing insight into Toll-like receptor signal
transduction. Details are published in this article:

.. pubmed:: 33972532 MyD88

DIALS was not used during the original work, but the datasets have been made
public. So we can show how to use DIALS to reproduce the published results.

Data
====

The datasets for this tutorial are available online at SBGrid:

`Clabbers, MT.B., H Xu, and X Zou. 2021. “Microcrystal Electron
Diffraction Data for: Myeloid Differentiation Primary Response
Protein MyD88. PDB Code 7beq.” SBGrid Data Bank. 814`__.

.. __: https://doi.org/10.15785/SBGRID/814

The downloaded data consists of a directory with eighteen incrementally numbered
subdirectories, ``814/{1..18}``, each of which contains diffraction images from
a single crystal. We recommend doing data processing for each dataset in its
own directory, separate from the images. For example, we can create a structure
that mirrors the path to the images. From a BASH command prompt we can do that
like this:

.. code-block:: bash

  mkdir -p dials-proc/{1..18}

.. _section-label-exploratory-analysis:

Exploratory analysis
====================

Before we dive in and try to process all eighteen datasets we should perform some
initial investigations to get a feel for the data and ensure the metadata
describing the diffraction experiments is being read correctly. We might as well
take the first dataset for this. To start, we change our working directory to
our processing area, and then list the contents of the dataset directory.

.. code-block:: bash

   cd dials-proc/1
   ls ../../814/1

We see that dataset 1 consists of 72 images with the extension ``.img``, with
filenames numbered sequentially from ``00001.img`` to ``00072.img``. If the DIALS
programs are in the ``PATH`` we can now check if DIALS can interpret one of these
images.

.. code-block:: bash

   dials.show ../../814/1/00001.img

The output shows us that indeed DIALS recognises this image. Near the top of the
output we see ``Format class: FormatSMVTimePix_SU_516x516``, indicating that there
is even a specific image reader for this image rather than DIALS falling back on a
generic reader. That's a good start, but if we look carefully through the rest
of the output we might spot a problem.

Detector distance
-----------------

.. code-block::

   distance: 2.193

Those with experience with DIALS might know that laboratory space distances in
the program are measured in millimetres [*]_. This is evident from an earlier line,
``pixel_size:{0.055,0.055}`` in which we see that the pixel size of the Timepix
detector is 0.055 mm, i.e. 55 µm, as expected.

.. [*] except for the wavelength, which is in Ångströms.

Here we see a general problem with the SMV_ file format used for these images. The
format definition is not rigorous, so there is no way to *be sure* what units are in
use. Despite that fact that almost all known SMV formats use millimetres for
detector distance, this is not enforced by any standard, and we have clearly
found an exception.

.. _SMV: https://strucbio.biologie.uni-konstanz.de/ccp4wiki/index.php/SMV_file_format

Let's assume the distance should be 2193 millimetres and import the full dataset:

.. code-block:: bash

   dials.import ../../814/1/*.img distance=2193

Looking at the metadata with ``dials.show show.imported.expt`` shows that the
distance from the headers is now overwritten to be 2193 mm. Let's now view the
images:

.. code-block:: bash

   dials.image_viewer imported.expt

There is a good description of functions available in the image viewer in other
tutorials, such as :doc:`Processing in Detail<../processing_in_detail_betalactamase>`.
Feel free to play with the settings. Nothing you do here will alter the experimental
geometry or affect further processing.

We see from the position of the blue cross in the centre of the region of low
angle scatter that the beam centre seems to be correctly recorded in the image
headers.

Tilt axis orientation
---------------------

The diffraction geometry metadata printed by ``dials.show imported.expt``
suggests that the orientation of the rotation axis is given by

.. code-block::

    Rotation axis:   {0.782563,-0.622571,0}

But should we trust this? It is not as immediately visible a thing as the beam centre,
but sometimes we *can* get a visual clue from the image viewer.
Thinking of the geometry of the diffraction experiment, we realise
that spots that are perpendicular to the rotation axis appear and disappear
rapidly during rotation of the sample. Conversely, spots located along the rotation
axis remain in the diffracting condition for a long time. Therefore, by clicking
through the diffraction images we can look for a direction in the images in which
spots seem to persist for a long time. Doing this quickly should produce a view similar
to this animation:

.. image:: https://dials.github.io/images/MyD88/diffraction_movie.gif
   :width: 50%
   :align: center

If we look carefully we see that the spots in the lower left
and upper right persist on more images than spots in the upper left and lower
right. Therefore, we expect the rotation axis to be approximately along the
lower left to upper right diagonal.

We can also get some idea of the axis orientation by stacking the images. It is helpful to alter
the ``Stack type`` on the ``Settings`` window first to select ``max``, and then
in the main image viewer window change the value of ``Stack`` from ``1`` to
``72``. The view now shows a composite image consisting of the maximum value at
each pixel position through the whole dataset.

.. image:: https://dials.github.io/images/MyD88/stack.png
   :width: 80%
   :align: center

It is possible to make out a diagonal line from bottom left to top right along
which there are few spots. Those spots that are present are those that persisted
for a long time in the animation. We'll take the position of one of these as an
anchor point to estimate the orientation of the rotation axis. At the bottom of
the image viewer is a status bar from which we can read information
like the pixel position of the cursor. Reading out the pixel position gives us

.. code-block::

   Readout 0: slow=112.660 / fast=334.812 pixels

Now hovering over the beam centre we see this is located at about 235 pixels in
the slow direction and 228 in the fast direction. Therefore the line from the
beam centre to the pixel position we found before is approximately
:math:`113 - 235 = -122` pixels in the slow direction and :math:`335 - 228 = 107`
pixels in the fast direction. Returning to the output of ``dials.show imported.expt``
we see

.. code-block::

   fast_axis: {1,0,0}
   slow_axis: {0,-1,0}

from which we can construct a putative rotation axis direction

.. math::

  \vec{\textrm{axis}} = 107 \times \begin{pmatrix}1\\0\\0\end{pmatrix}
   - 122 \times \begin{pmatrix}0\\-1\\0\end{pmatrix}
  = \begin{pmatrix}107\\122\\0\end{pmatrix}

We don't know if the direction of right handed rotation is along
this axis or its inverse though. This is something we'll come back to during indexing.
At the moment, we'll just re-import the data with this axis and move on to
spot-finding.

.. code-block:: bash

   dials.import ../../814/1/*.img distance=2193 goniometer.axis=107,122,0.0

.. note::
   There is currently no easy way to determine the rotation axis using the
   :doc:`dials.image_viewer<../../programs/dials_image_viewer>`, hence these manual
   steps. As with the detector distance and beam centre it is best if these
   things are carefully calibrated for the data collection and recorded with
   the images.

Spot-finding
------------

Finding appropriate spot-finding settings can be challenging for some electron diffraction
datasets with commonly-used types of integrating detectors. However, in this case
the Timepix is a counting detector with a gain of 1.0, avoiding issues with
improperly-modelled detector response. The default spot-finding settings as used
for X-ray photon counting detectors are also appropriate here. We can view the
effect of these settings in the :doc:`dials.image_viewer<../../programs/dials_image_viewer>`.
In this case we have ticked the box next to ``Threshold pixels`` in the ``Settings``
panel to show a diffraction image with the strong pixels marked up with a red
overlay.

.. image:: https://dials.github.io/images/MyD88/spot-finding.png
   :width: 80%
   :align: center

The strong diffraction spots are clearly being found, but we also see that
a big region of inelastic scatter around the direct beam
is picked up by the spot-finding algorithm. This might cause problems with
indexing. By using the ``Resolution`` reading at the bottom of the image viewer
we see that most of this occurs within a d-spacing of 20 Å or higher, so we'll
exclude that, but otherwise leave spot-finding settings as default:

.. code-block:: bash

   dials.find_spots imported.expt d_max=20

The log shows the number of strong pixels found on each image, and then the
number of spots these form, followed by how many make it through various
filtering steps. The log ends with an ASCII-art histogram:

.. code-block::

    Histogram of per-image spot count for imageset 0:
    633 spots found on 72 images (max 26 / bin)
                                        *
                                        *
        **                              *
        ***        * *  *   ** * *      **
        ***        * *  *  ***** *      **      *
        ***        * ** *  ***** * ** * ****  * *
        ***      *** ** * ****** * ** * ****  *** **    **   *
    ** ****     ****************** *****************  * ** * *
    ********  ************************************** **********
    ************************************************ ***********
    1                          image                          72

    --------------------------------------------------------------------------------

With X-ray datasets this histogram can often be used as a quick assessment of radiation
damage. With electron diffraction that can be more difficult, both because the
total rotation angle for the scan is usually smaller and because the Ewald sphere
for electron diffraction is very flat. This means that the variability of this
plot is rather high. Some orientations are close to zone axes and have many
spots in the diffracting condition, whereas other orientations produce
bands of fewer spots instead. We can explore the found spots now in the image
viewer with this command:

.. code-block:: bash

   dials.image_viewer imported.expt strong.refl

Another useful viewer is the 3D
:doc:`dials.reciprocal_lattice_viewer<../../programs/dials_reciprocal_lattice_viewer>`,
which we can now launch like this

.. code-block:: bash

   dials.reciprocal_lattice_viewer imported.expt strong.refl

.. image:: https://dials.github.io/images/MyD88/rlv-strong.png
   :width: 100%

Using the mouse we can rotate and zoom this view, and can easily find directions
showing a well-ordered lattice. This gives us confidence that indexing will
be successful.

.. _section-label-indexing:

Indexing
--------

During indexing we will find a lattice and then refine the diffraction geometry
to better fit the observed spots. One major difference between electron and X-ray
diffraction is that the wavelength is much shorter (0.0251 Å in this case
compared to ~1 Å typical for X-rays). As a result, the Ewald sphere is rather
flat for electron diffraction and can be approximated by an "Ewald plane". A
side effect of this is that simultaneous refinement of the detector distance and
the unit cell parameters is hardly possible. Changes in the distance can be offset
by a scaling of the cell volume with negligible
differences in the predicted spot positions. To avoid refinement wandering off
to give unreasonable values for both the cell and distance, we typically fix the latter
by adding the option ``detector.fix=distance`` to jobs that include
geometry refinement.

We may also have some reasonable doubts about the accuracy of our estimated
rotation axis orientation. Usually in DIALS the rotation axis is assumed to be
known very well and fixed to the laboratory frame, but we can change that behaviour by setting the
option ``goniometer.fix=None``. We will otherwise do indexing using the standard 3D FFT
algorithm and other parameters as default, so the command we need is:

.. code-block:: bash

   dials.index imported.expt strong.refl detector.fix=distance goniometer.fix=None

At the end of the log for this job we see a high proportion of indexed spots:

.. code-block::

  +------------+-------------+---------------+-------------+
  |   Imageset |   # indexed |   # unindexed | % indexed   |
  |------------+-------------+---------------+-------------|
  |          0 |         562 |            70 | 88.9%       |
  +------------+-------------+---------------+-------------+


and a little further up we see that the diffraction geometry model fits the
observed spots quite nicely:

.. code-block::

  RMSDs by experiment:
  +-------+--------+----------+----------+------------+
  |   Exp |   Nref |   RMSD_X |   RMSD_Y |     RMSD_Z |
  |    id |        |     (px) |     (px) |   (images) |
  |-------+--------+----------+----------+------------|
  |     0 |    500 |  0.50338 |  0.68044 |    0.57072 |
  +-------+--------+----------+----------+------------+

But before we continue, we must remember that we haven't resolved the handedness
of the rotation axis yet!

Direction of rotation
^^^^^^^^^^^^^^^^^^^^^

The flat Ewald sphere in electron diffraction not only causes increased
correlation between parameters in refinement, but it also makes it possible to
index a lattice with either direction of the rotation axis! This is explained in
detail in this paper:

.. pubmed:: 29872002 Electron diffraction

To investigate the alternative we need to invert the rotation axis on import and try
indexing again. We don't need to repeat the spot-finding though, as the experimental
geometry does not affect the spots that are found on each image.

.. code-block:: bash

  dials.import ../../814/1/*.img distance=2193 goniometer.axis=-107,-122,0.0
  dials.index imported.expt strong.refl detector.fix=distance goniometer.fix=None

In this case there is one more indexed spot:

.. code-block::

  +------------+-------------+---------------+-------------+
  |   Imageset |   # indexed |   # unindexed | % indexed   |
  |------------+-------------+---------------+-------------|
  |          0 |         563 |            69 | 89.1%       |
  +------------+-------------+---------------+-------------+

but most tellingly, the ``RMSD_Z`` value has significantly decreased:

.. code-block::

  RMSDs by experiment:
  +-------+--------+----------+----------+------------+
  |   Exp |   Nref |   RMSD_X |   RMSD_Y |     RMSD_Z |
  |    id |        |     (px) |     (px) |   (images) |
  |-------+--------+----------+----------+------------|
  |     0 |    528 |   0.5165 |  0.67284 |    0.28888 |
  +-------+--------+----------+----------+------------+

The RMSDs are the root mean square deviation between observed and predicted
spot positions for reflections used in refinement. For both this job and the
previous one, the positional RMSDs are less than 1 pixel
in both X (fast) and Y (slow) directions on the image. However, the RMSD in
the tilt angle direction dropped from about 0.6 images to less than 0.3. This is
the best indication we have to show that the rotation axis is actually inverted
from our original guess.

.. note::

  The effect of inverting the rotation axis in electron diffraction is rather
  subtle. As an aside we can investigate this with
  ``dials.reciprocal_lattice_viewer indexed.expt indexed.refl`` and then ticking
  the box next to ``Invert rotation axis``. We see that the orientation of the
  spots changes, but they still make a clear lattice either way.

In this case, we will continue with the inverted axis. As a recap and to ensure
we are working with the right files, the correct steps up to this point are:

.. code-block:: bash

  dials.import ../../814/1/*.img distance=2193 goniometer.axis=-107,-122,0.0
  dials.find_spots imported.expt d_max=20
  dials.index imported.expt strong.refl detector.fix=distance goniometer.fix=None

Determining lattice symmetry
----------------------------

Unless a space group is explicitly specified,
:doc:`dials.index<../../programs/dials_index>` will return the
best fitting triclinic (:math:`P\ 1`) solution. A separate program,
:doc:`dials.refine_bravais_settings<../../programs/dials_refine_bravais_settings>`,
can then be used to analyse the lattice symmetry and suggest a higher-symmetry point
group. As this also does geometry refinement, we need to ensure the detector
distance remains fixed:

.. code-block:: bash

   dials.refine_bravais_settings indexed.expt indexed.refl detector.fix=distance

In the output we see two solutions: the original triclinic solution and
a centred monoclinic lattice:

.. code-block::

  Chiral space groups corresponding to each Bravais lattice:
  aP: P1
  mC: C2
  +------------+--------------+--------+--------------+----------+-----------+------------------------------------------+----------+--------------+
  |   Solution |   Metric fit |   rmsd | min/max cc   |   #spots | lattice   | unit_cell                                |   volume | cb_op        |
  |------------+--------------+--------+--------------+----------+-----------+------------------------------------------+----------+--------------|
  |   *      2 |       1.4551 |  0.118 | 0.950/0.950  |      531 | mC        | 64.31  37.02 115.99  90.00 104.37  90.00 |   267525 | c,a,-a+2*b-c |
  |   *      1 |       0      |  0.047 | -/-          |      529 | aP        | 37.07  61.62  64.40  72.89  89.26  73.47 |   134380 | a,b,c        |
  +------------+--------------+--------+--------------+----------+-----------+------------------------------------------+----------+--------------+
  * = recommended solution

If we knew nothing about the crystal structure beforehand we might continue and process with
this solution. However, here we are going to "cheat" slightly and look
at the cell from the published paper. There it is given as

  =============== =====
  Cell dimensions
  =============== =====
  a, b, c (Å)     99.06, 31.01, 54.30
  α, β, γ (°)     90.00, 107.70, 90.00
  =============== =====

Clearly we have found a different cell! In fact, here ``dials.refine_bravais_settings`` found an
:math:`I\ 2` setting but erroneously reported it as :math:`C\ 2`. This
issue has been fixed (https://github.com/dials/dials/pull/1825) so that from the
next release of DIALS the table above will correctly read ``mI``. The reason
DIALS selects the :math:`I\ 2` setting here is because by default it favours
monoclinic centred cells that have β angles closer to 90°. However, we can change
that behaviour and run again:

.. code-block:: bash

   dials.refine_bravais_settings indexed.expt indexed.refl detector.fix=distance best_monoclinic_beta=False

and with that we get

.. code-block::

  Chiral space groups corresponding to each Bravais lattice:
  aP: P1
  mC: C2
  +------------+--------------+--------+--------------+----------+-----------+-------------------------------------------+----------+-----------+
  |   Solution |   Metric fit |   rmsd | min/max cc   |   #spots | lattice   | unit_cell                                 |   volume | cb_op     |
  |------------+--------------+--------+--------------+----------+-----------+-------------------------------------------+----------+-----------|
  |   *      2 |       1.4551 |  0.118 | 0.950/0.950  |      531 | mC        | 117.84  37.02  64.31  90.00 107.54  90.00 |   267525 | a-2*b,a,c |
  |   *      1 |       0      |  0.047 | -/-          |      529 | aP        | 37.07  61.62  64.40  72.89  89.26  73.47  |   134380 | a,b,c     |
  +------------+--------------+--------+--------------+----------+-----------+-------------------------------------------+----------+-----------+
  * = recommended solution

Now at least the β angle is about what we expect! However, the a, b and c axes
are all a bit too long. In fact, they are all about 20% higher than the published values.
The cell volume is too large. At the introduction to :ref:`section-label-indexing` we noted
that the cell volume and detector distance are highly correlated. It looks like
the presumed detector distance of 2193 mm is still not correct.

We will fix that in the next section, but first we have to reindex the reflections
to match the chosen Bravais lattice solution, ``bravais_setting_2.expt``.
To do that we take the change-of-basis operator from the solution table
and pass that into the :doc:`dials.reindex<../../programs/dials_reindex>`
program:

  .. code-block:: bash

    dials.reindex indexed.refl change_of_basis_op=a-2*b,a,c

This creates a file ``reindexed.refl`` that is compatible with our chosen
solution ``bravais_setting_2.expt``.

Refining the detector distance
------------------------------

In situations where the correct unit cell is known it *is* possible to refine
the detector distance. We can do this by providing a restraint to the known
unit cell. This allows refinement of the unit cell and detector parameters
simultaneously, while pushing the cell towards its ideal values, thus breaking the
degeneracy between these parameters. The strength of this "push" is adjustable,
so we have control over how much we want the data or the external target to
affect the refined unit cell values.

To set up a restraint we must write a file using
`PHIL syntax <https://cci.lbl.gov/docs/cctbx/doc_low_phil/>`_. The interface
to restraints is a bit awkward, but having written this once, we can then copy-and-paste
it for other uses, with
changes required only to the ``values`` and the ``sigmas``. So, in a text
editor, copy these lines and save the file as ``restraint.phil``:

.. code-block::

  refinement
  {
    parameterisation
    {
      crystal
      {
        unit_cell
        {
          restraints
          {
            tie_to_target
            {
              values=99.06,31.01,54.30,90,107.7,90
              sigmas=0.01,0.01,0.01,0.01,0.01,0.01
            }
          }
        }
      }
    }
  }

This describes a restraint to the known cell with reasonably strong ties, given
by fairly low ``sigma`` values. Once we have this set up, we can run refinement
like this:

.. code-block:: bash

   dials.refine bravais_setting_2.expt reindexed.refl restraint.phil scan_varying=False

We added ``scan_varying=False`` to ensure that only "scan static" refinement is
performed, otherwise we get one round of scan static refinement followed by a
round of scan-varying refinement. The latter allows the crystal parameters to
vary across the scan, but those smoothly-changing parameters are not affected by
the restraint. In this situation, we are not trying to get a sophisticated,
varying model for the diffraction geometry, but are just trying to correct the
wrong detector distance.

At the end of the log we see that the unit cell now looks correct:

.. code-block::

  Final refined crystal model:
  Crystal:
      Unit cell: 99.080(18), 31.010(18), 54.264(17), 90.0, 107.706(17), 90.0
      Space group: C 1 2 1

and the RMSDs still look about as good as before, so refinement appears to have
been successful. We can see the new detector distance by showing the output experiments
file:

.. code-block:: bash

  dials.show refined.expt

which contains the line:

.. code-block::

  distance: 1772.4

.. note::

  Repeating this procedure with all 18 datasets produces a distribution of
  detector distances between about 1740 and 1810 mm. The median of that distribution
  is at 1762 mm, so we might use *that* value for further processing.
  Of course it is much preferred if the detector distance is properly calibrated
  and stored along with the diffraction images!

We can now return to import with the correct detector distance, followed by the
other steps to get back to the correctly indexed cell. We will skip the
``dials.refine_bravais_settings`` and ``dials.reindex`` steps now that we have
determined the lattice symmetry, by selecting ``space_group=C2`` during the
indexing job. This will make it easier to script these steps for the other
datasets.

.. code-block:: bash

  dials.import ../../814/1/*.img distance=1762 goniometer.axis=-107,-122,0.0
  dials.find_spots imported.expt d_max=20
  dials.index imported.expt strong.refl detector.fix=distance goniometer.fix=None\
    space_group=C2 output.experiments=C2.expt output.reflections=C2.refl

Further refinement
------------------

After indexing we usually run :doc:`dials.refine<../../programs/dials_refine>` to
construct a more sophisticated model of the diffraction geometry prior to
integration. In particular, by default this will perform a round of scan-varying
refinement, in which the crystal model (unit cell and orientation) is allowed to
vary as a function of image number. For some electron diffraction datasets, for
which the direct beam position appears to drift during data collection, we
can also try to model scan-varying beam orientation, however we are not going
to try that in this tutorial. Following from the last ``dials.index`` job, our refinement
command is:

.. code-block:: bash

  dials.refine C2.expt C2.refl detector.fix=distance

The crystal model is printed at the end of the log:

.. code-block::

  Crystal:
      Unit cell: 94.71(13), 29.746(16), 51.73(12), 90.0, 107.61(14), 90.0
      Space group: C 1 2 1
      U matrix:  {{ 0.2692, -0.6800,  0.6820},
                  {-0.2568, -0.7332, -0.6296},
                  { 0.9282, -0.0057, -0.3720}}
      B matrix:  {{ 0.0106,  0.0000,  0.0000},
                  {-0.0000,  0.0336,  0.0000},
                  { 0.0034, -0.0000,  0.0203}}
      A = UB:    {{ 0.0051, -0.0229,  0.0138},
                  {-0.0048, -0.0246, -0.0128},
                  { 0.0086, -0.0002, -0.0075}}
      A sampled at 73 scan points

The unit cell printed here is the *static* cell refined during the first
macrocycle of refinement. We can tell that there is a *scan-varying* cell model
as well though, from the final line, ``A sampled at 73 scan points``.

Integration
-----------

Now we have a suitable model for the experiment, we can go ahead and integrate
the reflections. We won't use any special options here.

.. code-block:: bash

  dials.integrate refined.expt refined.refl

There is a table of output at the end of the log that provides some insight into
how well this proceeded. In particular, if there are very large numbers of
reflections that failed to integrate by either summation integration or profile
fitting then we should investigate. In this case though, everything looks okay.

.. code-block::

  +---------------------------------------+-----------+--------+--------+
  | Item                                  |   Overall |    Low |   High |
  |---------------------------------------+-----------+--------+--------|
  | dmin                                  |      2.03 |   5.5  |   2.03 |
  | dmax                                  |     45.15 |  45.15 |   2.06 |
  | number fully recorded                 |   4320    | 539    |   9    |
  | number partially recorded             |    960    | 113    |   2    |
  | number with invalid background pixels |    994    |   0    |  11    |
  | number with invalid foreground pixels |    402    |   0    |  11    |
  | number with overloaded pixels         |      0    |   0    |   0    |
  | number in powder rings                |      0    |   0    |   0    |
  | number processed with summation       |   4854    | 644    |   0    |
  | number processed with profile fitting |   4918    | 631    |   2    |
  | number failed in background modelling |      0    |   0    |   0    |
  | number failed in summation            |    402    |   0    |  11    |
  | number failed in profile fitting      |    338    |  13    |   9    |
  | ibg                                   |     35.73 |  94.92 |  13.84 |
  | i/sigi (summation)                    |      3.93 |  18.51 |   0    |
  | i/sigi (profile fitting)              |      6.02 |  25.58 |   0    |
  | cc prf                                |      0.99 |   0.99 |   0.99 |
  | cc_pearson sum/prf                    |      0.99 |   0.99 |   0    |
  | cc_spearman sum/prf                   |      0.76 |   0.91 |   0    |
  +---------------------------------------+-----------+--------+--------+

We can open the integration results in the image viewer, showing how well the
reflections are centred in their integration boxes.

.. code-block:: bash

  dials.image_viewer integrated.expt integrated.refl

.. image:: https://dials.github.io/images/MyD88/integration_shoeboxes.png
  :width: 80%
  :align: center

Here we see that the fairly large rocking curve for the reflections means there
are many overlapping shoeboxes. However, this is fine as long as the peak regions
do not overlap. We cannot easily see this here, but the limited number of failed
profile fits tells us the peak regions are well separated. The size of the
reflection shoebox model is given towards the top of ``dials.integrate.log``:

.. code-block::

  Using 456 / 457 reflections for sigma calculation
  Calculating E.S.D Beam Divergence.
  Calculating E.S.D Reflecting Range (mosaicity).
   sigma b: 0.004093 degrees
   sigma m: 1.104106 degrees

At this stage we would learn much more about the quality of the data from scaling
and merging. However this single dataset is very incomplete. We should try
integrating all the other datasets first and including them in the scaling job.

Scripted processing
===================

During the :ref:`section-label-exploratory-analysis` we came up with a reasonable set of
processing commands for dataset 1, while figuring out incorrect or missing
diffraction geometry metadata such as the detector distance and the rotation
axis orientation and direction. Rather than repeat all those steps manually
for the remaining datasets we will write a script to process them all in the
same way. This example uses a BASH shell script on Linux, but we could do
similar on other systems.

First we change to the directory above where we have been working on dataset 1,
and if not already done, we'll make separate processing directories for all the datasets

.. code-block:: bash

  cd ..
  mkdir -p {1..18}

Now we'll use a text editor to make a file, called ``process.sh`` for example,
and enter these lines:

.. code-block:: bash

  #!/bin/bash
  set -e
  for i in {1..18}
  do
    cd "$i"

    dials.import ../../814/"$i"/*.img goniometer.axis=-107,-122,0.0 distance=1762
    dials.find_spots imported.expt d_max=20 d_min=2.5
    dials.index imported.expt strong.refl detector.fix=distance goniometer.fix=None space_group=C2
    dials.refine indexed.expt indexed.refl detector.fix=distance
    dials.plot_scan_varying_model refined.expt
    dials.integrate refined.expt refined.refl prediction.d_min=2.5

    cd ..
  done

We can make this executable and run it like this:

.. code-block:: bash

  chmod +x process.sh
  ./process.sh

This should run through to completion (``set -e`` ensures that it will stop if
any one of the jobs fails with an error) and produce 18 integrated datasets.

Scaling
=======

We now want to combine the 18 datasets, scale them together, and calculate merging
statistics. We will use :doc:`dials.scale<../../programs/dials_scale>` for
this. A particularly helpful feature of this program is its capability to
automatically filter out bad parts of the combined dataset using the ΔCC½ metric.
By default this removes complete datasets, which is useful for snapshot and small
wedge serial crystallography where we tend to consider each dataset a unit and
are not paying much attention to changes *within* each dataset caused, for example,
by radiation damage. By contrast, in this case, we'd like to try filtering
individual bad images rather than complete datasets. The command we need for that
is the following:

.. code-block:: bash

  mkdir -p scale
  cd scale

  dials.scale ../{1..18}/integrated.{expt,refl}\
    filtering.method=deltacchalf\
    deltacchalf.mode=image_group\
    deltacchalf.group_size=1\
    d_min=3.0

  cd ..

The first three options passed to ``dials.scale`` set up the ΔCC½ filtering,
in ``image_group`` mode rather than the default ``dataset``, and we also set
the ``group_size`` to 1 rather than the
default 10 as these datasets are rather wide-sliced. The next option sets a
resolution limit. We have chosen 3.0 Å here to match the processing published in
the paper.

At the end of the ``dials.scale.log`` we see a table summarising the merging
statistics:

.. code-block::

                                               Overall    Low     High
  High resolution limit                           3.00    8.13    3.00
  Low resolution limit                           50.87   50.88    3.05
  Completeness                                   77.7    76.3    81.7
  Multiplicity                                   14.4    14.9    10.2
  I/sigma                                        12.6    28.6     4.4
  Rmerge(I)                                     0.385   0.247   1.238
  Rmerge(I+/-)                                  0.381   0.247   1.222
  Rmeas(I)                                      0.399   0.258   1.305
  Rmeas(I+/-)                                   0.408   0.264   1.344
  Rpim(I)                                       0.100   0.069   0.387
  Rpim(I+/-)                                    0.136   0.087   0.532
  CC half                                       0.993   0.994   0.326
  Anomalous completeness                         79.6    82.5    80.9
  Anomalous multiplicity                          7.6     8.6     5.3
  Anomalous correlation                        -0.267  -0.413  -0.047
  Anomalous slope                               1.337
  dF/F                                          0.158
  dI/s(dI)                                      1.058
  Total observations                            34189    1926    1184
  Total unique                                   2380     129     116

These, and many more, details are also saved to a HTML format report page. On Linux you
can usually open this up with the command

.. code-block:: bash

  xdg-open scale/dials.scale.html

This contains many useful plots, for example various statistics as a function of
resolution:

.. image:: https://dials.github.io/images/MyD88/cchalf.png
   :width: 49%

.. image:: https://dials.github.io/images/MyD88/ioversigi.png
   :width: 49%

.. image:: https://dials.github.io/images/MyD88/completeness.png
   :width: 49%

.. image:: https://dials.github.io/images/MyD88/multiplicity.png
   :width: 49%

Further processing
==================

Now we have processed the 18 datasets we want to export the data for structure
solution and refinement. We can export the unmerged but scaled intensities in
MTZ format (``scaled.mtz``) like this

.. code-block:: bash

  dials.export scale/scaled.expt scale/scaled.refl

However, many downstream steps require the merged intensities, which we can
create using the ``dials.merge`` program.

.. code-block:: bash

  dials.merge scale/scaled.expt scale/scaled.refl

This creates the file ``merged.mtz``, which we can use for structure solution
by molecular replacement.
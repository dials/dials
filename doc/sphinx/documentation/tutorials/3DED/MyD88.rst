##############################
MyD88\ :sup:`TIR` small wedges
##############################

.. highlight:: none

Introduction
============

The structures of MyD88\ :sup:`TIR` domain higher-order assemblies was resolved
by MicroED and SFX methods, providing insight into Toll-like recptor signal
transduction. Details are published in this article:

.. pubmed:: 33972532 MyD88

DIALS was not used during the original work, but the datasets have been made
public, so we can show how to use DIALS to reproduce those results.

Data
====

The data for this tutorial are available online at SBGrid:

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

Exploratory analysis
====================

Before we dive in and try to process all eighteen datasets we should perform some
initial investigations to get a feel for the data and ensure the metadata
describing the diffraction experiments are being read correctly. We might as well
take the first dataset for this. To start, we change our working directory to
our processing area, and then list the contents of the dataset directory.

.. code-block:: bash

   cd dials-proc/1
   ls ../../814/1

We see that dataset 1 consists of 16 images with the extension ``.img``, with
filenames numbered sequentially from ``00001.img`` to ``00016.img``. If the DIALS
programs are in the ``PATH`` we can check if DIALS can interpret one of these
images.

.. code-block:: bash

   dials.show ../../814/1/00001.img

The output shows us that indeed DIALS recognises this image. Near the top of the
output we see ``Format class: FormatSMVTimePix_SU_516x516``, indicating that there
is a specific image reader for this image rather than DIALS falling back on a
generic reader. That's a good start, but if we look carefully through the rest
of the output we might spot a problem.

Detector distance
-----------------

.. code-block::

   distance: 2.193

Those with experience with DIALS might know that distances in laboratory space
in the program are measured in millimetres [*]_. This is evident from an earlier line,
``pixel_size:{0.055,0.055}`` in which we see that the pixel size of the Timepix
detector is 0.055 mm, or 55 µm, as expected.

.. [*] except for the wavelength, which is in Ångströms.

Here we see a general problem with the SMV_ format used for these images. The
format definition is not rigorous, so there is no way to know what units are in
use. Despite that fact that almost all known SMV formats use millimetres for detector
distance, this is not enforced, and we have clearly found an exception.

.. _SMV: https://strucbio.biologie.uni-konstanz.de/ccp4wiki/index.php/SMV_file_format

Let's assume the distance should be 2193 millimetres and import the full dataset:

.. code-block:: bash

   dials.import ../../814/1/*.img distance=2193

Looking at the metadata with ``dials.show show.imported.expt`` shows that the
distance is now overwritten to be 2193 mm. At this point we can now view the images:

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

An important piece of metadata that is not so immediately obvious
is the orientation of the rotation axis. Sometimes we can get a visual clue about
this though. Thinking of the geometry of the diffraction experiment, we realise
that spots that are perpendicular to the rotation axis appear and disappear
rapidly during rotation of the sample. Conversely, spots located along the rotation
axis remain in the diffracting condition for a long time. Therefore, by clicking
through the diffraction images we can look for a direction in the images in which
spots seem to persist for a long time. Doing this should produce a view similar
to this animation:

.. image:: https://dials.github.io/images/MyD88/diffraction_movie.gif
   :width: 50%
   :align: center

It is a little tricky to see, but we notice that the spots in the lower left and
upper right do seem to persist for longer than spots in the upper left and lower
right. Therefore, we may expect the rotation axis to be approximately along the
lower left to upper right diagonal.

We can also get some idea of this by stacking the images. It is helpful to alter
the ``Stack type`` on the ``Settings`` window first to select ``max``, and then
in the main image viewer window change the value of ``Stack`` from ``1`` to
``16``. The view now shows a composite image consisting of the maximum value at
each pixel position through the whole dataset.

.. image:: https://dials.github.io/images/MyD88/stack.png
   :width: 80%
   :align: center

At the bottom of the image viewer is a status bar from which we can read information
like the pixel position of the cursor. Here the cursor is positioned around the
the diagonal line through the beam centre where spot density seems to be low. It
is not very obvious to see, but we will use this as a starting point to determine
the rotation axis. Reading out the pixel position gives us:

.. code-block::

   Readout 0: slow=75.660 / fast=359.000 pixels

Now hovering over the beam centre we see this is located at about 235 pixels in
the slow direction and 228 in the fast direction. Therefore the line from the
beam centre to the pixel position we found before is approximately
:math:`76 - 235 = -159` pixels in the slow direction and :math:`359 - 228 = 131`
pixels in the fast direction.

.. note::
   There is currently no easy way to determine the rotation axis using the
   :doc:`dials.image_viewer<../programs/dials_image_viewer>`, hence these manual
   steps. As with the detector distance and beam centre it is best if these
   things are carefully calibrated for the data collection and recorded with
   the images.
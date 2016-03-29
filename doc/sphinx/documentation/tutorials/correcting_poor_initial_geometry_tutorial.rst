Correcting poor initial geometry
================================

Introduction
------------

One of the most common problems faced by users of any data processing program
is incorrect information about the experimental geometry in the input. This
will often cause outright failures in indexing, but sometimes more subtle
effects are possible, such as misindexing. In such cases, it may be possible
to index and refine a lattice that on the face of it looks reasonable but is
actually shifted so that H, K or L (or some combination of these) are off by
some integer value (often +/- 1).

DIALS uses the information written to the image headers to form its initial
model for the experiment geometry. That is, we trust that the beamline has been
set up correctly and the files produced at the beamline are a good
description of the experiment that was performed. Using the
:program:`dxtbx` library within :program:`cctbx` means we even recognise
specific beamlines in many cases, and interpret the files according to
custom format readers that are appropriate for those beamlines.

Unfortunately it is not always possible to be this trusting. Beamlines are
often changing and sometimes the metadata recorded with your diffraction images
are out of date. Sometimes values for things like the beam centre depend on
other factors such as the wavelength and detector distance, and your experiment
might have been performed outside the range that was used during beamline
calibration. Whatever the cause, incorrect image headers are a reality that
we have to be aware of. Some programs go as far as to ignore the image
headers entirely and pass the responsiblity on to the supposed existence of an
accurate 'site file' describing the experimental set up. Although it may be
easier for a user to edit an incorrect site file rather than broken image
headers, the problem still remains that if the experiment description is
sufficiently wrong, the processing will not be straightforward.

In this tutorial we describe some of the ways in which bad geometry
derived from wrong image headers can be worked around within DIALS, and we
use an example in which the initially incorrect geometry tricks the indexing
program into misindexing, which, if we were being careless, could have lead
to the integration of a useless data set. This tutorial is a cautionary tale,
the moral of which is that the user is still required to read the output of
the programs they run!

Tutorial data
-------------

The following example uses a dataset kindly provided by Wolfram Tempel, which
was collected at beamline 19-ID at the APS. This dataset is available for
download from |DPF3|.

.. |DPF3| image:: https://zenodo.org/badge/doi/10.5281/zenodo.45756.svg
          :target: http://dx.doi.org/10.5281/zenodo.45756

This dataset consists of a tar archive of bz2-compressed images. Very recent
versions of DIALS can read these directly, however here we shall be using
DIALS 1.1, as included in CCP4 7.0. In that case, we need to uncompress the
image files first. To do that, we first extract the archive::

  tar xvf DPF3_247398.tar

then uncompress each of the 200 images. This will probably take a couple of
minutes if the files are on a local hard drive::

  bunzip2 *.bz2

Import
^^^^^^

At this point we have no reason not to trust the image headers. We shall just
go ahead and import the whole sweep as normal::

  dials.import x247398/t1.0*.img

From the output we can see that beamline 19-ID at the APS is one that is
specifically recognised by DIALS. As expected we found a single sweep of
data. So all looks good so far.::

  --------------------------------------------------------------------------------
  DataBlock 0
    format: <class 'dxtbx.format.FormatSMVADSCSNAPSID19.FormatSMVADSCSNAPSID19'>
    num images: 200
    num sweeps: 1
    num stills: 0
  --------------------------------------------------------------------------------

We can get some more information about the expected experiment details using
:program:`dials.show`::

  dials.show datablock.json

which produces output including the experimental geometry::

  Detector:
  Panel:
    pixel_size:{0.1024,0.1024}
    image_size: {3072,3072}
    trusted_range: {1,65535}
    thickness: 0
    material:
    mu: 0
    fast_axis: {1,0,0}
    slow_axis: {0,-1,0}
    origin: {-159.98,154.501,-150}
  Max resolution: 1.355231

  Beam:
      wavelength: 1.28215
      sample to source direction : {0,0,1}
      divergence: 0
      sigma divergence: 0
      polarization normal: {0,1,0}
      polarization fraction: 0.999
  Beam centre: (159.98,154.50)

  Scan:
      image range:   {1,200}
      oscillation:   {-100,1}

  Goniometer:
      Rotation axis:   {-1,0,0}
      Fixed rotation:  {1,0,0,0,1,0,0,0,1}
      Setting rotation:{1,0,0,0,1,0,0,0,1}

At the moment we don't know that any of this is wrong. Happily, the 19-ID-specific
format has recognised the 'inverse phi' rotation of the goniometer at this
beamline, and thus produced a rotation axis of ``{-1,0,0}`` rather than
``{1,0,0}``. These inverse phi settings can the cause of problems with
processing data from currently unrecognised beamlines. As an aside, in such
a case we could force the rotation axis to be whatever we want like this::

  dials.import t1.0*.img geometry.goniometer.rotation_axis=-1,0,0

We can fix any aspect of the experimental geometry in this way, as long as we
know in advance what it should be. This information could all be included in
a file, say :file:`site.phil` and passed to :program:`dials.import` thus
combining the freedom of a site file with the ability to read image headers.
However, in general we would prefer to produce a new format in such cases.
More information about this is available in the :program:`dxtbx`
`paper <http://dx.doi.org/10.1107/S1600576714011996>`_

Find Spots
^^^^^^^^^^

Spot-finding in DIALS usually works well for Pilatus detectors, where default
assumptions about Poisson statistics of pixel counts, unity gain and no point
spread are accurate. These assumptions are not correct for CCD detectors and
this can be another source of problems with data processing. This may be the
subject of a future tutorial! In this case though, the defaults do a
reasonable, though possibly non-optimal job. We continue on regardless,
requesting only a larger number of processes to speed the job up::

  dials.find_spots datablock.json nproc=4


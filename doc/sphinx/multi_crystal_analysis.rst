Multi-crystal analysis with DIALS and BLEND
===========================================

Introduction
------------

BLEND_ is a CCP4 program for analysis of multiple data sets. It attempts to
identify isomorphous clusters that may be scaled and merged together to form a
more complete multi-crystal dataset. Clustering in blend is based on the refined
cell parameters from integration, so it is important that these are determined
accurately. Unfortunately, for narrow wedges of data (where BLEND is most
important) cell refinement may be complicated by issues such as the high
correlation between e.g. the detector distance and the cell volume.

.. _BLEND: http://www.ccp4.ac.uk/html/blend.html

One solution is to fix detector parameters during cell refinement so that the
detector is the same for every dataset processed. However, this is only feasible
if an accurate detector model is determined ahead of time, which might require
the use of a well-diffracting sacrificial crystal. If we only have the narrow
wedges of data available then it is more complicated to determine what the best
detector model would be to fix.

If we can make the assumption that the beam and detector *did not move* during
and between all of the datasets we collected then we can use DIALS joint
refinement to refine all of the crystal cells *at the same time*, simultaneously
with shared detector and beam models.

In this tutorial, we will attempt to do that for 73 sweeps of data collected
from crystals of TehA, a well-diffracting integral membrane protein measured
using *in situ* diffraction from crystallisation plates at room temperature.
Each sweep provides between 4 and 10 degrees of data.

Individual processing
---------------------

Before we perform the joint analysis, we want to do the individual analysis
to compare to. This will also give us intermediate files so that we don't have
to start from scratch when setting up the joint refinement job. Essentially
we just want to run a sequence of DIALS commands to process each recorded sweep.
This is an ideal job for scripting, which could be done in BASH, tcsh, perl,
ruby - whatever you feel most comfortable with. However here we will use Python,
or more specifically :program:`dials.python` because we will take advantage of
features in the cctbx to make it easy to write scripts that take advantage
of `parallel execution <http://cctbx.sourceforge.net/current/python/libtbx.easy_mp.html>`_.

.. highlight:: python
The script used to do this is reproduced here::

  #!/bin/env dials.python
  import os
  import sys
  import glob
  from libtbx import easy_run, easy_mp
  from dials.test import cd

  def process_sweep(task):
    """Process a single sweep of data. The parameter 'task' will be a
    tuple, the first element of which is an integer job number and the
    second is the path to the images to process"""

    num = task[0] + 1 # change first dataset number from 00 to 01
    template = task[1]

    # create directory
    with cd("sweep_%02d" % num):
      cmd = "dials.import {0}".format(template)
      easy_run.fully_buffered(command=cmd)
      easy_run.fully_buffered(command="dials.find_spots datablock.json")

      # initial indexing in P 1
      cmd = "dials.index datablock.json strong.pickle " +\
            "output.experiments=P1_experiments.json"
      easy_run.fully_buffered(command=cmd)
      if not os.path.isfile("P1_experiments.json"):
        print "Job %02d failed in initial indexing" % num
        return

      # bootstrap from the refined P 1 cell
      cmd = "dials.index P1_experiments.json strong.pickle space_group='H 3'"
      easy_run.fully_buffered(command=cmd)
      if not os.path.isfile("experiments.json"):
        print "Job %02d failed in indexing" % num
        return

      # static model refinement
      cmd = "dials.refine experiments.json indexed.pickle " + \
            "do_outlier_rejection=true"
      easy_run.fully_buffered(command=cmd)
      if not os.path.isfile("refined_experiments.json"):
        print "Job %02d failed in refinement" % num
        return

      # WARNING! Fast and dirty integration.
      # Do not use the result for scaling/merging!
      cmd = "dials.integrate refined_experiments.json indexed.pickle " + \
            "intensity.algorithm=sum prediction.dmin=3 prediction.dmax=8"
      easy_run.fully_buffered(command=cmd)
      if not os.path.isfile("integrated.pickle"):
        print "Job %02d failed during integration" % num
        return

      # create MTZ
      cmd = "dials.export_mtz refined_experiments.json integrated.pickle " +\
            "hklout=integrated.mtz"
      easy_run.fully_buffered(command=cmd)
      if not os.path.isfile("integrated.mtz"):
        print "Job %02d failed during MTZ export" % num
        return

    # if we got this far, return the path to the MTZ
    return "sweep_%02d/integrated.mtz" % num

  if __name__ == "__main__":

    if len(sys.argv) != 2:
      sys.exit("Usage: dials.python process_TehA.py /path/to/images")
    data_dir = sys.argv[1]

    pathname = os.path.join(data_dir, "*.log")
    logfiles = glob.glob(pathname)

    templates = [f[:-8] + "*.cbf" for f in logfiles]
    tasklist = list(enumerate(sorted(templates)))

    from libtbx import Auto
    nproc = easy_mp.get_processes(Auto)
    print nproc

    print "Attempting to process the following datasets, with {} processes".format(nproc)
    for task in tasklist:
      print "%d: %s" % task

    results = easy_mp.parallel_map(
      func=process_sweep,
      iterable=tasklist,
      processes=nproc,
      preserve_order=True)

    good_results = [e for e in results if e is not None]
    print "Successfully created the following MTZs:"
    for result in good_results:
      print result


We will now briefly describe what is in this script. The first lines are
just imports to bring in modules from the Python standard library as well as
:samp:`easy_run` and :samp:`easy_mp` from :samp:`libtbx` (part of cctbx) and
a class from the :samp:`dials.test` package that simplifies running commands in
a new directory. Following that is a definition for the function
:samp:`process_sweep` which will perform all the steps required to process one
dataset from images to unmerged MTZ. The code block under::

  if __name__ == "__main__":

are the lines that are executed when the script starts. First we check that the
script has been passed a path to images. Having looked at the directory
containing images we realised that each dataset is associated with a log file,
so a quick way to identify all the distinct datasets is just to list the
:file:`*.log` files in the data directory. As an alternative we could have run::

  dials.import /path/to/images/*.cbf

As this would have created a datablock listing all of the individual datasets
found, from which we could have extracted the ImageSweep templates. This would
have been a more general solution, but for this case the existence of the
:file:`.log` files gave us a simple alternative.

After manipulating the :file:`.log` filenames we have templates for each of the
datasets. We want to pass each of these into :samp:`process_sweep`, but instead
of doing this in serial we can use :samp:`easy_mp` to run in parallel. This will
be okay because inside :samp:`process_sweep`, we ensure that all results are
written into a new directory. First we use a facility of the :samp:`easy_mp`
module to determine the number of processes to run in parallel and then we submit
the job with :samp:`parallel_map`.

Within :samp:`process_sweep` all external commands are run within a :samp:`with`
block where execution is controlled by the *context manager* :samp:`cd`. If you
want the gory details, they are `here <https://docs.python.org/2/reference/datamodel.html#context-managers>`_.
Essentially this is a way to write clean code that tidies up after itself
properly. In this case, we will create a new directory, execute commands in that
directory, then change back to the old directory afterwards. If the directory
already exists, this will fail with an error.

The commands that are run inside the managed block are usual dials commands,
familiar from the earlier tutorial. There are a couple of interesting points
to note though. We know that the correct space group is *H* 3, but it turns out
that if we ask :program:`dials.index` to find an *H* 3 cell right from the start
then many of the sweeps fail to index. This is simply because the initial models
contained in :samp:`datablock.json` are too poor to locate a cell with the
symmetry constraints. However, for many of the sweeps the indexing program will
refine the *P* 1 solution to the correct cell. For this reason we first run
indexing in *P* 1::

  dials.index datablock.json strong.pickle output.experiments=P1_experiments.json

and then we feed the refined :file:`P1_experiments.json` back into
:program:`dials.index` specifying the correct symmetry::

  dials.index P1_experiments.json strong.pickle space_group='H 3'

When :program:`dials.index` is passed an :file:`experiments.json` containing
a crystal model rather than just a :file:`databock.json` then it automatically
uses a :samp:`known_orientation` indexer, which avoids doing the basis vector
search again. It uses the basis of the refined *P* 1 cell and just reassigns
indices under the assumption of *H* 3 symmetry. The symmetry constraints are
then enforced during the refinement steps carried out by :program:`dials.index`.
This procedure gives us a greater success rate of indexing in *H* 3, and required
no manual intervention.

Following indexing we do scan-static cell refinement::

  dials.refine experiments.json indexed.pickle do_outlier_rejection=true

Outlier rejection was switched on in an attempt to avoid any zingers or other
errant spots from affecting our refined cells. Without analysing the data closer
it is not clear whether there are any particularly bad outliers here. We could repeat
the whole analysis with this switched off if we want to investigate more closely,
or look through all the :file:`dials.refine.log` files to see results of the
outlier rejection step.

We don't bother with the time-consuming step of scan-varying refinement, because
it is the scan-static cell that will be written into the MTZ header. Scan-
varying refinement would give us better models for integration but as we will
only be running blend in 'analysis' mode we are in the unusual situation of not
actually caring what the intensities are. In this case, the MTZ file is just a
carrier for the globally refined unit cell!

Following refinement we integrate the data in a very quick and dirty way, simply
to get an MTZ file as quickly as possible. This is a terrible way to integrate
data usually!::

  dials.integrate refined_experiments.json indexed.pickle intensity.algorithm=sum prediction.dmin=3 prediction.dmax=8

The :samp:`intensity.algorithm=sum` option ensures we only do summation integration,
no profile fitting, while the :samp:`prediction.dmin=3` and
:samp:`prediction.dmax=8` options only integrate data between 3 and 8 Angstroms.

Finally we use :program:`dials.export` to create an MTZ file::

  dials.export_mtz refined_experiments.json integrated.pickle hklout=integrated.mtz

After each of these major steps we check whether the last command ran successfully
by checking for the existence of an expected output file. If the file does not
exist we make no effort to rescue that dataset, we just return early from the
:samp:`process_sweep` function, freeing up a process so that
:samp:`parallel_map` can start up the next.

We saved this as :samp:`process_TehA.py` and then ran it as follows::

  time dials.python process_TehA.py /path/to/images/

On a Linux desktop with a Core i7 CPU this script took XXXX minutes to run
and successfully processed 41 datasets. Here is the output::

  Attempting to process the following datasets, with 7 processes
  0: /home/david/xray/TehA/xta30_1_*.cbf
  1: /home/david/xray/TehA/xta31_1_*.cbf
  2: /home/david/xray/TehA/xta32_1_*.cbf
  3: /home/david/xray/TehA/xta33_1_*.cbf
  4: /home/david/xray/TehA/xta34_1_*.cbf
  5: /home/david/xray/TehA/xta9_1_*.cbf
  6: /home/david/xray/TehA/xta9_2_*.cbf
  7: /home/david/xray/TehA/xtal10_1_*.cbf
  8: /home/david/xray/TehA/xtal11_1_*.cbf
  9: /home/david/xray/TehA/xtal12_1_*.cbf
  10: /home/david/xray/TehA/xtal12_2_*.cbf
  11: /home/david/xray/TehA/xtal13_1_*.cbf
  12: /home/david/xray/TehA/xtal14_1_*.cbf
  13: /home/david/xray/TehA/xtal15_1_*.cbf
  14: /home/david/xray/TehA/xtal16_1_*.cbf
  15: /home/david/xray/TehA/xtal17_1_*.cbf
  16: /home/david/xray/TehA/xtal18_1_*.cbf
  17: /home/david/xray/TehA/xtal19_1_*.cbf
  18: /home/david/xray/TehA/xtal1_1_*.cbf
  19: /home/david/xray/TehA/xtal20_1_*.cbf
  20: /home/david/xray/TehA/xtal21_1_*.cbf
  21: /home/david/xray/TehA/xtal22_1_*.cbf
  22: /home/david/xray/TehA/xtal23_1_*.cbf
  23: /home/david/xray/TehA/xtal24_1_*.cbf
  24: /home/david/xray/TehA/xtal25_1_*.cbf
  25: /home/david/xray/TehA/xtal26_1_*.cbf
  26: /home/david/xray/TehA/xtal26_2_*.cbf
  27: /home/david/xray/TehA/xtal27_1_*.cbf
  28: /home/david/xray/TehA/xtal28_1_*.cbf
  29: /home/david/xray/TehA/xtal29_1_*.cbf
  30: /home/david/xray/TehA/xtal2_1_*.cbf
  31: /home/david/xray/TehA/xtal35_1_*.cbf
  32: /home/david/xray/TehA/xtal36_1_*.cbf
  33: /home/david/xray/TehA/xtal37_1_*.cbf
  34: /home/david/xray/TehA/xtal37_2_*.cbf
  35: /home/david/xray/TehA/xtal38_1_*.cbf
  36: /home/david/xray/TehA/xtal39_1_*.cbf
  37: /home/david/xray/TehA/xtal3_2_*.cbf
  38: /home/david/xray/TehA/xtal40_1_*.cbf
  39: /home/david/xray/TehA/xtal40_2_*.cbf
  40: /home/david/xray/TehA/xtal40_3_*.cbf
  41: /home/david/xray/TehA/xtal40_4_*.cbf
  42: /home/david/xray/TehA/xtal41_1_*.cbf
  43: /home/david/xray/TehA/xtal42_1_*.cbf
  44: /home/david/xray/TehA/xtal43_1_*.cbf
  45: /home/david/xray/TehA/xtal44_1_*.cbf
  46: /home/david/xray/TehA/xtal45_1_*.cbf
  47: /home/david/xray/TehA/xtal46_1_*.cbf
  48: /home/david/xray/TehA/xtal47_1_*.cbf
  49: /home/david/xray/TehA/xtal48_1_*.cbf
  50: /home/david/xray/TehA/xtal49_1_*.cbf
  51: /home/david/xray/TehA/xtal4_3_*.cbf
  52: /home/david/xray/TehA/xtal50_1_*.cbf
  53: /home/david/xray/TehA/xtal50_2_*.cbf
  54: /home/david/xray/TehA/xtal51_1_*.cbf
  55: /home/david/xray/TehA/xtal52_1_*.cbf
  56: /home/david/xray/TehA/xtal53_1_*.cbf
  57: /home/david/xray/TehA/xtal54_1_*.cbf
  58: /home/david/xray/TehA/xtal55_1_*.cbf
  59: /home/david/xray/TehA/xtal55_2_*.cbf
  60: /home/david/xray/TehA/xtal56_1_*.cbf
  61: /home/david/xray/TehA/xtal56_2_*.cbf
  62: /home/david/xray/TehA/xtal57_1_*.cbf
  63: /home/david/xray/TehA/xtal58_1_*.cbf
  64: /home/david/xray/TehA/xtal58_2_*.cbf
  65: /home/david/xray/TehA/xtal58_3_*.cbf
  66: /home/david/xray/TehA/xtal59_1_*.cbf
  67: /home/david/xray/TehA/xtal5_1_*.cbf
  68: /home/david/xray/TehA/xtal60_1_*.cbf
  69: /home/david/xray/TehA/xtal60_2_*.cbf
  70: /home/david/xray/TehA/xtal6_1_*.cbf
  71: /home/david/xray/TehA/xtal7_1_*.cbf
  72: /home/david/xray/TehA/xtal8_1_*.cbf
  Job 07 failed in initial indexing
  Job 05 failed in indexing
  Job 08 failed in indexing
  Job 09 failed in indexing
  Job 12 failed in indexing
  Job 13 failed in indexing
  Job 11 failed in indexing
  Job 16 failed in initial indexing
  Job 22 failed in initial indexing
  Job 21 failed in initial indexing
  Job 33 failed in initial indexing
  Job 38 failed in indexing
  Job 36 failed in indexing
  Job 39 failed in indexing
  Job 40 failed in indexing
  Job 41 failed in indexing
  Job 42 failed in indexing
  Job 45 failed in indexing
  Job 46 failed in indexing
  Job 53 failed in initial indexing
  Job 48 failed in indexing
  Job 50 failed in indexing
  Job 56 failed in initial indexing
  Job 58 failed in initial indexing
  Job 62 failed in indexing
  Job 63 failed in indexing
  Job 71 failed in indexing
  Job 67 failed in indexing
  Job 70 failed in indexing
  Job 69 failed in indexing
  Job 72 failed in initial indexing
  Job 73 failed in indexing
  Successfully created the following MTZs:
  sweep_01/integrated.mtz
  sweep_02/integrated.mtz
  sweep_03/integrated.mtz
  sweep_04/integrated.mtz
  sweep_06/integrated.mtz
  sweep_10/integrated.mtz
  sweep_14/integrated.mtz
  sweep_15/integrated.mtz
  sweep_17/integrated.mtz
  sweep_18/integrated.mtz
  sweep_19/integrated.mtz
  sweep_20/integrated.mtz
  sweep_23/integrated.mtz
  sweep_24/integrated.mtz
  sweep_25/integrated.mtz
  sweep_26/integrated.mtz
  sweep_27/integrated.mtz
  sweep_28/integrated.mtz
  sweep_29/integrated.mtz
  sweep_30/integrated.mtz
  sweep_31/integrated.mtz
  sweep_32/integrated.mtz
  sweep_34/integrated.mtz
  sweep_35/integrated.mtz
  sweep_37/integrated.mtz
  sweep_43/integrated.mtz
  sweep_44/integrated.mtz
  sweep_47/integrated.mtz
  sweep_49/integrated.mtz
  sweep_51/integrated.mtz
  sweep_52/integrated.mtz
  sweep_54/integrated.mtz
  sweep_55/integrated.mtz
  sweep_57/integrated.mtz
  sweep_59/integrated.mtz
  sweep_60/integrated.mtz
  sweep_61/integrated.mtz
  sweep_64/integrated.mtz
  sweep_65/integrated.mtz
  sweep_66/integrated.mtz
  sweep_68/integrated.mtz

  real	8m31.718s
  user	21m49.950s
  sys	1m46.923s

Check the unit cells

  dials.show_models sweep_*/refined_experiments.json | grep "Unit cell"

We see one dataset has a hugely outlying cell. Let's remove that.

Acknowledgements
----------------

Danny Axford, Nien-Jen Hu, James Foadi, Hassanul Ghani Choudhury, So Iwata, Konstantinos Beis, Gwyndaf Evans & Yilmaz Alguel

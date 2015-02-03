Multi-crystal analysis with DIALS and BLEND
===========================================

Introduction
------------

BLEND_ is a CCP4 program for analysis of multiple data sets. It attempts to
identify isomorphous clusters that may be scaled and merged together to form a
more complete multi-crystal dataset. Clustering in blend is based on the refined
cell parameters from integration, so it is important that these are determined
accurately. Unfortunately, for narrow wedges of data (where BLEND is most
useful) cell refinement may be complicated by issues such as the high
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
The script we used to do this is reproduced below. You can copy this into a file,
save it as :samp:`process_TehA.py` and then run it as follows::

  time dials.python process_TehA.py /path/to/images/

On a Linux desktop with a Core i7 CPU running at 3.07GHz it script took about 8
minutes to run (though file i/o is a significant factor)
and successfully processed 41 datasets. If time is short, you
might like to start running it now before reading the description of what the
script does. If time is *really* short then try uncommenting the line
:samp:`tasklist = tasklist[0:35]` to reduce the number of datasets processed.::

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

    num = task[0]
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
            "do_outlier_rejection=true use_all_reflections=true"
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

    if len(tasklist) == 0: sys.exit("No images found!")

    # uncomment the following line if short on time!
    #tasklist = tasklist[0:35]

    from libtbx import Auto
    nproc = easy_mp.get_processes(Auto)

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

We will now describe what is in this script. The first lines are
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
search again. It uses the basis of the refined *P* 1 cell and just assigns
indices under the assumption of *H* 3 symmetry. The symmetry constraints are
then enforced during the refinement steps carried out by :program:`dials.index`.
This procedure gives us a greater success rate of indexing in *H* 3, and required
no manual intervention.

Following indexing we do scan-static cell refinement::

  dials.refine experiments.json indexed.pickle do_outlier_rejection=true use_all_reflections=true

Outlier rejection was switched on in an attempt to avoid any zingers or other
errant spots from affecting our refined cells. Without analysing the data closer
it is not clear whether there are any particularly bad outliers here. We could repeat
the whole analysis with this switched off if we want to investigate more closely,
or look through all the :file:`dials.refine.log` files to see results of the
outlier rejection step.

We elected use all reflections rather than taking a random subset because these
are narrow wedges and there are few reflections anyway. Taking a random subset
is only a time-saving procedure, and it won't provide much benefit here anyway.

We don't bother with the time-consuming step of scan-varying refinement, because
it is the scan-static cell that will be written into the MTZ header. Scan-
varying refinement would give us better models for integration but as we will
only be running blend in 'analysis' mode we are in the unusual situation of not
actually caring what the intensities are. In this case, the MTZ file is just a
carrier for the globally refined unit cell!

Following refinement we integrate the data in a very quick and dirty way, simply
to get an MTZ file as fast as possible. This is a terrible way to integrate
data usually!::

  dials.integrate refined_experiments.json indexed.pickle intensity.algorithm=sum prediction.dmin=3 prediction.dmax=8

The :samp:`intensity.algorithm=sum` option ensures we only do summation integration,
no profile fitting, while the :samp:`prediction.dmin=3` and
:samp:`prediction.dmax=8` options only integrate data between 3 and 8 Angstroms.

.. warning::

  Do not use the data produced by this script for scaling and merging. More
  careful processing should be done first!

Finally we use :program:`dials.export` to create an MTZ file::

  dials.export_mtz refined_experiments.json integrated.pickle hklout=integrated.mtz

After each of these major steps we check whether the last command ran successfully
by checking for the existence of an expected output file. If the file does not
exist we make no effort to rescue the dataset, we just return early from the
:samp:`process_sweep` function, freeing up a process so that
:samp:`parallel_map` can start up the next.

Here is the output of a run of the script::

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
  Job 06 failed in initial indexing
  Job 04 failed in indexing
  Job 07 failed in indexing
  Job 08 failed in indexing
  Job 12 failed in indexing
  Job 11 failed in indexing
  Job 10 failed in indexing
  Job 15 failed in initial indexing
  Job 20 failed in initial indexing
  Job 21 failed in initial indexing
  Job 32 failed in initial indexing
  Job 37 failed in indexing
  Job 35 failed in indexing
  Job 38 failed in indexing
  Job 39 failed in indexing
  Job 40 failed in indexing
  Job 41 failed in indexing
  Job 44 failed in indexing
  Job 45 failed in indexing
  Job 47 failed in indexing
  Job 52 failed in initial indexing
  Job 49 failed in indexing
  Job 55 failed in initial indexing
  Job 57 failed in initial indexing
  Job 61 failed in indexing
  Job 62 failed in indexing
  Job 69 failed in indexing
  Job 66 failed in indexing
  Job 68 failed in indexing
  Job 70 failed in indexing
  Job 71 failed in initial indexing
  Job 72 failed in indexing
  Successfully created the following MTZs:
  sweep_00/integrated.mtz
  sweep_01/integrated.mtz
  sweep_02/integrated.mtz
  sweep_03/integrated.mtz
  sweep_05/integrated.mtz
  sweep_09/integrated.mtz
  sweep_13/integrated.mtz
  sweep_14/integrated.mtz
  sweep_16/integrated.mtz
  sweep_17/integrated.mtz
  sweep_18/integrated.mtz
  sweep_19/integrated.mtz
  sweep_22/integrated.mtz
  sweep_23/integrated.mtz
  sweep_24/integrated.mtz
  sweep_25/integrated.mtz
  sweep_26/integrated.mtz
  sweep_27/integrated.mtz
  sweep_28/integrated.mtz
  sweep_29/integrated.mtz
  sweep_30/integrated.mtz
  sweep_31/integrated.mtz
  sweep_33/integrated.mtz
  sweep_34/integrated.mtz
  sweep_36/integrated.mtz
  sweep_42/integrated.mtz
  sweep_43/integrated.mtz
  sweep_46/integrated.mtz
  sweep_48/integrated.mtz
  sweep_50/integrated.mtz
  sweep_51/integrated.mtz
  sweep_53/integrated.mtz
  sweep_54/integrated.mtz
  sweep_56/integrated.mtz
  sweep_58/integrated.mtz
  sweep_59/integrated.mtz
  sweep_60/integrated.mtz
  sweep_63/integrated.mtz
  sweep_64/integrated.mtz
  sweep_65/integrated.mtz
  sweep_67/integrated.mtz

  real	7m46.071s
  user	22m19.016s
  sys	1m47.299s

Analysis of individually processed datasets
-------------------------------------------

The paths to :file:`integrated.mtz` files can be copied directly into a file,
say :file:`individual_mtzs.dat`, and passed to blend for analysis::

  echo "END" | blend -a individual_mtzs.dat

The dendrogram resulting from clustering is shown here:

  .. image:: figures/tree_01.png

Immediately the dendrogram shows that datasets 7 and 28 are extreme outliers.
From :file:`FINAL_list_of_files.dat` we can see that these refer to
:file:`sweep_13/integrated.mtz` and :file:`sweep_46/integrated.mtz`.
As we kept all the dials :file:`.log` files
from DIALS processing we could investigate this further, however as these are
only two sweeps out of 41, our time is better spent throwing them away and
moving on. So, edit :file:`individual_mtzs.dat` to remove
the lines :file:`sweep_13/integrated.mtz` and :file:`sweep_46/integrated.mtz`
and rerun blend.

Now the dendrogram looks better:

  .. image:: figures/tree_02.png

The Linear Cell Variation (LCV) is now less than 1%, with an absolute value
of 0.42 Angstroms, indicating good isomorphism amongst all the remaining
datasets.

Joint refinement
----------------

Now that we have done the BLEND analysis for individually processed datasets,
we would like to do joint refinement of the crystals to reduce correlations
between the detector or beam parameters with individual crystals. As motivation
we may look at these correlations for one of these datasets. For example::

  cd sweep_00
  dials.refine experiments.json indexed.pickle \
    track_parameter_correlation=true correlation_plot.filename=corrplot.png
  cd ..

The new file :file:`sweep_00/corrplot.png` shows correlations between parameters
refined with this single 8 degree dataset. Clearly parameters like the
detector distance and the crystal metrical matrix parameters are highly
correlated.

 .. image:: figures/sweep_00_corrplot.png

Although the DIALS toolkit has a sophisticated mechanism for modelling
multi-experiment data, the user interface for handling such data is still
rather limited. In order to do joint refinement of the sweeps we need to combine them
into a single multi-experiment :file:`experiments.json` and corresponding
:file:`reflections.pickle`. Whilst doing this we want to reduce the separate
detector, beam and goniometer models for each experiment into a single shared
model of each type. The program :program:`dials.combine_experiments` can
be used for this, but first we have to prepare an input file with a text editor
listing the individual sweeps in order. We can use
:file:`individual_mtzs.dat` as a template to start with. In our case the final
file looks like this::

  input {
    experiments = "/home/david/TehA_processing/sweep_00/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_01/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_02/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_03/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_05/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_09/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_13/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_14/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_16/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_17/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_18/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_19/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_22/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_23/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_24/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_25/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_26/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_27/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_28/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_29/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_30/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_31/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_33/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_34/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_36/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_42/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_43/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_46/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_48/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_50/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_51/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_53/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_54/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_56/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_58/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_59/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_60/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_63/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_64/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_65/refined_experiments.json"
    experiments = "/home/david/TehA_processing/sweep_67/refined_experiments.json"
    reflections = "/home/david/TehA_processing/sweep_00/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_01/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_02/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_03/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_05/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_09/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_13/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_14/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_16/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_17/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_18/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_19/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_22/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_23/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_24/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_25/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_26/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_27/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_28/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_29/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_30/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_31/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_33/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_34/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_36/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_42/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_43/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_46/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_48/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_50/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_51/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_53/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_54/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_56/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_58/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_59/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_60/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_63/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_64/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_65/indexed.pickle"
    reflections = "/home/david/TehA_processing/sweep_67/indexed.pickle"
  }

We called this file :file:`experiments_and_reflections.phil` then run
:program:`dials.combine_experiments` like this::

  dials.combine_experiments experiments_and_reflections.phil \
    reference_from_experiment.beam=0 \
    reference_from_experiment.goniometer=0 \
    reference_from_experiment.detector=0

The :samp:`reference_from_experiment` options tell the program to replace all
beam, goniometer and detector models in the input experiments with those
models taken from the first experiment, i.e. experiment '0' using 0-based
indexing. The output lists the number of reflections in each sweep contributing
to the final :file:`combined_reflections.pickle`::

  ---------------------
  | Experiment | Nref |
  ---------------------
  | 0          | 1446 |
  | 1          | 1422 |
  | 2          | 1209 |
  | 3          | 1376 |
  | 4          | 452  |
  | 5          | 1663 |
  | 6          | 1528 |
  | 7          | 1445 |
  | 8          | 1275 |
  | 9          | 239  |
  | 10         | 1614 |
  | 11         | 1052 |
  | 12         | 1845 |
  | 13         | 1495 |
  | 14         | 2041 |
  | 15         | 1308 |
  | 16         | 1839 |
  | 17         | 1828 |
  | 18         | 1644 |
  | 19         | 243  |
  | 20         | 1061 |
  | 21         | 2416 |
  | 22         | 1884 |
  | 23         | 949  |
  | 24         | 3569 |
  | 25         | 2967 |
  | 26         | 935  |
  | 27         | 1329 |
  | 28         | 650  |
  | 29         | 1324 |
  | 30         | 633  |
  | 31         | 1231 |
  | 32         | 2131 |
  | 33         | 2094 |
  | 34         | 2141 |
  | 35         | 1661 |
  | 36         | 2543 |
  | 37         | 2227 |
  | 38         | 1138 |
  ---------------------
  Saving combined experiments to combined_experiments.json
  Saving combined reflections to combined_reflections.pickle

We may also inspect the contents of :file:`combined_experiments.json`, by using
:program:`dials.show_models`, for example::

  dials.show_models combined_experiments.json

Useful though this is, it is clear how this could become unwieldy as the number
of experiments increases. Work on better interfaces to multi-crystal (or
generally, multi-experiment) data is ongoing within the DIALS project.
Suggestions are always welcome!

Now we have the joint experiments and reflections files we can run our multi-
crystal refinement job. First we try outlier rejection, so that the refinement
run is similar to the jobs we ran on individual datasets::

  dials.refine combined_experiments.json combined_reflections.pickle \
    do_outlier_rejection=true

::

  The following parameters have been modified:

  refinement {
    reflections {
      do_outlier_rejection = true
    }
  }
  input {
    experiments = combined_experiments.json
    reflections = combined_reflections.pickle
  }

  Configuring refiner

  Summary statistics for observations matched to predictions:
  ----------------------------------------------------------------------
  |                   | Min    | Q1      | Med        | Q3     | Max   |
  ----------------------------------------------------------------------
  | Xc - Xo (mm)      | -14.61 | -0.8011 | -0.08364   | 0.7517 | 15.76 |
  | Yc - Yo (mm)      | -21.55 | -0.4907 | -0.01917   | 0.4474 | 16.99 |
  | Phic - Phio (deg) | -16.99 | -0.2279 | -0.0006402 | 0.2305 | 28.72 |
  | X weights         | 108.4  | 129.6   | 132.2      | 133.8  | 135.2 |
  | Y weights         | 114.8  | 133.8   | 134.7      | 135.1  | 135.2 |
  | Phi weights       | 81.19  | 99.99   | 100        | 100    | 100   |
  ----------------------------------------------------------------------

  15921 reflections have been rejected as outliers
  Traceback (most recent call last):
    File "/home/david/bsx/cctbx-svn/build/../sources/dials/command_line/refine.py", line 370, in <module>
      halraiser(e)
    File "/home/david/bsx/cctbx-svn/build/../sources/dials/command_line/refine.py", line 368, in <module>
      script.run()
    File "/home/david/bsx/cctbx-svn/build/../sources/dials/command_line/refine.py", line 274, in run
      reflections, experiments)
    File "/home/david/bsx/cctbx-svn/sources/dials/algorithms/refinement/refiner.py", line 336, in from_parameters_data_experiments
      verbosity=verbosity)
    File "/home/david/bsx/cctbx-svn/sources/dials/algorithms/refinement/refiner.py", line 581, in _build_components
      target = cls.config_target(params, experiments, refman, pred_param, do_stills)
    File "/home/david/bsx/cctbx-svn/sources/dials/algorithms/refinement/refiner.py", line 1004, in config_target
      options.jacobian_max_nref)
    File "/home/david/bsx/cctbx-svn/sources/dials/algorithms/refinement/target.py", line 404, in __init__
      self._reflection_manager.finalise()
    File "/home/david/bsx/cctbx-svn/sources/dials/algorithms/refinement/reflection_manager.py", line 237, in finalise
      self._check_too_few()
    File "/home/david/bsx/cctbx-svn/sources/dials/algorithms/refinement/reflection_manager.py", line 262, in _check_too_few
      raise RuntimeError(msg)
  RuntimeError: Please report this error to dials-support@lists.sourceforge.net: Remaining number of reflections = 6, for experiment 19, which is below the configured limit for this reflection manager

Oops! That wasn't good. Looking at the error we see that experiment 19 provides
only 6 reflections to refinement, which is disallowed by a default
parameters of :program:`dials.refine`, namely `minimum_number_of_reflections=20`.
But from the output of :program:`dials.combine_experiments` we see that experiment
19 has 243 indexed reflections. What happened? Well, forcing the individual
experiments to share the beam and detector models of experiment 0 has led to some
very poor predictions for some of these experiments. See the ``Summary statistics``
table, where the worst positional residuals are greater than 20 mm! We may put this
down to the very narrow wedges of data we have. Experiment 19 is one of the
narrowest, with only 4 degrees of data. Outlier rejection is not a good idea here
because it selectively removes reflections from the worst fitting experiments.

Instead we try without outlier rejection::

  dials.refine combined_experiments.json combined_reflections.pickle \
    use_all_reflections=true \
    output.experiments=refined_combined_experiments.json

This worked much better::

  The following parameters have been modified:

  output {
    experiments = refined_combined_experiments.json
  }
  refinement {
    reflections {
      use_all_reflections = true
    }
  }
  input {
    experiments = combined_experiments.json
    reflections = combined_reflections.pickle
  }

  Configuring refiner

  Summary statistics for observations matched to predictions:
  ----------------------------------------------------------------------
  |                   | Min    | Q1      | Med        | Q3     | Max   |
  ----------------------------------------------------------------------
  | Xc - Xo (mm)      | -14.61 | -0.8011 | -0.08364   | 0.7517 | 15.76 |
  | Yc - Yo (mm)      | -21.55 | -0.4907 | -0.01917   | 0.4474 | 16.99 |
  | Phic - Phio (deg) | -16.99 | -0.2279 | -0.0006402 | 0.2305 | 28.72 |
  | X weights         | 108.4  | 129.6   | 132.2      | 133.8  | 135.2 |
  | Y weights         | 114.8  | 133.8   | 134.7      | 135.1  | 135.2 |
  | Phi weights       | 81.19  | 99.99   | 100        | 100    | 100   |
  ----------------------------------------------------------------------

  Performing refinement...

  Refinement steps:
  -----------------------------------------------
  | Step | Nref  | RMSD_X  | RMSD_Y  | RMSD_Phi |
  |      |       | (mm)    | (mm)    | (deg)    |
  -----------------------------------------------
  | 0    | 56703 | 1.6811  | 1.3938  | 1.3119   |
  | 1    | 56703 | 1.3728  | 1.0393  | 0.70978  |
  | 2    | 56703 | 1.1418  | 0.86757 | 0.65172  |
  | 3    | 56703 | 0.87359 | 0.66465 | 0.57709  |
  | 4    | 56703 | 0.60635 | 0.47194 | 0.44672  |
  | 5    | 56703 | 0.37995 | 0.31262 | 0.28325  |
  | 6    | 56703 | 0.22145 | 0.19743 | 0.16597  |
  | 7    | 56703 | 0.17484 | 0.16522 | 0.12868  |
  | 8    | 56703 | 0.17164 | 0.16306 | 0.12515  |
  | 9    | 56703 | 0.1714  | 0.16287 | 0.12503  |
  | 10   | 56703 | 0.1713  | 0.16277 | 0.12496  |
  | 11   | 56703 | 0.17131 | 0.16274 | 0.12491  |
  | 12   | 56703 | 0.17132 | 0.16273 | 0.12489  |
  | 13   | 56703 | 0.17132 | 0.16273 | 0.12489  |
  -----------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 1374 | 0.63135 | 0.40973 | 0.35223  |
  | 1   | 1326 | 0.65259 | 0.39367 | 0.34253  |
  | 2   | 1138 | 0.90566 | 0.85055 | 0.75363  |
  | 3   | 1294 | 0.67156 | 0.5088  | 0.27957  |
  | 4   | 406  | 0.76238 | 0.50361 | 0.3676   |
  | 5   | 1578 | 1.0475  | 1.5447  | 0.93663  |
  | 6   | 1452 | 0.64011 | 0.33055 | 0.34482  |
  | 7   | 1372 | 1.0639  | 1.116   | 0.89393  |
  | 8   | 1203 | 1.0557  | 1.4787  | 0.6994   |
  | 9   | 213  | 2.0415  | 2.0383  | 1.3647   |
  | 10  | 1543 | 0.7825  | 0.47977 | 0.5151   |
  | 11  | 980  | 0.96061 | 1.1603  | 0.72562  |
  | 12  | 1783 | 0.74111 | 0.84793 | 0.67643  |
  | 13  | 1424 | 0.73923 | 0.51892 | 0.37183  |
  | 14  | 1937 | 1.1602  | 1.4408  | 0.84359  |
  | 15  | 1237 | 0.92553 | 0.50867 | 0.42323  |
  | 16  | 1751 | 0.71129 | 0.37352 | 0.34289  |
  | 17  | 1742 | 0.66178 | 0.40449 | 0.29842  |
  | 18  | 1550 | 0.84153 | 1.2567  | 0.71992  |
  | 19  | 222  | 1.1245  | 0.77295 | 0.95415  |
  ---------------------------------------------
  Table truncated to show the first 20 experiments only
  Re-run with verbosity >= 2 to show all experiments
  Saving refined experiments to refined_combined_experiments.json

The overall final RMSDs are 0.17 mm in X, 0.16 mm in Y and 0.12 degrees in
:math:`\phi`. The RMSDs per experiment are also shown, but only for the first
20 experiments. Rerunning with :samp:`verbosity=2` does give the full table,
but also produces a great deal more log output, so it would be easier to find
in the file :file:`dials.refine.log` rather than scrolling up pages in your
terminal.

We can compare the RMSDs from individually refined experiments to those from
the joint experiments. For example, look at the RSMDs for experiment 0, in the
logfile :file:`sweep_00/dials.refine.log`::

  RMSDs by experiment:
  --------------------------------------------
  | Exp | Nref | RMSD_X | RMSD_Y  | RMSD_Z   |
  |     |      | (px)   | (px)    | (images) |
  --------------------------------------------
  | 0   | 1342 | 0.534  | 0.30643 | 0.2561   |
  --------------------------------------------

Clearly allowing the detector and beam to refine only against this data lets
the model better fit the observations, but is it a more accurate description of
reality? Given that we *know* or can comfortably assume that the detector and
beam did not move between data collections, then the constraints applied by
joint refinement seem appropriate. For better parity with the original results
perhaps we should use outlier rejection though. Now the models are close enough
it is safe to do so::

  dials.refine refined_combined_experiments.json combined_reflections.pickle \
    use_all_reflections=true \
    do_outlier_rejection=true \
    output.experiments=refined_combined_experiments_outrej.json

The RMSD tables resulting from this::

  Refinement steps:
  ------------------------------------------------
  | Step | Nref  | RMSD_X  | RMSD_Y   | RMSD_Phi |
  |      |       | (mm)    | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 50112 | 0.10315 | 0.062074 | 0.058395 |
  | 1    | 50112 | 0.10292 | 0.061742 | 0.057896 |
  | 2    | 50112 | 0.10271 | 0.061592 | 0.057869 |
  | 3    | 50112 | 0.1024  | 0.061383 | 0.057734 |
  | 4    | 50112 | 0.10213 | 0.061227 | 0.057411 |
  | 5    | 50112 | 0.10197 | 0.061185 | 0.057029 |
  | 6    | 50112 | 0.10186 | 0.061202 | 0.056831 |
  | 7    | 50112 | 0.10178 | 0.061214 | 0.056807 |
  | 8    | 50112 | 0.10173 | 0.061164 | 0.056806 |
  | 9    | 50112 | 0.10168 | 0.061055 | 0.056777 |
  | 10   | 50112 | 0.10167 | 0.060948 | 0.056713 |
  | 11   | 50112 | 0.1017  | 0.060897 | 0.05664  |
  | 12   | 50112 | 0.10172 | 0.060884 | 0.056602 |
  | 13   | 50112 | 0.10172 | 0.060882 | 0.056594 |
  | 14   | 50112 | 0.10172 | 0.060882 | 0.056593 |
  ------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 1302 | 0.57135 | 0.34799 | 0.30443  |
  | 1   | 1275 | 0.59907 | 0.34379 | 0.31076  |
  | 2   | 1008 | 0.68104 | 0.4229  | 0.29659  |
  | 3   | 1213 | 0.61056 | 0.4238  | 0.27042  |
  | 4   | 373  | 0.6637  | 0.41751 | 0.28468  |
  | 5   | 1425 | 0.53209 | 0.30844 | 0.25475  |
  | 6   | 1426 | 0.51294 | 0.28226 | 0.23702  |
  | 7   | 1236 | 0.65703 | 0.32861 | 0.27816  |
  | 8   | 1091 | 0.54379 | 0.34609 | 0.25901  |
  | 9   | 137  | 1.2479  | 0.48073 | 0.31642  |
  | 10  | 1484 | 0.5417  | 0.33476 | 0.2514   |
  | 11  | 906  | 0.56075 | 0.39302 | 0.26312  |
  | 12  | 1697 | 0.53371 | 0.33843 | 0.25628  |
  | 13  | 1353 | 0.59367 | 0.32434 | 0.27128  |
  | 14  | 1765 | 0.55622 | 0.30903 | 0.25697  |
  | 15  | 1101 | 0.67655 | 0.35542 | 0.31188  |
  | 16  | 1633 | 0.56375 | 0.32634 | 0.30048  |
  | 17  | 1654 | 0.53093 | 0.3281  | 0.26622  |
  | 18  | 1401 | 0.51477 | 0.37377 | 0.27729  |
  | 19  | 171  | 0.89704 | 0.38654 | 0.39885  |
  ---------------------------------------------
  Table truncated to show the first 20 experiments only
  Re-run with verbosity >= 2 to show all experiments
  Saving refined experiments to refined_combined_experiments_outrej.json

Now we have RMSDs in X down to 0.1 mm, in Y to 0.06 mm and 0.06 degrees in
:math:`\phi`. The RMSDs for experiment 0 are not so much worse than from the
individual refinement job. We are happy with this result and move on to
re-integrating the data to create MTZs for BLEND.

Analysis of jointly refined datasets
------------------------------------

:program:`dials.integrate` will not work with our :file:`refined_combined_experiments_outrej.json`
and :file:`combined_reflections.pickle` directly, so we have to separate these
into individual files for each experiment. It is best to do this inside a new
directory::

  mkdir joint
  cd !$
  dials.split_experiments ../refined_combined_experiments_outrej.json ../combined_reflections.pickle

This fills the directory with 39 individual :file:`experiments_##.json` and
:file:`reflections_##.pickle` files. To integrate these quickly we want a script
to run in parallel, similar to the one used previously::

  #!/bin/env dials.python
  import os
  import sys
  import glob
  from libtbx import easy_run, easy_mp
  from dials.test import cd

  def process_sweep(task):
    """Process a single sweep of data. The parameter 'task' will be a
    tuple, the first element of which is an integer job number and the
    second is the path to the directory containing the data"""

    num = task[0]
    datadir = task[1]

    experiments_file = "experiments_%02d.json" % num
    reflections_file = "reflections_%02d.pickle" % num
    experiments_path = os.path.join(datadir, experiments_file)
    reflections_path = os.path.join(datadir, reflections_file)

    # create directory
    with cd("sweep_%02d" % num):
      # WARNING! Fast and dirty integration.
      # Do not use the result for scaling/merging!
      cmd = "dials.integrate %s %s " + \
            "intensity.algorithm=sum prediction.dmin=3 prediction.dmax=8"
      cmd = cmd % (experiments_path, reflections_path)
      easy_run.fully_buffered(command=cmd)
      if not os.path.isfile("integrated.pickle"):
        print "Job %02d failed during integration" % num
        return

      # create MTZ
      cmd = "dials.export_mtz %s integrated.pickle hklout=integrated.mtz"
      cmd = cmd % experiments_path
      easy_run.fully_buffered(command=cmd)
      if not os.path.isfile("integrated.mtz"):
        print "Job %02d failed during MTZ export" % num
        return

    # if we got this far, return the path to the MTZ
    return "sweep_%02d/integrated.mtz" % num

  if __name__ == "__main__":

    if len(sys.argv) != 2:
      sys.exit("Usage: dials.python integrate_joint_TehA.py ..")
    data_dir = os.path.abspath(sys.argv[1])

    pathname = os.path.join(data_dir, "experiments_*.json")
    experiments = glob.glob(pathname)

    templates = [data_dir for f in experiments]
    tasklist = list(enumerate(sorted(templates)))

    from libtbx import Auto
    nproc = easy_mp.get_processes(Auto)

    print "Attempting to process the following datasets, with {} processes".format(nproc)
    for task in tasklist:
      print "%d: %s/experiments%02d" % (task[0], task[1], task[0])

    results = easy_mp.parallel_map(
      func=process_sweep,
      iterable=tasklist,
      processes=nproc,
      preserve_order=True)

    good_results = [e for e in results if e is not None]
    print "Successfully created the following MTZs:"
    for result in good_results:
      print result

This, if saved as :file:`integrate_joint_TehA.py` in the new :file:`joint`
directory can be run as follows::

  dials.python integrate_joint_TehA.py .

As expected this creates all 39 MTZs for the jointly refined sweeps without any
problem. We can copy the paths to these into a new file, say
:file:`joint_mtzs.dat`, and run blend::

  echo "END" | blend -a joint_mtzs.dat

The :file:`tree.png` resulting from this is very interesting.

  .. image:: figures/tree_03.png

The LCV is now as low as 0.36% (aLCV 0.27 Angstroms). This indicates an even
higher degree of isomorphism than detected during after individual processing.
So although joint refinement leads to slightly higher RMSDs for each experiment
(as we expected) the resulting unit cells are more similar. It is worth
remembering that no restraints were applied between unit cells in refinement.
Given that we know that the detector and beam did not move between the data
collections we might like to think that the joint refinement analysis is a more
accurate depiction of reality, and thus the unit cells are closer to the truth.

What to do next?
----------------

FIXME e.g.

* go back and fix datasets that didn't index properly
* integrate data properly for blend synthesis mode

Acknowledgements
----------------

Danny Axford, Nien-Jen Hu, James Foadi, Hassanul Ghani Choudhury, So Iwata, Konstantinos Beis, Gwyndaf Evans & Yilmaz Alguel

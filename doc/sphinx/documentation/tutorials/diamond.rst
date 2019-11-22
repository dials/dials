.. raw:: html

  <a href="https://dials.github.io/dials-1.14/documentation/tutorials/diamond.html" class="new-documentation">
  This tutorial requires a DIALS 2.0 installation.<br/>
  Please click here to go to the tutorial for DIALS 1.14.
  </a>

###################################
Using DIALS at Diamond Light Source
###################################

.. _introduction:

DIALS is centrally installed on the Linux computer systems at Diamond. You can
load the default version of DIALS by typing at the command line::

   $ module load dials

DIALS 1.14 is available by typing::

   $ module load dials/1

The default version of DIALS (and xia2, etc.) is usually fixed at the start of
a run, however a more recent nightly build of DIALS may be available::

   $ module load dials/latest

We suggest using this latest version of DIALS when trying the tutorials.

General
=======

* `Main documentation page`_

.. _Main documentation page: https://dials.github.io/

* `Tutorials`_

.. _Tutorials: https://dials.github.io/documentation/tutorials/index.html

* `Program documentation`_

.. _Program documentation: https://dials.github.io/documentation/programs/index.html

Using DIALS with xia2
=====================

As an alternative to running DIALS directly, DIALS can be run within xia2_,
using :samp:`pipeline=dials`::

   $ xia2 pipeline=dials /dls/mx-scratch/dials_tutorial/i04_thaumatin/
   Environment configuration...
   Python => /dls/science/groups/scisoft/DIALS/CD/latest/dials-dev20191028/build/../base/bin/python2.7
   CCTBX => /dls/science/groups/scisoft/DIALS/CD/latest/dials-dev20191028/modules
   CCP4 => /dls_sw/apps/ccp4/7.0.075/ccp4-7.0
   CLIBD => /dls_sw/apps/ccp4/7.0.075/ccp4-7.0/lib/data
   CCP4_SCR => /tmp/tmp4Xxusn
   Starting directory: /dls/tmp/thaumatin
   Working directory: /dls/tmp/thaumatin
   Free space:        4229.16 GB
   Host: cs03r-sc-serv-16
   Contact: xia2.support@gmail.com
   XIA2 0.6.257-g73f70f18-remotes/origin/dials-2.0
   DIALS 2.dev.1038-gf37b0db7a
   CCP4 7.0.075
   Command line: xia2 pipeline=dials /dls/mx-scratch/dials_tutorial/i04_thaumatin
   -------------------- Spotfinding SWEEP1 --------------------
   120152 spots found on 540 images (max 2197 / bin)
   *                                       *  * **    *  **  **
   ************************************************************
   ************************************************************
   ************************************************************
   ************************************************************
   ************************************************************
   ************************************************************
   ************************************************************
   ************************************************************
   ************************************************************
   1                         image                          540
   ------------------- Autoindexing SWEEP1 --------------------
   All possible indexing solutions:
   tP  57.78  57.78 150.00  90.00  90.00  90.00
   oC  81.72  81.74 150.02  90.00  90.00  90.00
   oP  57.76  57.79 149.99  90.00  90.00  90.00
   mC  81.74  81.75 150.04  90.00  89.99  90.00
   mP  57.77  57.80 150.01  90.00  89.98  90.00
   aP  57.78  57.81 150.03  89.99  89.99  90.01
   Indexing solution:
   tP  57.78  57.78 150.00  90.00  90.00  90.00
   -------------------- Integrating SWEEP1 --------------------
   Processed batches 1 to 540
   Standard Deviation in pixel range: 0.41 0.51
   Integration status per image (60/record):
   %ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
   "o" => good        "%" => ok        "!" => bad rmsd
   "O" => overloaded  "#" => many bad  "." => weak
   "@" => abandoned
   Mosaic spread: 0.043 < 0.043 < 0.043
   -------------------- Preparing DEFAULT ---------------------
   --------------------- Scaling DEFAULT ----------------------
   Completed a round of scaling using dials.scale
   Resolution for sweep NATIVE/SWEEP1: 1.17 (cc_half > 0.3, unmerged <I/sigI> > 0.25)
   --------------------- Scaling DEFAULT ----------------------
   Completed a round of scaling using dials.scale
   ---------------- Systematic absences check -----------------
   Most likely space group: P 41 21 2
   ------------------- Unit cell refinement -------------------
   Overall:  57.78  57.78 149.99  90.00  90.00  90.00
   Overall twinning score: 2.09
   Your data do not appear twinned
   Project: AUTOMATIC
   Crystal: DEFAULT
   Sequence:
   Wavelength name: NATIVE
   Wavelength 0.97625
   Sweeps:
   SWEEP SWEEP1 [WAVELENGTH NATIVE]
   TEMPLATE th_8_2_####.cbf
   DIRECTORY /dls/mx-scratch/dials_tutorial/i04_thaumatin
   IMAGES (USER) 1 to 540
   MTZ file: /dls/tmp/thaumatin/DEFAULT/NATIVE/SWEEP1/integrate/13_integrated.refl
   For AUTOMATIC/DEFAULT/NATIVE                 Overall    Low     High
   High resolution limit                           1.17    3.18    1.17
   Low resolution limit                           53.92   53.98    1.19
   Completeness                                   82.2   100.0     7.8
   Multiplicity                                    4.5     5.4     1.0
   I/sigma                                        10.0    46.6     0.4
   Rmerge(I)                                     0.070   0.032   0.356
   Rmerge(I+/-)                                  0.062   0.029   0.000
   Rmeas(I)                                      0.078   0.035   0.503
   Rmeas(I+/-)                                   0.077   0.035   0.000
   Rpim(I)                                       0.034   0.015   0.356
   Rpim(I+/-)                                    0.043   0.020   0.000
   CC half                                       0.998   0.998   0.700
   Wilson B factor                               9.005
   Anomalous completeness                         71.3    99.2     0.3
   Anomalous multiplicity                          2.6     3.1     1.0
   Anomalous correlation                        -0.012  -0.060   0.000
   Anomalous slope                               0.705
   Total observations                           321472   25415     339
   Total unique                                  70801    4729     328
   Assuming spacegroup: P 41 21 2
   Unit cell (with estimated std devs):
   57.78346(8) 57.78346(8) 149.9943(3)
   90.0        90.0         90.0
   mtz format:
   Scaled reflections: /dls/tmp/thaumatin/DataFiles/AUTOMATIC_DEFAULT_free.mtz
   mtz_unmerged format:
   Scaled reflections (NATIVE): /dls/tmp/thaumatin/DataFiles/AUTOMATIC_DEFAULT_scaled_unmerged.mtz
   Processing took 00h 06m 20s
   XIA2 used... ccp4, dials, xia2
   Here are the appropriate citations (BIBTeX in xia2-citations.bib.)
   Winn, M. D. et al. (2011) Acta Cryst. D67, 235-242.
   Winter, G. (2010) J. Appl. Cryst. 43, 186-190.
   Winter, G. et al. (2018) Acta Cryst. D74, 85-97.
   Status: normal termination


Feedback
========

If you encounter any problems using DIALS or xia2, or indeed have any other
feedback (positive or negative - we love to hear from you either way),
please contact dials-support@lists.sourceforge.net or xia2.support@gmail.com
respectively.

.. _xia2: https://xia2.github.io/

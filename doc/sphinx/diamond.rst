
+++++++++++++++++++++++++++++++++++
Using DIALS at Diamond Light Source
+++++++++++++++++++++++++++++++++++

.. _introduction:

DIALS is centrally installed on the Linux computer systems at Diamond. You can
load the default version of DIALS by typing at the command line::

  $ module load dials

The default version of DIALS (and xia2, etc.) is usually fixed at the start of
a run, however more recent nightly builds of DIALS may be available::

  $ module avail dials

  ------------------------------------------------- /dls_sw/apps/Modules/modulefiles --------------------------------------------------
  dials/2014-01-16       dials/2014-07-03       dials/dev-105          dials/dev-208          dials/svn
  dials/2014-01-17       dials/2014-07-21       dials/dev-124          dials/dev-83
  dials/2014-04-10       dials/2014_09_19_0127  dials/dev-127          dials/dev-95
  dials/2014-04-28       dials/dev-100          dials/dev-194          dials/dls-434(default)

Tutorials
=========

For the following tutorials we suggest using the latest version of DIALS::

  $ module load dials/dev-208

The thaumatin dataset used in the :doc:`advanced_tutorial` is centrally
available on the Diamond filesystem at
:samp:`/dls/mx-scratch/dials_tutorial/i04_thaumatin/`.

.. toctree::
   :maxdepth: 1

   advanced_tutorial
   multi_crystal_analysis
   metrology_corrections

General
=======

.. toctree::
   :maxdepth: 1

   Main documentation page <index>
   programs/index
   conventions
   library_reference/index

Using DIALS with xia2
=====================

.. note::

   Due to the rapidly evolving nature of the DIALS project, more recent
   versions of DIALS may not be compatible with the default version of xia2.
   We plan to bundle xia2 with DIALS in the near future, so this problem will
   go away.

As an alternative to running DIALS directly, support for DIALS has been added
within xia2_. You can use xia2_ just as before, simply replacing
:samp:`-2d`, :samp:`-3d`, etc. with :samp:`-dials`::

  $ module load xia2
  $ module load dials
  $ xia2 -dials /dls/mx-scratch/dials_tutorial/i04_thaumatin/
  Environment configuration...
  XIA2_ROOT => /dls_sw/apps/xia2/5006/xia2
  XIA2CORE_ROOT => /dls_sw/apps/xia2/5006/xia2/core
  Python => /dls_sw/apps/dials/dls-434/dials-dev-2015_01_05_0757/build/../base/bin/python2.7
  CCTBX => /dls_sw/apps/dials/dls-434/dials-dev-2015_01_05_0757/modules
  CCP4 => /dls_sw/apps/ccp4/x86_64/6.4.0/11nov2014b/ccp4-6.4.0
  CLIBD => /dls_sw/apps/ccp4/x86_64/6.4.0/11nov2014b/ccp4-6.4.0/lib/data
  CCP4_SCR => /tmp/tmpfCl8HK
  Working directory: /dls/mx-scratch/rjgildea/tmp/insulin
  Free space:        113794.78 GB
  Host: ws133
  Build: 4992
  Contact: xia2.support@gmail.com
  XIA2 0.3.8.0
  Command line: xia2 -dials /dls/mx-scratch/dials_tutorial/i04_thaumatin/
  ------------------- Autoindexing SWEEP1 --------------------
  All possible indexing solutions:
  tP  57.79  57.79 150.01  90.00  90.00  90.00
  oC  81.72  81.73 150.01  90.00  90.00  90.00
  oP  57.80  57.78 150.01  90.00  90.00  90.00
  mC  81.72  81.73 150.01  90.00  89.99  90.00
  mP  57.80  57.78 150.01  90.00  89.99  90.00
  aP  57.80  57.78 150.01  89.99  89.99  90.00
  Indexing solution:
  tP  57.79  57.79 150.01  90.00  90.00  90.00
  -------------------- Integrating SWEEP1 --------------------
  Integration status per image (60/record):
  oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
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
  Mosaic spread: 0.069 < 0.069 < 0.069
  -------------------- Preparing DEFAULT ---------------------
  Likely spacegroups:
  P 41 21 2
  P 43 21 2
  Reindexing to first spacegroup setting: P 41 21 2 (h,k,l)
  --------------------- Scaling DEFAULT ----------------------
  Resolution limit for NATIVE/SWEEP1:  1.35
  --------------------- Scaling DEFAULT ----------------------
  Resolution limit for NATIVE/SWEEP1:  1.35
  Overall twinning score: 2.30
  Your data do not appear to be twinned
  Project: AUTOMATIC
  Crystal: DEFAULT
  Sequence:
  Wavelength name: NATIVE
  Wavelength 0.97625
  Sweeps:
  SWEEP SWEEP1 [WAVELENGTH NATIVE]
  TEMPLATE th_8_2_####.cbf
  DIRECTORY /dls/mx-scratch/dials_tutorial/data
  IMAGES (USER) 1 to 540
  MTZ file: /dls/mx-scratch/rjgildea/tmp/insulin/DEFAULT/NATIVE/SWEEP1/integrate/10_integrated.mtz
  For AUTOMATIC/DEFAULT/NATIVE
  High resolution limit                   	  1.35	  5.23	  1.35
  Low resolution limit                    	 50.00	 50.00	  1.40
  Completeness                            	 99.7	 99.9	 97.5
  Multiplicity                            	  5.2	  4.9	  3.0
  I/sigma                                 	 12.8	 32.9	  2.1
  Rmerge                                  	0.063	0.026	0.371
  Rmeas(I)                                	0.078	0.033	0.524
  Rmeas(I+/-)                             	0.077	0.032	0.497
  Rpim(I)                                 	0.033	0.014	0.283
  Rpim(I+/-)                              	0.044	0.018	0.328
  CC half                                 	0.999	0.999	0.810
  Wilson B factor                         	8.469
  Anomalous completeness                  	 96.9	 99.9	 75.7
  Anomalous multiplicity                  	  2.6	  3.0	  1.8
  Anomalous correlation                   	-0.007	 0.089	 0.014
  Anomalous slope                         	0.956	0.000	0.000
  Total observations                      	293318	5742	16079
  Total unique                            	56672	1160	5308
  Assuming spacegroup: P 41 21 2
  Other likely alternatives are:
  P 43 21 2
  Unit cell:
  57.788  57.788 150.010
  90.000  90.000  90.000
  mtz format:
  Scaled reflections: /dls/mx-scratch/rjgildea/tmp/insulin/DataFiles/AUTOMATIC_DEFAULT_free.mtz
  sca format:
  Scaled reflections (NATIVE): /dls/mx-scratch/rjgildea/tmp/insulin/DataFiles/AUTOMATIC_DEFAULT_scaled.sca
  sca_unmerged format:
  Scaled reflections (NATIVE): /dls/mx-scratch/rjgildea/tmp/insulin/DataFiles/AUTOMATIC_DEFAULT_scaled_unmerged.sca
  Processing took 00h 07m 34s
  XIA2 used...  aimless ccp4 dials distl labelit pointless xia2
  Here are the appropriate citations (BIBTeX in xia-citations.bib.)
  (1994) Acta Crystallogr. D 50, 760--763
  Evans, Philip (2006) Acta Crystallographica Section D 62, 72--82
  Sauter, Nicholas K. and Grosse-Kunstleve, Ralf W. and Adams, Paul D. (2004) Journal of Applied Crystallography 37, 399--409
  Winter, G. (2010) Journal of Applied Crystallography 43
  Zhang, Z. and Sauter, N.K. and van den Bedem, H. and Snell, G. and Deacon, A.M. (2006) J. Appl. Cryst 39, 112--119
  Status: normal termination


Feedback
========

If you encounter any problems using DIALS or xia2, or indeed have any other
feedback (positive or negative - we love to hear from you either way),
please contact dials-support@lists.sourceforge.net or xia2.support@gmail.com
respectively.

.. _xia2: http://www.ccp4.ac.uk/xia/

+++++++++++++++++++++++++++++++++++
Using DIALS at Diamond Light Source
+++++++++++++++++++++++++++++++++++

.. _introduction:

DIALS is centrally installed on the Linux computer systems at Diamond. You can
load the default version of DIALS by typing at the command line::

   $ module load dials

The default version of DIALS (and xia2, etc.) is usually fixed at the start of
a run, however a more recent nightly build of DIALS may be available::

   $ module load dials/latest

Tutorials
=========

For the following tutorials we suggest using the latest version of DIALS::

   $ module load dials/latest

The thaumatin dataset used in the
:doc:`processing_in_detail_tutorial` is centrally
available on the Diamond filesystem at
:samp:`/dls/mx-scratch/dials_tutorial/i04_thaumatin/`.

.. toctree::
   :maxdepth: 1

   processing_in_detail_tutorial
   multi_crystal_analysis
   metrology_corrections

General
=======

* `Main documentation page`_

.. _Main documentation page: https://dials.github.io/

* `Program documentation`_

.. _Program documentation: https://dials.github.io/documentation/programs/index.html

Using DIALS with xia2
=====================

As an alternative to running DIALS directly, support for DIALS has been added
within xia2_. You can use xia2_ just as before, simply replacing
:samp:`-2d`, :samp:`-3d`, etc. with :samp:`-dials`::

   $ xia2 -dials /dls/mx-scratch/dials_tutorial/i04_thaumatin/
   Environment configuration...
   XIA2_ROOT => /dls_sw/apps/dials/latest/dials-dev20151124/modules/xia2
   XIA2CORE_ROOT => /dls_sw/apps/dials/latest/dials-dev20151124/modules/xia2/core
   Python => /dls_sw/apps/dials/latest/dials-dev20151124/build/../base/bin/python2.7
   CCTBX => /dls_sw/apps/dials/latest/dials-dev20151124/modules
   CCP4 => /dls_sw/apps/ccp4/64/6.5/update19/ccp4-6.5
   CLIBD => /dls_sw/apps/ccp4/64/6.5/update19/ccp4-6.5/lib/data
   CCP4_SCR => /tmp/tmpPJJ4I7
   Working directory: /dls/mx-scratch/rjgildea/tmp/insulin
   Free space:        175586.36 GB
   Host: cs03r-sc-serv-16
   Build: v0.4.0.0-94-gf62dc5d
   Contact: xia2.support@gmail.com
   XIA2 0.4.0.0
   DIALS 1.0-pre-210-gdd8f759
   Command line: xia2 -dials /dls/mx-scratch/dials_tutorial/i04_thaumatin/
   ------------------- Autoindexing SWEEP1 --------------------
   All possible indexing solutions:
   tP  57.79  57.79 150.01  90.00  90.00  90.00
   oC  81.72  81.74 150.02  90.00  90.00  90.00
   oP  57.80  57.77 150.01  90.00  90.00  90.00
   mC  81.73  81.75 150.03  90.00  89.99  90.00
   mP  57.77  57.80 150.02  90.00  90.02  90.00
   aP  57.81  57.78 150.03  90.01  89.99  89.99
   Indexing solution:
   tP  57.79  57.79 150.01  90.00  90.00  90.00
   -------------------- Integrating SWEEP1 --------------------
   Processed batches 1 to 540
   Standard Deviation in pixel range: 0.370000 0.450000
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
   Mosaic spread: 0.066 < 0.066 < 0.066
   -------------------- Preparing DEFAULT ---------------------
   Likely spacegroups:
   P 41 21 2
   P 43 21 2
   Reindexing to first spacegroup setting: P 41 21 2 (h,k,l)
   --------------------- Scaling DEFAULT ----------------------
   Resolution limit for NATIVE/SWEEP1:  1.36
   --------------------- Scaling DEFAULT ----------------------
   Resolution limit for NATIVE/SWEEP1:  1.36
   Overall twinning score: 2.02
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
   MTZ file: /dls/mx-scratch/rjgildea/tmp/insulin/DEFAULT/NATIVE/SWEEP1/integrate/dials_integrated.mtz
   For AUTOMATIC/DEFAULT/NATIVE
   High resolution limit                   	  1.36	  6.08	  1.36
   Low resolution limit                    	 50.00	 50.00	  1.40
   Completeness                            	 99.8	 99.8	 97.8
   Multiplicity                            	  5.3	  4.9	  3.1
   I/sigma                                 	 12.1	 31.7	  2.1
   Rmerge                                  	0.060	0.025	0.364
   Rmeas(I)                                	0.074	0.031	0.508
   Rmeas(I+/-)                             	0.074	0.031	0.488
   Rpim(I)                                 	0.031	0.014	0.272
   Rpim(I+/-)                              	0.042	0.017	0.321
   CC half                                 	0.999	0.999	0.826
   Wilson B factor                         	8.873
   Anomalous completeness                  	 97.5	100.0	 77.3
   Anomalous multiplicity                  	  2.6	  3.1	  1.8
   Anomalous correlation                   	-0.001	 0.299	-0.024
   Anomalous slope                         	0.960	0.000	0.000
   Total observations                      	292195	3740	12288
   Total unique                            	55487	767	3923
   Assuming spacegroup: P 41 21 2
   Other likely alternatives are:
   P 43 21 2
   Unit cell:
   57.784  57.784 150.002
   90.000  90.000  90.000
   mtz format:
   Scaled reflections: /dls/mx-scratch/rjgildea/tmp/insulin/DataFiles/AUTOMATIC_DEFAULT_free.mtz
   mtz_unmerged format:
   Scaled reflections (NATIVE): /dls/mx-scratch/rjgildea/tmp/insulin/DataFiles/AUTOMATIC_DEFAULT_scaled_unmerged.mtz
   sca format:
   Scaled reflections (NATIVE): /dls/mx-scratch/rjgildea/tmp/insulin/DataFiles/AUTOMATIC_DEFAULT_scaled.sca
   sca_unmerged format:
   Scaled reflections (NATIVE): /dls/mx-scratch/rjgildea/tmp/insulin/DataFiles/AUTOMATIC_DEFAULT_scaled_unmerged.sca
   Processing took 00h 13m 16s
   XIA2 used...  aimless ccp4 dials pointless xia2
   Here are the appropriate citations (BIBTeX in xia-citations.bib.)
   (1994) Acta Crystallogr. D 50, 760--763
   Evans, Philip (2006) Acta Crystallographica Section D 62, 72--82
   Winter, G. (2010) Journal of Applied Crystallography 43
   Status: normal termination


Feedback
========

If you encounter any problems using DIALS or xia2, or indeed have any other
feedback (positive or negative - we love to hear from you either way),
please contact dials-support@lists.sourceforge.net or xia2.support@gmail.com
respectively.

.. _xia2: https://xia2.github.io/

.. raw:: html

  <a href="https://dials.github.io/dials-2.2/documentation/tutorials/multi_lattice_tutorial.html" class="new-documentation">
  This tutorial requires a DIALS 3 installation.<br/>
  Please click here to go to the tutorial for DIALS 2.2.
  </a>

Multi-lattice Tutorial
======================

.. highlight:: none

Introduction
------------

The following example uses semi-synthetic multi-lattice trypsin datasets
collected using beamline I04 at Diamond Light Source which is available for
download from |semisynthetic|. This tutorial uses one of the 2-lattice
datasets which are contained in the file `semisynthetic_multilattice_data_2.tar.bz2`_.

.. _semisynthetic_multilattice_data_2.tar.bz2: https://zenodo.org/record/10820/files/semisynthetic_multilattice_data_2.tar.bz2

.. |semisynthetic| image:: https://zenodo.org/badge/doi/10.5281/zenodo.10820.svg
               :target: https://doi.org/10.5281/zenodo.10820

In this tutorial we shall focus on the processing steps that diverge from the
regular single-lattice processing discussed in :doc:`processing_in_detail_betalactamase`.

Import and Spotfinding
^^^^^^^^^^^^^^^^^^^^^^

As for single-lattice processing, the first steps are to import the data and
find spots using the following commands::

  dials.import semisynthetic_multilattice_data/2/ag/trp_ag_*.cbf

  dials.find_spots imported.expt min_spot_size=3

During import, all that happens here is that the image headers are read, and a
file describing their contents (:samp:`imported.expt`) is written. The output
just describes what the software understands of the images it was
passed, in this case one sequence of data containing 100 images::

  The following parameters have been modified:

  input {
    experiments = <image files>
  }

  --------------------------------------------------------------------------------
    format: <class 'dxtbx.format.FormatCBFMiniPilatusDLS6MSN100.FormatCBFMiniPilatusDLS6MSN100'>
    num images: 100
    sequences:
      still:    0
      sweep:    1
    num stills: 0
  --------------------------------------------------------------------------------
  Writing experiments to imported.expt

For the spot finding, we tweak the minimum spot size (min_spot_size=3) to improve
the results for this dataset::

  Extracted 46332 spots
  Removed 8863 spots with size < 3 pixels
  Removed 3 spots with size > 1000 pixels
  Calculated 37466 spot centroids
  Calculated 37466 spot intensities
  Filtered 35422 of 37466 spots by peak-centroid distance

  Histogram of per-image spot count for imageset 0:
  35422 spots found on 100 images (max 1190 / bin)
                                                             *
  *                                                          *
  *                                                          *
  *                                                          *
  **                                                         *
  **** ******************* * ************* ******** ** *******
  ************************************************************
  ************************************************************
  ************************************************************
  ************************************************************
  1                         image                          100

  --------------------------------------------------------------------------------
  Saved 35422 reflections to strong.refl


Indexing
^^^^^^^^
The next step is the indexing of the strong spots. By default only one
lattice is searched for, but if there are sufficient unindexed reflections
remaining after indexing the first lattice, we can switch on indexing of
multiple lattices using the parameter max_lattices=2 (e.g.)::

  dials.index imported.expt strong.refl max_lattices=2

::

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  | id  |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 1000 | 0.41832 | 0.25232 | 0.15582  |
  | 1   | 1000 | 0.33596 | 0.24331 | 0.17531  |
  ---------------------------------------------

  Refined crystal models:
  model 1 (16636 reflections):
  Crystal:
      Unit cell: (54.063(4), 58.2475(18), 66.494(2), 89.9778(13), 90.013(3), 90.012(3))
      Space group: P 1
      U matrix:  {{ 0.1870,  0.7632, -0.6185},
                  { 0.0427,  0.6227,  0.7813},
                  { 0.9814, -0.1725,  0.0838}}
      B matrix:  {{ 0.0185,  0.0000,  0.0000},
                  { 0.0000,  0.0172,  0.0000},
                  { 0.0000, -0.0000,  0.0150}}
      A = UB:    {{ 0.0035,  0.0131, -0.0093},
                  { 0.0008,  0.0107,  0.0117},
                  { 0.0182, -0.0030,  0.0013}}
  model 2 (17247 reflections):
  Crystal:
      Unit cell: (54.080(3), 58.263(2), 66.498(3), 90.0060(18), 90.021(3), 90.045(3))
      Space group: P 1
      U matrix:  {{ 0.0094,  0.6714, -0.7410},
                  { 0.3813, -0.6875, -0.6180},
                  {-0.9244, -0.2768, -0.2625}}
      B matrix:  {{ 0.0185,  0.0000,  0.0000},
                  { 0.0000,  0.0172,  0.0000},
                  { 0.0000,  0.0000,  0.0150}}
      A = UB:    {{ 0.0002,  0.0115, -0.0111},
                  { 0.0070, -0.0118, -0.0093},
                  {-0.0171, -0.0048, -0.0039}}
  --------------------------------------------------
  | Imageset | # indexed | # unindexed | % indexed |
  --------------------------------------------------
  | 0        | 33883     | 1539        | 95.7%     |
  --------------------------------------------------
  Change of basis op: a,b,c
  Rotation matrix to transform crystal 1 to crystal 2:
  {{0.973, -0.160, -0.169},
  {-0.071, -0.895, 0.441},
  {-0.222, -0.417, -0.881}}
  Rotation of -154.399 degrees about axis (0.993, -0.061, -0.103)

  Saving refined experiments to indexed.expt
  Saving refined reflections to indexed.refl

Next we run
:doc:`dials.refine_bravais_settings </documentation/programs/dials_refine_bravais_settings>`
refining each indexing solution (separately) in all Bravais settings
consistent with the indexed unit cell. In this example we would continue
processing using bravais_setting_5.expt, i.e. solution number 5.

::

  dials.refine_bravais_settings indexed.expt indexed.refl crystal_id=0

  dials.refine_bravais_settings indexed.expt indexed.refl crystal_id=1

gives a table containing the metric fit, rmsds (in mm) and unit cell for
each Bravais setting...

::

  ----------------------------------------------------------------------------------------------------------------
  Solution Metric fit  rmsd  min/max cc #spots lattice                                 unit_cell volume      cb_op
  ----------------------------------------------------------------------------------------------------------------
        9     4.2490 1.579 0.384/0.763   1000      tP  60.31  60.31  69.15  90.00  90.00  90.00 251517      a,b,c
        8     4.2490 1.508 0.372/0.529   1000      oC  85.65  84.55  69.06  90.00  90.00  90.00 500054 a+b,-a+b,c
        7     4.2490 1.487 0.372/0.372   1000      mC  84.57  85.57  69.01  90.00  89.92  90.00 499425  a-b,a+b,c
        6     4.2490 1.560 0.529/0.529   1000      mC  85.71  84.31  69.03  90.00  89.87  90.00 498876 a+b,-a+b,c
  *     5     0.0000 0.095 0.245/0.904   1000      oP  54.10  58.27  66.51  90.00  90.00  90.00 209657      a,b,c
  *     4     0.0000 0.088 0.904/0.904   1000      mP  58.27  54.11  66.52  90.00  89.98  90.00 209735   -b,-a,-c
  *     3     0.0000 0.095 0.245/0.245   1000      mP  54.11  58.27  66.51  90.00  90.02  90.00 209715      a,b,c
  *     2     0.0000 0.093 0.384/0.384   1000      mP  54.11  66.52  58.28  90.00  90.02  90.00 209787   -a,-c,-b
  *     1     0.0000 0.090         -/-   1000      aP  54.11  58.27  66.51  89.98  90.02  90.00 209725      a,b,c
  ----------------------------------------------------------------------------------------------------------------
  ----------------------------------------------------------------------------------------------------------------
  Solution Metric fit  rmsd  min/max cc #spots lattice                                 unit_cell volume      cb_op
  ----------------------------------------------------------------------------------------------------------------
        9     4.2639 1.740 0.227/0.833   1000      tP  59.33  59.33  68.32  90.00  90.00  90.00 240463      a,b,c
        8     4.2639 1.709 0.242/0.833   1000      oC  84.84  83.85  68.54  90.00  90.00  90.00 487536 a+b,-a+b,c
        7     4.2639 1.343 0.631/0.631   1000      mC  83.78  82.53  68.28  90.00  88.90  90.00 472000  a-b,a+b,c
        6     4.2639 1.569 0.242/0.242   1000      mC  82.84  83.47  67.49  90.00  91.40  90.00 466492 a+b,-a+b,c
  *     5     0.0497 0.076 0.658/0.833   1000      oP  54.10  58.29  66.52  90.00  90.00  90.00 209775      a,b,c
  *     4     0.0497 0.078 0.658/0.658   1000      mP  58.29  54.10  66.52  90.00  89.99  90.00 209773   -b,-a,-c
  *     3     0.0453 0.075 0.811/0.811   1000      mP  54.09  58.29  66.51  90.00  90.02  90.00 209673      a,b,c
  *     2     0.0221 0.075 0.833/0.833   1000      mP  54.09  66.51  58.27  90.00  90.03  90.00 209642   -a,-c,-b
  *     1     0.0000 0.075         -/-   1000      aP  54.08  58.27  66.50  90.01  90.02  90.03 209577      a,b,c
  ----------------------------------------------------------------------------------------------------------------

Now we re-run the indexing, this time imposing the lattice constraints for
the chosen Bravais setting, in this case number 5, i.e. oP, or point group
P222.

::

  dials.index imported.expt strong.refl max_lattices=2 space_group=P222

::

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  | id  |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 1000 | 0.45694 | 0.26413 | 0.16573  |
  | 1   | 1000 | 0.35593 | 0.27804 | 0.20011  |
  ---------------------------------------------

  Refined crystal models:
  model 1 (16635 reflections):
  Crystal:
      Unit cell: (54.100(4), 58.2684(17), 66.517(2), 90.0, 90.0, 90.0)
      Space group: P 2 2 2
      U matrix:  {{ 0.1873,  0.7630, -0.6186},
                  { 0.0429,  0.6228,  0.7812},
                  { 0.9814, -0.1728,  0.0839}}
      B matrix:  {{ 0.0185,  0.0000,  0.0000},
                  {-0.0000,  0.0172,  0.0000},
                  {-0.0000,  0.0000,  0.0150}}
      A = UB:    {{ 0.0035,  0.0131, -0.0093},
                  { 0.0008,  0.0107,  0.0117},
                  { 0.0181, -0.0030,  0.0013}}
  model 2 (17249 reflections):
  Crystal:
      Unit cell: (54.117(3), 58.2882(15), 66.526(2), 90.0, 90.0, 90.0)
      Space group: P 2 2 2
      U matrix:  {{ 0.0091,  0.6714, -0.7410},
                  { 0.3810, -0.6874, -0.6182},
                  {-0.9245, -0.2768, -0.2621}}
      B matrix:  {{ 0.0185,  0.0000,  0.0000},
                  {-0.0000,  0.0172,  0.0000},
                  { 0.0000,  0.0000,  0.0150}}
      A = UB:    {{ 0.0002,  0.0115, -0.0111},
                  { 0.0070, -0.0118, -0.0093},
                  {-0.0171, -0.0047, -0.0039}}
  --------------------------------------------------
  | Imageset | # indexed | # unindexed | % indexed |
  --------------------------------------------------
  | 0        | 33884     | 1538        | 95.7%     |
  --------------------------------------------------
  Change of basis op: -a,b,-c
  Rotation matrix to transform crystal 1 to crystal 2:
  {{0.052, 0.997, -0.063},
  {-0.978, 0.038, -0.203},
  {-0.200, 0.072, 0.977}}
  Rotation of -88.056 degrees about axis (-0.138, -0.069, 0.988)

  Saving refined experiments to indexed.expt
  Saving refined reflections to indexed.refl


Refinement and Integration
^^^^^^^^^^^^^^^^^^^^^^^^^^

After indexing, processing proceeds similarly to the single-lattice case.
First, the crystal models can be further refined with a scan varying model,
in this example also using the tukey outlier rejection algorithm::

  dials.refine indexed.expt indexed.refl outlier.algorithm=tukey

::

  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 2000 | 0.079758 | 0.046104 | 0.018187 |
  | 1    | 2000 | 0.066176 | 0.042452 | 0.017411 |
  | 2    | 2000 | 0.065727 | 0.042236 | 0.016897 |
  | 3    | 2000 | 0.065412 | 0.042413 | 0.016653 |
  | 4    | 2000 | 0.065282 | 0.042591 | 0.016592 |
  | 5    | 2000 | 0.065257 | 0.042631 | 0.016585 |
  | 6    | 2000 | 0.065253 | 0.042632 | 0.016585 |
  ------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  ----------------------------------------------
  | Exp | Nref  | RMSD_X  | RMSD_Y  | RMSD_Z   |
  | id  |       | (px)    | (px)    | (images) |
  ----------------------------------------------
  | 0   | 13833 | 0.48806 | 0.27504 | 0.1626   |
  | 1   | 14752 | 0.35492 | 0.2951  | 0.20587  |
  ----------------------------------------------
  Updating predictions for indexed reflections
  Saving refined experiments to refined.expt
  Saving reflections with updated predictions to refined.refl

Next, we integrate the data::

  dials.integrate refined.expt refined.refl

This program outputs a lot of information as integration progresses,
concluding with a summary of the integration results for each lattice::

  Summary for experiment 0
  ---------------------------------------------------------------
  Item                                  | Overall | Low    | High
  ---------------------------------------------------------------
  dmin                                  | 1.06    | 2.87   | 1.06
  dmax                                  | 43.87   | 43.87  | 1.08
  number fully recorded                 | 25027   | 1859   | 29
  number partially recorded             | 9850    | 760    | 6
  number with invalid background pixels | 10448   | 623    | 29
  number with invalid foreground pixels | 5595    | 394    | 17
  number with overloaded pixels         | 4       | 4      | 0
  number in powder rings                | 0       | 0      | 0
  number processed with summation       | 29114   | 2208   | 18
  number processed with profile fitting | 23507   | 1796   | 6
  number failed in background modelling | 20      | 5      | 0
  number failed in summation            | 5595    | 394    | 17
  number failed in profile fitting      | 11202   | 806    | 29
  ibg                                   | 17.81   | 49.24  | 4.31
  i/sigi (summation)                    | 26.91   | 151.98 | 1.32
  i/sigi (profile fitting)              | 35.42   | 199.63 | 1.28
  cc prf                                | 0.94    | 0.85   | 0.98
  cc_pearson sum/prf                    | 0.89    | 0.86   | 0.93
  cc_spearman sum/prf                   | 0.97    | 0.98   | 0.37
  ---------------------------------------------------------------

  Summary for experiment 1
  ---------------------------------------------------------------
  Item                                  | Overall | Low    | High
  ---------------------------------------------------------------
  dmin                                  | 1.06    | 2.87   | 1.06
  dmax                                  | 25.49   | 25.49  | 1.08
  number fully recorded                 | 24816   | 1800   | 38
  number partially recorded             | 10172   | 818    | 9
  number with invalid background pixels | 9013    | 539    | 34
  number with invalid foreground pixels | 5147    | 369    | 23
  number with overloaded pixels         | 3       | 3      | 0
  number in powder rings                | 0       | 0      | 0
  number processed with summation       | 29676   | 2237   | 24
  number processed with profile fitting | 24069   | 1812   | 11
  number failed in background modelling | 85      | 29     | 0
  number failed in summation            | 5147    | 369    | 23
  number failed in profile fitting      | 10754   | 794    | 36
  ibg                                   | 17.79   | 49.75  | 4.36
  i/sigi (summation)                    | 24.71   | 139.98 | 1.44
  i/sigi (profile fitting)              | 32.19   | 182.76 | 1.77
  cc prf                                | 0.94    | 0.84   | 0.96
  cc_pearson sum/prf                    | 0.67    | 0.60   | 0.95
  cc_spearman sum/prf                   | 0.97    | 0.97   | 0.70
  ---------------------------------------------------------------

Symmetry, Scaling and Merging
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Again, we can proceed as standard, with the programs handling the multiple
lattices found in the datafiles::

  dials.symmetry integrated.expt integrated.refl

::

  Scoring all possible sub-groups

  ---------------------------------------------------------------------------------------------
  Patterson group       Likelihood  NetZcc  Zcc+   Zcc-   CC     CC-    delta  Reindex operator
  ---------------------------------------------------------------------------------------------
  P m m m          ***  0.988        9.87    9.87   0.00   0.99   0.00  0.0    a,b,c
  P 1 2/m 1             0.004        0.13    9.93   9.80   1.00   0.98  0.0    -b,-a,-c
  P 1 2/m 1             0.004        0.10    9.92   9.81   1.00   0.98  0.0    -a,-c,-b
  P 1 2/m 1             0.004        0.03    9.88   9.85   1.00   0.99  0.0    a,b,c
  P -1                  0.000        0.17    9.99   9.82   1.00   0.98  0.0    a,b,c
  ---------------------------------------------------------------------------------------------

  Best solution: P m m m
  Unit cell: (54.1104, 58.2822, 66.5198, 90, 90, 90)
  Reindex operator: a,b,c
  Laue group probability: 0.988
  Laue group confidence: 0.986

  ...

  Laue group: P m m m
  ---------------------------------------------------------------------------------------------------------------
  | Screw axis | Score | No. present | No. absent | <I> present | <I> absent | <I/sig> present | <I/sig> absent |
  ---------------------------------------------------------------------------------------------------------------
  | 21a        | 0.000 | 0           | 0          | 0.000       | 0.000      | 0.000           | 0.000          |
  | 21b        | 1.000 | 9           | 10         | 16592.455   | 13.316     | 141.588         | 0.412          |
  | 21c        | 1.000 | 12          | 13         | 28257.845   | 2.137      | 150.792         | 0.055          |
  ---------------------------------------------------------------------------------------------------------------
  ------------------------
  | Space group | score  |
  ------------------------
  | P 2 2 2     | 0.0000 |
  | P 2 2 21    | 0.0000 |
  | P 2 21 2    | 0.0000 |
  | P 21 2 2    | 0.0000 |
  | P 21 21 2   | 0.0000 |
  | P 21 2 21   | 0.0000 |
  | P 2 21 21   | 1.0000 |
  | P 21 21 21  | 0.0000 |
  ------------------------
  Recommended space group: P 2 21 21

The symmetry analysis suggested space group P 2 21 21, however it is worth
noting that no reflections were available to test the 21a screw axis, so this
possibility should also be tested during structure solution.

Next we scale the data and inspect the results from the log output or the
:samp:`dials.scale.html` generated html report::

  dials.scale symmetrized.expt symmetrized.refl

::

                       ----------Merging statistics----------

  Resolution: 28.89 - 1.06
  Observations: 46152
  Unique reflections: 34551
  Redundancy: 1.3
  Completeness: 36.52%
  Mean intensity: 2288.5
  Mean I/sigma(I): 15.4
  R-merge: 0.037
  R-meas:  0.051
  R-pim:   0.035


  Statistics by resolution bin:
  d_max  d_min   #obs  #uniq   mult.  %comp       <I>  <I/sI>    r_mrg   r_meas    r_pim   cc1/2   cc_ano
  28.90   2.89   3375   2366    1.43  46.98   15165.9    43.1    0.030    0.042    0.028   0.993*   0.725
   2.89   2.29   3400   2285    1.49  47.13    5801.2    36.6    0.031    0.042    0.028   0.992*  -0.239
   2.29   2.00   3369   2270    1.48  47.57    4059.8    31.6    0.035    0.046    0.030   0.993*  -0.235
   2.00   1.82   3297   2143    1.54  44.86    2512.9    26.2    0.039    0.053    0.036   0.986*  -0.064
   1.82   1.69   3226   2255    1.43  47.61    1506.1    19.6    0.049    0.066    0.045   0.983*   1.000
   1.69   1.59   3246   2372    1.37  50.22    1133.3    15.7    0.059    0.081    0.055   0.978*   0.000
   1.59   1.51   3229   2387    1.35  50.31     840.0    12.4    0.072    0.100    0.068   0.969*  -0.345
   1.51   1.45   3289   2481    1.33  52.63     659.5    10.0    0.081    0.112    0.077   0.962*  -0.561
   1.45   1.39   3177   2494    1.27  53.10     499.8     7.7    0.082    0.114    0.079   0.959*  -1.000
   1.39   1.34   3082   2480    1.24  52.71     409.9     6.3    0.095    0.132    0.091   0.937*   0.000
   1.34   1.30   3207   2540    1.26  54.03     371.5     5.6    0.089    0.124    0.085   0.973*   0.000
   1.30   1.26   2878   2306    1.25  49.25     331.8     4.9    0.097    0.136    0.095   0.956*   0.000
   1.26   1.23   2272   1845    1.23  39.23     313.8     4.6    0.077    0.108    0.076   0.977*   0.000
   1.23   1.20   1632   1293    1.26  27.46     286.4     4.1    0.075    0.105    0.074   0.979*   0.000
   1.20   1.17   1217    997    1.22  21.36     267.1     3.8    0.072    0.102    0.072   0.982*   0.000
   1.17   1.15    893    775    1.15  16.66     235.7     3.3    0.069    0.097    0.069   0.981*   0.000
   1.15   1.12    651    589    1.11  12.60     211.3     2.8    0.074    0.104    0.074   0.981*   0.000
   1.12   1.10    425    400    1.06   8.55     166.5     2.2    0.134    0.189    0.134   0.943*   0.000
   1.10   1.08    237    225    1.05   4.83     171.0     2.2    0.125    0.176    0.125   0.912*   0.000
   1.08   1.06     50     48    1.04   1.02     151.9     1.9    0.197    0.278    0.197   1.000   0.000
  28.89   1.06  46152  34551    1.34  36.52    2288.5    15.4    0.037    0.051    0.035   0.995*   0.445*


If required, we can rerun scaling with a resolution limit using the option
:samp:`d_min=`, however in this case the CC1/2 and <I/sI> are reasonable
to the highest resolution measured. The "Analysis by image number" plots
in the :samp:`dials.scale.html` report also indicate that both datasets are of similar
quality.

Once we are happy with the scaled dataset, a merged MTZ file can be generated::

  dials.merge scaled.expt scaled.refl

::

  Writing reflections to merged.mtz
  Title: From dials.merge
  Space group symbol from file: P22121
  Space group number from file: 18
  Space group from matrices: P 2 21 21 (No. 18)
  Point group symbol from file: 222
  Number of crystals: 1
  Number of Miller indices: 32658
  Resolution range: 28.8871 1.06499
  History:
    From DIALS 2.dev.1041-gf88516da7, run on 2019-10-28 at 15:40:44 GMT
  Crystal 1:
    Name: XTAL
    Project: AUTOMATIC
    Id: 1
    Unit cell: (54.1104, 58.2822, 66.5198, 90, 90, 90)
    Number of datasets: 1
    Dataset 1:
      Name: NATIVE
      Id: 1
      Wavelength: 0.97949
      Number of columns: 17
      label    #valid  %valid     min       max type
      H         32658 100.00%    0.00     37.00 H: index h,k,l
      K         32658 100.00%    0.00     53.00 H: index h,k,l
      L         32658 100.00%    0.00     60.00 H: index h,k,l
      IMEAN     32658 100.00% -109.31 213880.50 J: intensity
      SIGIMEAN  32658 100.00%    5.93   5046.71 Q: standard deviation
      I(+)      24848  76.09% -109.31 213880.50 K: I(+) or I(-)
      SIGI(+)   24848  76.09%    5.93   5046.71 M: standard deviation
      I(-)      12695  38.87%  -93.28 180658.78 K: I(+) or I(-)
      SIGI(-)   12695  38.87%   11.23   4263.08 M: standard deviation
      N(+)      24848  76.09%    2.00      4.00 I: integer
      N(-)      12695  38.87%    4.00      4.00 I: integer
      F         32658 100.00%    2.08    460.34 F: amplitude
      SIGF      32658 100.00%    0.25     27.28 Q: standard deviation
      F(+)      24848  76.09%    2.08    460.34 G: F(+) or F(-)
      SIGF(+)   24848  76.09%    0.25      5.48 L: standard deviation
      F(-)      12695  38.87%    2.51    421.81 G: F(+) or F(-)
      SIGF(-)   12695  38.87%    0.26      5.05 L: standard deviation

This program also performs a truncation, giving a set of merged intensities and
strictly-positive structure factors (Fs) suitable for downstream structure solution.
In the mtz file, both lattices are combined to give a single dataset. To generate
individual mtz files, one could first use :samp:`dials.split_experiments` before
merging with :samp:`dials.merge`::

  dials.split_experiments scaled.expt scaled.refl

::

  input {
    experiments = scaled.expt
    reflections = scaled.refl
  }

  Saving experiment 0 to split_0.expt
  Saving reflections for experiment 0 to split_0.refl
  Saving experiment 1 to split_1.expt
  Saving reflections for experiment 1 to split_1.refl

It is worth noting that the splitting of experiments can be performed on the
multi-lattice datafiles at any point during the processing if desired.

Refining multi-tile detector metrology with DIALS
=================================================

Introduction
------------

At the end of the :doc:`advanced_tutorial`, we showed plots from
:program:`dials.analyse_output`, including images of the positional
residuals as a function of position on the detector face. There were
clearly some systematic effects, suggesting whole-tile shifts or rotations.
Generally, these are small (typically much less than the pixel size) and are
expected. DECTRIS provides calibration tables from factory metrology
experiments, which XDS can use. Currently we are not applying the correction
tables to reflection positions in DIALS. In some situations we may have to work
with a much less well calibrated or characterised detector than the Pilatus
instruments installed at Diamond. Nevertheless, the multi-panel detector models
in DIALS allow us to *discover* the appropriate corrections assuming we have
good diffraction data to refine against.

Here we have access to a thaumatin dataset collected with low dose to avoid
radiation damage and with the detector 400 mm from the sample to ensure the
coverage of reflections extends right out to the corner of the images.

Preparing for refinement
------------------------

Usually multi-panel detectors like the Pilatus P6M are treated as if they were
single panels (with dead regions) in DIALS. Clearly, if we are to make corrections
to the relative positions and orientations of the panels, we need a more
detailed detector model in which the panels are treated separately. Currently
there isn't really a neat way to switch between the single and multiple panel
models, but until we have a native mechanism for doing this in DIALS there is
a hidden workaround that does this. We just need to set an environment variable::

  export P6M_60_PANEL=1

With that in place, we start as usual by importing the dataset::

  dials.import Thaum_M10S15_3_*.cbf

The fact that the multiple panel detector model was used is clear from the
output of::

  dials.show_models datablock.json

We can now inspect images from the dataset with::

  dials.image_viewer datablock.json

However, the dataset is at such low dose that it can be difficult to find the
spots. It might be helpful to us a more mature image viewer such as ADXV or
Albula for this.

We have to run spot finding throughout the dataset. It is a 360 degree sweep
so this will take a few minutes. We can use more processes to move a little
quicker::

  dials.find_spots datablock.json nproc=8 min_spot_size=3

Here we chose to set ``min_spot_size=3``, which overrides the default of 6 used
for this detector model. We did this because otherwise this weakly diffracting
dataset produces rather few strong spots, whereas we want as much coverage as
we can get, right out to the corners


Now to index the data. Although it is a well-diffracting crystal, running indexing
with defaults finds an approximate supercell with the *a* and *b* lengths slighly
more than doubled. FIXME WHAT CAN WE DO HERE? We don't worry too much about that
though, because if we pass in the known unit cell then it works just fine. We
also chose to apply tetragonal symmetry immediately::

  dials.index datablock.json strong.pickle space_group="P 4" unit_cell="58.7 58.7 151.59 90 90 90"

The output of refinement in the highest resolution macrocycle is as follows::

  ################################################################################
  Starting refinement (macro-cycle 5)
  ################################################################################


    Summary statistics for observations matched to predictions:
    --------------------------------------------------------------------------
    |                   | Min     | Q1       | Med        | Q3      | Max    |
    --------------------------------------------------------------------------
    | Xc - Xo (mm)      | -2.492  | -0.04261 | -0.0004832 | 0.0424  | 2.131  |
    | Yc - Yo (mm)      | -1.697  | -0.03481 | 0.001733   | 0.03682 | 1.656  |
    | Phic - Phio (deg) | -0.9645 | -0.06454 | -0.0007655 | 0.06379 | 0.8353 |
    | X weights         | 80.55   | 127.8    | 131.8      | 133.9   | 135.2  |
    | Y weights         | 77.02   | 127.6    | 131.8      | 134     | 135.2  |
    | Phi weights       | 26.58   | 32.59    | 32.65      | 32.65   | 32.65  |
    --------------------------------------------------------------------------


    Refinement steps:
    -------------------------------------------------
    | Step | Nref  | RMSD_X   | RMSD_Y   | RMSD_Phi |
    |      |       | (mm)     | (mm)     | (deg)    |
    -------------------------------------------------
    | 0    | 18025 | 0.067273 | 0.063827 | 0.086195 |
    | 1    | 18025 | 0.067369 | 0.063706 | 0.086197 |
    | 2    | 18025 | 0.067439 | 0.063629 | 0.086155 |
    | 3    | 18025 | 0.067531 | 0.063529 | 0.086115 |
    | 4    | 18025 | 0.067587 | 0.063456 | 0.086122 |
    | 5    | 18025 | 0.067607 | 0.063419 | 0.086155 |
    | 6    | 18025 | 0.067617 | 0.063401 | 0.086176 |
    | 7    | 18025 | 0.06762  | 0.063395 | 0.086181 |
    -------------------------------------------------
    RMSD no longer decreasing

    RMSDs by experiment:
    ----------------------------------------------
    | Exp | Nref  | RMSD_X  | RMSD_Y  | RMSD_Z   |
    |     |       | (px)    | (px)    | (images) |
    ----------------------------------------------
    | 0   | 18025 | 0.39314 | 0.36858 | 0.24623  |
    ----------------------------------------------

    RMSDs by panel:
    -----------------------------------------------
    | Panel | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
    |       |      | (px)    | (px)    | (images) |
    -----------------------------------------------
    | 0     | 10   | 0.53322 | 0.82846 | 0.14379  |
    | 1     | 140  | 0.37822 | 0.43153 | 0.26887  |
    | 2     | 297  | 0.41585 | 0.40848 | 0.26118  |
    | 3     | 146  | 0.5563  | 0.40867 | 0.27493  |
    | 4     | 17   | 0.49115 | 0.43483 | 0.21656  |
    | 5     | 44   | 0.41713 | 0.46071 | 0.27804  |
    | 6     | 405  | 0.44932 | 0.51832 | 0.26569  |
    | 7     | 758  | 0.33153 | 0.28766 | 0.25461  |
    | 8     | 478  | 0.38555 | 0.40525 | 0.26244  |
    | 9     | 96   | 0.47338 | 0.45054 | 0.24652  |
    | 10    | 152  | 0.40191 | 0.959   | 0.2976   |
    | 11    | 701  | 0.31004 | 0.3035  | 0.24597  |
    | 12    | 1128 | 0.37366 | 0.24303 | 0.24016  |
    | 13    | 802  | 0.33404 | 0.32094 | 0.24614  |
    | 14    | 201  | 0.50746 | 0.47124 | 0.27169  |
    | 15    | 231  | 0.40185 | 0.56122 | 0.3029   |
    | 16    | 745  | 0.21926 | 0.30075 | 0.21426  |
    | 17    | 831  | 0.17864 | 0.17706 | 0.2176   |
    | 18    | 778  | 0.36795 | 0.21167 | 0.20696  |
    | 19    | 269  | 0.44757 | 0.40028 | 0.23809  |
    | 20    | 205  | 0.48078 | 0.65283 | 0.34241  |
    | 21    | 467  | 0.37577 | 0.35789 | 0.20992  |
    | 22    | 370  | 0.22329 | 0.29697 | 0.20943  |
    | 23    | 447  | 0.36708 | 0.22722 | 0.20176  |
    | 24    | 231  | 0.28466 | 0.52795 | 0.31729  |
    | 25    | 3    | 0.36377 | 1.1547  | 0.70452  |
    | 26    | 28   | 0.14842 | 0.29869 | 0.23437  |
    | 27    | 22   | 0.13375 | 0.2666  | 0.17175  |
    | 28    | 17   | 0.1222  | 0.14108 | 0.17261  |
    | 29    | 4    | 0.40725 | 0.85708 | 0.42044  |
    | 30    | 145  | 0.60805 | 0.57195 | 0.43388  |
    | 31    | 281  | 0.18353 | 0.21351 | 0.2061   |
    | 32    | 250  | 0.19984 | 0.15417 | 0.1993   |
    | 33    | 294  | 0.15006 | 0.27665 | 0.21823  |
    | 34    | 158  | 0.44084 | 0.50436 | 0.38365  |
    | 35    | 203  | 0.86526 | 0.47248 | 0.32568  |
    | 36    | 525  | 0.54992 | 0.31614 | 0.19036  |
    | 37    | 521  | 0.22611 | 0.12494 | 0.21894  |
    | 38    | 624  | 0.23105 | 0.20793 | 0.22246  |
    | 39    | 259  | 0.43408 | 0.43312 | 0.2724   |
    | 40    | 134  | 0.55289 | 0.3935  | 0.26933  |
    | 41    | 557  | 0.29543 | 0.24331 | 0.23127  |
    | 42    | 771  | 0.46252 | 0.19655 | 0.2348   |
    | 43    | 658  | 0.24901 | 0.31232 | 0.24217  |
    | 44    | 199  | 0.44994 | 0.43547 | 0.28695  |
    | 45    | 54   | 0.5679  | 0.44719 | 0.23983  |
    | 46    | 333  | 0.44117 | 0.44019 | 0.24586  |
    | 47    | 592  | 0.3056  | 0.44653 | 0.24795  |
    | 48    | 411  | 0.67885 | 0.39219 | 0.2588   |
    | 49    | 93   | 0.45355 | 0.5659  | 0.2611   |
    | 50    | 5    | 1.0009  | 0.46232 | 0.21423  |
    | 51    | 162  | 0.72162 | 0.36619 | 0.23884  |
    | 52    | 324  | 0.51574 | 0.48569 | 0.2488   |
    | 53    | 224  | 0.41145 | 0.40151 | 0.2671   |
    | 54    | 16   | 0.5362  | 0.50568 | 0.21454  |
    | 55    | 1    | 0.85351 | 0.57076 | 0.24389  |
    | 56    | 37   | 0.72428 | 0.38373 | 0.27121  |
    | 57    | 106  | 0.80208 | 0.98469 | 0.27939  |
    | 58    | 61   | 0.5144  | 1.3281  | 0.30346  |
    | 59    | 4    | 0.33371 | 0.50058 | 0.24425  |
    -----------------------------------------------
    Final refined crystal models:
    model 1 (192715 reflections):
    Crystal:
        Unit cell: (57.834, 57.834, 150.022, 90.000, 90.000, 90.000)
        Space group: P 4
        U matrix:  {{ 0.4122, -0.9018,  0.1299},
                    { 0.2361, -0.0320, -0.9712},
                    { 0.8800,  0.4310,  0.1997}}
        B matrix:  {{ 0.0173,  0.0000,  0.0000},
                    { 0.0000,  0.0173,  0.0000},
                    { 0.0000,  0.0000,  0.0067}}
        A = UB:    {{ 0.0071, -0.0156,  0.0009},
                    { 0.0041, -0.0006, -0.0065},
                    { 0.0152,  0.0075,  0.0013}}


This refinement was performed moving all the panels as a rigid block, as usual.
With overall positional RMSDs within 40% of the pixel size and a
quarter of the image width in :math:`\phi` we can see straight away that we are
dealing with a good quality
dataset. There are a few outliers of well over 1 mm on the detector surface and nearly
1 degree in :math:`\phi` though, which we would prefer not to include in
refinement. The outliers are not as bad if we had kept :samp:`min_spot_size=6`,
but the detector coverage is worse in that case. At the moment we see that
coverage of reflections on the outer panels is rather low, but so far
we let refinement take a random subset of the data in order to index quicker,
so there's no need to worry about that yet.

Now we will refine the detector as a rigid block again, turning on outlier
rejection and requesting to use all reflections to get the best we can out
of the dataset::

  dials.refine indexed.pickle experiments.json do_outlier_rejection=true use_all_reflections=true

Here is the output::

  The following parameters have been modified:

  refinement {
    reflections {
      use_all_reflections = true
      do_outlier_rejection = true
    }
  }
  input {
    experiments = experiments.json
    reflections = indexed.pickle
  }

  Configuring refiner

  Summary statistics for observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -2.487  | -0.0424  | -0.0001972 | 0.04293 | 2.128  |
  | Yc - Yo (mm)      | -1.697  | -0.03456 | 0.001844   | 0.03669 | 1.654  |
  | Phic - Phio (deg) | -0.9622 | -0.06588 | -0.001969  | 0.06271 | 0.8323 |
  | X weights         | 80.55   | 127.8    | 131.8      | 133.9   | 135.2  |
  | Y weights         | 77.02   | 127.6    | 131.8      | 134     | 135.2  |
  | Phi weights       | 26.58   | 32.59    | 32.65      | 32.65   | 32.65  |
  --------------------------------------------------------------------------

  5849 reflections have been rejected as outliers

  Summary statistics for observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.2933 | -0.04194 | -0.0004627 | 0.04169 | 0.3725 |
  | Yc - Yo (mm)      | -0.2991 | -0.03369 | 0.001993   | 0.03605 | 0.339  |
  | Phic - Phio (deg) | -0.3981 | -0.06553 | -0.002115  | 0.06205 | 0.4137 |
  | X weights         | 80.55   | 127.9    | 131.9      | 134     | 135.2  |
  | Y weights         | 77.02   | 127.8    | 131.8      | 134     | 135.2  |
  | Phi weights       | 26.58   | 32.59    | 32.65      | 32.65   | 32.65  |
  --------------------------------------------------------------------------

  Performing refinement...

  Refinement steps:
  --------------------------------------------------
  | Step | Nref   | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |        | (mm)     | (mm)     | (deg)    |
  --------------------------------------------------
  | 0    | 181640 | 0.064502 | 0.056782 | 0.083589 |
  | 1    | 181640 | 0.064332 | 0.056888 | 0.083668 |
  | 2    | 181640 | 0.064055 | 0.057074 | 0.083787 |
  | 3    | 181640 | 0.063711 | 0.057327 | 0.083944 |
  | 4    | 181640 | 0.0635   | 0.057483 | 0.084045 |
  | 5    | 181640 | 0.063432 | 0.057512 | 0.084096 |
  | 6    | 181640 | 0.063413 | 0.057511 | 0.084131 |
  | 7    | 181640 | 0.063408 | 0.057511 | 0.084143 |
  | 8    | 181640 | 0.063408 | 0.057511 | 0.084144 |
  --------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  -----------------------------------------------
  | Exp | Nref   | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |        | (px)    | (px)    | (images) |
  -----------------------------------------------
  | 0   | 181640 | 0.36865 | 0.33436 | 0.24041  |
  -----------------------------------------------

  RMSDs by panel:
  ------------------------------------------------
  | Panel | Nref | RMSD_X   | RMSD_Y  | RMSD_Z   |
  |       |      | (px)     | (px)    | (images) |
  ------------------------------------------------
  | 0     | 63   | 0.43761  | 0.59343 | 0.22242  |
  | 1     | 1363 | 0.34575  | 0.46174 | 0.26593  |
  | 2     | 3076 | 0.41044  | 0.40681 | 0.26674  |
  | 3     | 1630 | 0.46819  | 0.40144 | 0.25618  |
  | 4     | 112  | 0.39168  | 0.43874 | 0.2404   |
  | 5     | 394  | 0.36022  | 0.51907 | 0.25565  |
  | 6     | 3489 | 0.32483  | 0.48727 | 0.2574   |
  | 7     | 6172 | 0.30681  | 0.27397 | 0.2534   |
  | 8     | 3908 | 0.36528  | 0.38764 | 0.25624  |
  | 9     | 680  | 0.45094  | 0.45795 | 0.25022  |
  | 10    | 1178 | 0.32726  | 0.5171  | 0.26474  |
  | 11    | 5690 | 0.2811   | 0.29208 | 0.24445  |
  | 12    | 9385 | 0.36543  | 0.23363 | 0.23893  |
  | 13    | 6878 | 0.31472  | 0.31484 | 0.23731  |
  | 14    | 1749 | 0.44638  | 0.47916 | 0.25978  |
  | 15    | 1870 | 0.28936  | 0.47177 | 0.26865  |
  | 16    | 6328 | 0.19245  | 0.29027 | 0.21606  |
  | 17    | 7599 | 0.17261  | 0.17213 | 0.22109  |
  | 18    | 7072 | 0.36515  | 0.19168 | 0.20884  |
  | 19    | 2609 | 0.38794  | 0.43612 | 0.24767  |
  | 20    | 1806 | 0.41482  | 0.58552 | 0.33021  |
  | 21    | 4247 | 0.35742  | 0.34731 | 0.19571  |
  | 22    | 3558 | 0.21577  | 0.30826 | 0.20734  |
  | 23    | 4103 | 0.36506  | 0.19519 | 0.1819   |
  | 24    | 2270 | 0.24861  | 0.4954  | 0.2843   |
  | 25    | 79   | 0.32992  | 0.70354 | 0.50852  |
  | 26    | 210  | 0.15254  | 0.24686 | 0.24963  |
  | 27    | 201  | 0.13616  | 0.18361 | 0.17548  |
  | 28    | 167  | 0.095196 | 0.22282 | 0.22815  |
  | 29    | 36   | 0.35193  | 0.50843 | 0.37494  |
  | 30    | 1283 | 0.57562  | 0.37178 | 0.33153  |
  | 31    | 2865 | 0.15601  | 0.20173 | 0.18891  |
  | 32    | 2535 | 0.20577  | 0.14842 | 0.19661  |
  | 33    | 2982 | 0.13275  | 0.31198 | 0.20958  |
  | 34    | 1754 | 0.40887  | 0.44466 | 0.3884   |
  | 35    | 1934 | 0.60375  | 0.36473 | 0.2557   |
  | 36    | 5875 | 0.55012  | 0.29471 | 0.19304  |
  | 37    | 5992 | 0.21294  | 0.11923 | 0.21505  |
  | 38    | 6589 | 0.22515  | 0.21212 | 0.21398  |
  | 39    | 2763 | 0.38378  | 0.41087 | 0.28153  |
  | 40    | 1512 | 0.50435  | 0.37702 | 0.25212  |
  | 41    | 6476 | 0.2759   | 0.22112 | 0.22689  |
  | 42    | 9212 | 0.40231  | 0.18054 | 0.22522  |
  | 43    | 7521 | 0.22574  | 0.29632 | 0.23973  |
  | 44    | 2392 | 0.39211  | 0.42886 | 0.2795   |
  | 45    | 611  | 0.56255  | 0.48595 | 0.26298  |
  | 46    | 4379 | 0.41975  | 0.4279  | 0.25309  |
  | 47    | 7749 | 0.28495  | 0.43514 | 0.24637  |
  | 48    | 5470 | 0.66489  | 0.37512 | 0.25272  |
  | 49    | 1210 | 0.42496  | 0.53357 | 0.27586  |
  | 50    | 127  | 0.89456  | 0.46188 | 0.24677  |
  | 51    | 2235 | 0.7014   | 0.33806 | 0.2579   |
  | 52    | 4379 | 0.39793  | 0.38196 | 0.25756  |
  | 53    | 2757 | 0.43811  | 0.39934 | 0.26732  |
  | 54    | 327  | 0.53208  | 0.50801 | 0.28399  |
  | 56    | 457  | 0.58832  | 0.43364 | 0.26501  |
  | 57    | 1507 | 0.69138  | 0.66988 | 0.25896  |
  | 58    | 786  | 0.52859  | 0.50638 | 0.27265  |
  | 59    | 39   | 0.45783  | 0.65016 | 0.25962  |
  ------------------------------------------------
  Saving refined experiments to refined_experiments.json

Outlier rejection has cleaned up the positional residuals so now the greatest
deviation is within 0.4 mm of the predicted position. The angular extreme is
now just over 0.4 degrees. Coverage of the outer and central panels (where
reflections are occluded by the backstop or thrown away for being too close
to the spindle) is still a little low. Notably, panel 55 (a corner panel) is
completely missing.

Acknowledgements
----------------

Dave Hall (Diamond Light Source) for collecting the data


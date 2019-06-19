Data files
==========

The DIALS programs read and write three main types of data file for storing the
experimental geometry, image data and processed reflection data. These are
summarized in the following table and described in more detail in the sections
below.

+------------------+-----------------------------------------------------+
| File type        | Contains                                            |
+==================+=====================================================+
| Datablock        | Experimental geometry and image data                |
+------------------+-----------------------------------------------------+
| Experiment list  | Experimental geometry (plus crystal) and image data |
+------------------+-----------------------------------------------------+
| Reflection table | Processed reflection data                           |
+------------------+-----------------------------------------------------+

The data files input and output from the main dials programs are described
below.

+------------------+-------------------------------+-------------------------------+
| Program          | Reads                         | Writes                        |
+==================+===============================+===============================+
| dials.import     | N/A                           | datablock.expt                |
+------------------+-------------------------------+-------------------------------+
| dials.find_spots | datablock.expt                | strong.refl                   |
+------------------+-------------------------------+-------------------------------+
| dials.index      | | datablock.expt              | | indexed.expt                |
|                  | | strong.refl                 | | indexed.refl                |
+------------------+-------------------------------+-------------------------------+
| dials.refine     | | indexed.expt                | | refined.expt                |
|                  | | indexed.refl                | | refined.refl                |
+------------------+-------------------------------+-------------------------------+
| dials.integrate  | | refined.expt                | | integrated.expt             |
|                  | | refined.refl                | | integrated.refl             |
|                  | |                             | | profile_model.phil          |
+------------------+-------------------------------+-------------------------------+
| dials.scale      | | integrated.expt             | | scaled.expt                 |
|                  | | integrated.refl             | | scaled.refl                 |
+------------------+-------------------------------+-------------------------------+

.. _datablock-json:

Datablock files
---------------

The datablock file is stored as a JSON file in ascii format. Whilst being human
readable (and editable), it is generally recommended to let dials.import
generate the datablock file. The file contains the location of any imported
imagesets and the initial experimental models (i.e. beam, detector, goniometer
and scan) for the sets of sweeps or still images. It also encodes the
relationship between models if multiple sweeps or sets of stills are imported.
An example of a short file is shown below.

.. code-block:: js

  [
    {
      "__id__": "DataBlock",
      "imageset": [
        {
          "__id__": "ImageSweep",
          "template": "image_####.cbf",
          "beam": 0,
          "detector": 0,
          "goniometer": 0,
          "scan": 0
        }
      ],
      "beam": [
        {
          "direction": [
            0.0,
            0.0,
            1.0
          ],
          "polarization_normal": [
            0.0,
            1.0,
            0.0
          ],
          "divergence": 0.0,
          "polarization_fraction": 0.999,
          "sigma_divergence": 0.0,
          "wavelength": 0.9795
        }
      ],
      "detector": [
        {
          "panels": [
            {
              "origin": [
                -212.47848,
                220.00176,
                -190.17999999999998
              ],
              "fast_axis": [
                1.0,
                0.0,
                0.0
              ],
              "name": "Panel",
              "slow_axis": [
                0.0,
                -1.0,
                0.0
              ],
              "trusted_range": [
                -1.0,
                495976.0
              ],
              "image_size": [
                2463,
                2527
              ],
              "px_mm_strategy": {
                "type": "ParallaxCorrectedPxMmStrategy"
              },
              "type": "SENSOR_PAD",
              "pixel_size": [
                0.17200000000000001,
                0.17200000000000001
              ]
            }
          ]
        }
      ],
      "goniometer": [
        {
          "fixed_rotation": [
            1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            1.0
          ],
          "rotation_axis": [
            1.0,
            0.0,
            0.0
          ]
        }
      ],
      "scan": [
        {
          "exposure_time": [
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2,
            0.2
          ],
          "epochs": [
            1360324992.0,
            1360324992.0,
            1360324993.0,
            1360324993.0,
            1360324993.0,
            1360324993.0,
            1360324993.0,
            1360324994.0,
            1360324994.0
          ],
          "image_range": [
            1,
            9
          ],
          "oscillation": [
            0.0,
            0.2
          ]
        }
      ]
    }
  ]

.. _experiments_json:

Experiment list files
---------------------

The experiment list file is stored as a JSON file in ascii format. Whilst being human
readable (and editable), editing the file directly is generally not recommended.
The file contains the location of any imported imagesets and the same
experimental models contained in the datablock file (i.e. beam, detector,
goniometer and scan). It also contains the crystal models found during indexing.
The data is modelled by a set of experiments which each contain a single set of
models. Experiments can share models. For example two experiments may share
detector models. This allows, for example, joint refinement of experiments.

An example of a short file is shown below.

.. code-block:: js

  {
    "__id__": "ExperimentList",
    "experiment": [
      {
        "__id__": "Experiment",
        "beam": 0,
        "detector": 0,
        "goniometer": 0,
        "scan": 0,
        "crystal": 0,
        "imageset": 0
      }
    ],
    "imageset": [
      {
        "__id__": "ImageSweep",
        "template": "centroid_####.cbf"
      }
    ],
    "beam": [
      {
        "direction": [
          -0.007852057721998333,
          3.772524827250213e-14,
          0.9999691721195861
        ],
        "polarization_normal": [
          0.0,
          1.0,
          0.0
        ],
        "divergence": 0.0,
        "polarization_fraction": 0.999,
        "sigma_divergence": 0.058,
        "wavelength": 0.9795
      }
    ],
    "detector": [
      {
        "panels": [
          {
            "origin": [
              -211.53596470096178,
              219.45303890619488,
              -192.7062494437063
            ],
            "fast_axis": [
              0.9999551354884303,
              0.0021159302715049923,
              0.009233084500921031
            ],
            "name": "Panel",
            "slow_axis": [
              0.0021250002879257116,
              -0.999997269169901,
              -0.0009726389448611214
            ],
            "trusted_range": [
              -1.0,
              495976.0
            ],
            "image_size": [
              2463,
              2527
            ],
            "px_mm_strategy": {
              "type": "ParallaxCorrectedPxMmStrategy"
            },
            "type": "SENSOR_UNKNOWN",
            "pixel_size": [
              0.172,
              0.172
            ]
          }
        ]
      }
    ],
    "goniometer": [
      {
        "fixed_rotation": [
          1.0,
          0.0,
          0.0,
          0.0,
          1.0,
          0.0,
          0.0,
          0.0,
          1.0
        ],
        "rotation_axis": [
          1.0,
          -1.5919306617286774e-16,
          -6.904199434387693e-16
        ]
      }
    ],
    "scan": [
      {
        "exposure_time": [
          0.2,
          0.2,
          0.2,
          0.2,
          0.2,
          0.2,
          0.2,
          0.2,
          0.2
        ],
        "epochs": [
          1360324992.0,
          1360324992.0,
          1360324993.0,
          1360324993.0,
          1360324993.0,
          1360324993.0,
          1360324993.0,
          1360324994.0,
          1360324994.0
        ],
        "image_range": [
          1,
          9
        ],
        "oscillation": [
          0.0,
          0.2
        ]
      }
    ],
    "crystal": [
      {
        "__id__": "crystal",
        "real_space_a": [
          35.23781811553089,
          -7.600614003857873,
          22.077690418635804
        ],
        "real_space_b": [
          -22.657129890916668,
          -1.4698317405529955,
          35.65693038892429
        ],
        "real_space_c": [
          -5.295803077552249,
          -38.99952334925477,
          -4.972795822746061
        ],
        "space_group_hall_symbol": " P 4 2",
        "mosaicity": 0.157
      }
    ]
  }

.. _reflection_pickle:

Reflection files
----------------

The reflection files are saved in python's "pickle" format. This is a binary
format that is convenient for serializing python classes. The reflection files
will contain a table with some or all of the following columns.


+-------------------------------+------------------------------------------------------+
| Column                        | Description                                          |
+===============================+======================================================+
| flags                         | bit mask status flags                                |
+-------------------------------+------------------------------------------------------+
| id                            | experiment id                                        |
+-------------------------------+------------------------------------------------------+
| panel                         | the detector panel index                             |
+-------------------------------+------------------------------------------------------+
| miller_index                  | miller indices                                       |
+-------------------------------+------------------------------------------------------+
| entering                      | reflection entering/exiting                          |
+-------------------------------+------------------------------------------------------+
| s1                            | the diffracted beam vector                           |
+-------------------------------+------------------------------------------------------+
| xyzcal.mm                     | the predicted location (mm, mm, rad)                 |
+-------------------------------+------------------------------------------------------+
| xyzcal.px                     | the predicted location (px, px, frame)               |
+-------------------------------+------------------------------------------------------+
| ub_matrix                     | predicted crystal setting                            |
+-------------------------------+------------------------------------------------------+
| xyzobs.px.value               | centroid pixel position  (px, px, frame)             |
+-------------------------------+------------------------------------------------------+
| xyzobs.px.variance            | centroid pixel variance                              |
+-------------------------------+------------------------------------------------------+
| xyzobs.mm.value               | centroid millimetre position (mm, mm, rad)           |
+-------------------------------+------------------------------------------------------+
| xyzobs.mm.variance            | centroid millimetre variance                         |
+-------------------------------+------------------------------------------------------+
| rlp                           | reciprocal lattice point                             |
+-------------------------------+------------------------------------------------------+
| intensity.sum.value           | raw intensity value                                  |
+-------------------------------+------------------------------------------------------+
| intensity.sum.variance        | raw intensity variance                               |
+-------------------------------+------------------------------------------------------+
| intensity.prf.value           | profile fitted intensity value                       |
+-------------------------------+------------------------------------------------------+
| intensity.prf.variance        | profile fitted intensity variance                    |
+-------------------------------+------------------------------------------------------+
| | intensity.scale.value       | | intensity value used for scaling                   |
| |                             | | (without scale factor applied)                     |
+-------------------------------+------------------------------------------------------+
| intensity.scale.variance      | variance of intensity value used for scaling         |
+-------------------------------+------------------------------------------------------+
| inverse_scale_factor          | scale factor determined by scaling (divisory)        |
+-------------------------------+------------------------------------------------------+
| inverse_scale_factor_variance | variance of inverse scale factor                     |
+-------------------------------+------------------------------------------------------+
| lp                            | LP correction (multiplicative)                       |
+-------------------------------+------------------------------------------------------+
| qe                            | detector quantum efficiency correction (divisory)    |
+-------------------------------+------------------------------------------------------+
| profile.correlation           | correlation in profile fitting                       |
+-------------------------------+------------------------------------------------------+
| | partiality                  | | fraction of reflection measured                    |
| |                             | | (i.e. I\ :sub:`full` = I\ :sub:`sum`\ /partiality) |
+-------------------------------+------------------------------------------------------+
| bbox                          | bounding box                                         |
+-------------------------------+------------------------------------------------------+
| shoebox                       | shoebox data/mask/background struct                  |
+-------------------------------+------------------------------------------------------+

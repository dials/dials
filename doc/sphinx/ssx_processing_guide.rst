Experimental SSX processing guide
=================================

This is a guide on how to process synchrotron serial crystallography (SSX) data
with DIALS, using tools that are currently under development. These tools should
be considered experimental and subject to change \& improvement as the tools
become more widely tested and user feedback is taken on board.

Indexing SSX data with dev.dials.ssx_index
------------------------------------------

A sequence of SSX images can be imported and spotfinding can be run in the same
way as for rotation data::

    dials.import /data/images*cbf
    dials.find_spots imported.expt

Note that for simplicity in this example we are not using a reference geometry when importing.
These commands produce an :samp:`imported.expt` experiments file containing the experimental
metadata and a :samp:`strong.refl` reflections file containing the found spots.
Then to index the data, we can use::

    dev.dials.ssx_index strong.refl imported.expt

This program wraps a call to the dials.index program with options suitable
for processing still images. As with indexing regular sweeps, if the unit cell
and/or space group are known, providing them as input to the program gives a
greater chance of successful indexing. These can be provided as additional options,
for example in this format::

    dev.dials.ssx_index strong.refl imported.expt space_group=P4 unit_cell=60,60,85,90,90,90

For :samp:`dials.ssx_index`, if the unit cell is given, then indexing will be attempted
on each image with the "fft1d" algorithm, followed by the "real space grid search"
algorithm if the fft1d indexing is not successful. If the unit cell is not given,
only the fft1d algorithm is used. To specify the use of only one particular method,
the :samp:`method` option can be set to just one of :samp:`fft1d real_space_grid_search`.
As each image may contain diffraction from more than one crystal, another useful
option is to set the :samp:`max_lattices` parameter to greater than one, to
enable multiple crystals to be found on each image.

As the indexing of each image is an independent process, the processing will
be split across the available computing cores, significantly speeding up the
processing when using computing clusters/high performance machines. After
indexing has been attempted on all images, unit cell clustering is performed and
reported. This is particularly useful if the unit cell is not currently known,
as it can be used as a guide for assessing the crystal parameters::

    2 clusters:

    Cluster_id       N_xtals  Med_a         Med_b         Med_c         Med_alpha    Med_beta     Med_gamma   Delta(deg)
    22 in P1.
    cluster_1       22       96.38 (0.02 ) 96.39 (0.02 ) 96.42 (0.03 ) 90.02 (0.02) 90.02 (0.03) 90.02 (0.03)
        P m -3 m (No. 221)  96.40         96.40         96.40         90.00        90.00        90.00         0.038

    7 in P1.
    cluster_2       7        136.26(0.10 ) 136.32(0.04 ) 166.99(0.13 ) 90.01 (0.05) 90.01 (0.03) 119.94(29.64)
        P 6/m m m (No. 191)  136.33        136.33        166.99        90.00        90.00        120.00        0.077

Repeating indexing with the unit cell/space group of the most populous cluster
can then improve the overall indexing result.

To help with the interpretation of the indexing results for a large number of
images, a :samp:`dials.ssx_index.html` report is generated which contains plots
of useful statistics such as the number of spots indexed on each image, the distribution
of rmsd values and unit cell clustering analysis. This data can also be output to
json format for further analysis, by providing a filename to the option :samp:`output.json`.

Integrating SSX data with dev.dials.ssx_integrate
-------------------------------------------------

After indexing, the experimental models can be further refined with dials.refine,
or the indexing output can also be integrated directly.
To integrate the data, we can use::

    dev.dials.ssx_integrate indexed.expt indexed.refl

This program wraps a call to parts of the :samp:`dials.integrate` program,
using either the :samp:`stills` integrator or the :samp:`ellipsoid` integration algorithm.
The stills integrator is the default algorithm used for integration in
:samp:`dials.stills_process`. The ellipsoid algorithm refines the unit cell,
orientation and a 3D ellipsoidal mosaicity parameterisation for each crystal,
by assessing the pixel-intensity distribution of the strong spots::

    dev.dials.ssx_integrate indexed.refl indexed.expt algorithm=stills
    dev.dials.ssx_integrate indexed.refl indexed.expt algorithm=ellipsoid

Processing will be split across the available computing cores for performance.
During processing, data files will be created after each batch of crystals has
been processed. The size of the batch for saving data can be set with the
:samp:`batch_size` option. This creates numbered output files such as
:samp:`integrated_0.refl, integrated_0.expt, integrated_1.refl, integrated_1.expt` etc.
After all images have been integated, unit cell clustering is performed and
reported, as this will have changed compared to at the end of indexing if
using the ellipsoid integration algorithm.

To help with the interpretation of the integration results for a large number of
crystals, a :samp:`dials.ssx_integrate.html` report is generated which contains plots
of useful statistics such as the number of spots integrated on each image,
the modelled mosaicity values and unit cell clustering analysis. This data can
also be output to json format for further analysis, by providing a filename to
the option :samp:`output.json`.

Automated processing of SSX data with dev.xia2.ssx
--------------------------------------------------
An experimental processing pipeline is being developed to provide automated
processing of SSX data. To process SSX data, it is important to have an accurate
reference geometry for the experiment. Also, processing with the correct space
group and unit cell will increase the success rates in indexing and generally
improve the data integration.

Assessing the unit cell and space group
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If no unit cell or no space group is provided, this is the first thing that will be assessed.
To run the program, at a minimum, the path to images must be provided::

    dev.xia2.ssx images=~/data/myimages_*cbf

This will import the images, creating an :samp:`imported.expt` file in the :samp:`initial_import` folder.
Then, an assessment of the space group and unit cell will be performed on the first
1000 images. This is done by running spotfinding and indexing, and inspecting
the indexing clustering results. The clustering results after indexing indicate
which is the highest-symmetry possible, based on the unit cell dimensions. To control
which images are used for crystal assessment, there are the additional options
:samp:`assess_crystal.n_images` and :samp:`assess_crystals.images_to_use`.
Once a decision has been made on the unit cell and space group, the program can
then be rerun with a unit cell and space group of choice to continue with the data
processing::

    dev.xia2.ssx images=~/data/myimages_*cbf space_group=P1 unit_cell=95,95,97,90,90,90

Using a reference geometry
^^^^^^^^^^^^^^^^^^^^^^^^^^

For the reference geometry, this can be provided initially as a separate :samp:`refined.expt`
file, or determined as part of the data processing.::

    dev.xia2.ssx images=~/data/myimages_*cbf reference_geometry=~/data/refined.expt

If a reference geometry is supplied, the images will be imported with this as a
reference, and processing will proceed from that point. If no reference geometry
is supplied, a reference geometry will first be determined by performing a joint
refinement on the first 1000 images, with the results saved to the
:samp:`geometry_refinement` folder. There are additional options of
:samp:`geometry_refinement.n_images` and :samp:`geometry_refinement.images_to_use` to control
which images are used for the geometry refinement. Once this has been determined,
the images will be reimported with this as a reference geometry, and processing
will proceed. For subsequent processing, this can then be explicity provided as the
reference::

    dev.xia2.ssx images=~/data/myimages_*cbf
    dev.xia2.ssx images=~/data/myimages_*cbf reference_geometry=geometry_refinement/refined.expt

Main processing
^^^^^^^^^^^^^^^
With a given unit cell, space group and reference geometry, the whole of the
dataset is then integrated in batches of 1000 images, although this can be adjusted
with the :samp:`batch_size` option. A subfolder is created for each batch, where
spotfinding, indexing and integration are performed. This creates a set of
integrated reflection and experiment files which can be used for data reduction.



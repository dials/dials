+++++
Processing Sweeps with Missing Images
+++++

Most of DIALS assumes you have a single sweep of data, and in some circumstances deviations from this will cause problems. A clear example of this is having one or more "bad" images in your data set - simply removing them *will not* give the expected outcome i.e. processing proceeding otherwise smoothly.

Importing
=====

Importing the images with ``template=blah_####.cbf`` will not work, as this checks for continuous images, and it will be necessary to use ``allow_multiple_sweeps=true``. After this, finding spots and indexing work as usual.

Refinement
=====

Refinement *without* ``scan_varying=true`` will work fine, but with scan varying an error will occur - you therefore need to split the experiment into a number of blocks i.e.

``
dials.split_experiments indexed.pickle experiments.json
``

Then (ideally in a new directory) process each block as per the usual tutorial instructions. Finally, sort together all exported MTZ files with POINTLESS.

############################################################################
# DIALS tutorial script - associated with data doi:10.5281/zenodo.10271
############################################################################

# first import the data - this reads the image headers and creates a file which
# describes the sweeps found therein. here there should only be one sweep found
# data should be passed on the command-line e.g. ./tutorial.sh /path/to/data
# or linked to the directory where the script is run as "data"

if test -z $1; then
  data_directory=./data
  echo "Using ./data as import parameter"
else
  data_directory=$1
  echo "Using $1 (as entered by CLI)"
fi

if test -z $2; then
  nproc=$((`libtbx.show_number_of_processors` / 2))
else
  nproc=$2
fi

if test ! -f ${data_directory}/th_8_2_0001.cbf; then
  echo "Data not found in directory: ${data_directory} -"
  echo "please download from doi:10.5281/zenodo.10271 and pass on command-line"
  echo "or create softlink from data directory to ./data"
  exit
fi

dials.import ${data_directory}/th_8_2_0*.cbf

# find spots on the images - this will take a while! also accept small spots
# in this case as the data are good (N.B. to explore DIALS options you can
# run dials.parameters - this shows you the entire parameter set which is
# available for DIALS)

dials.find_spots datablock.json nproc=4

# index these found spots - in this case the crystal is thaumatin which is known
# to be tetragonal, so impose appropriate lattice constraints in the indexing
# and refinement. can also specify unit cell here and also apply different
# indexing algorithms

dials.index datablock.json strong.pickle


# Take the results of the P1 autoindexing and run refinement
# with all of the possible Bravais settings applied

dials.refine_bravais_settings experiments.json indexed.pickle

# Sometimes it may be necessary to reindex the indexed.pickle file output
# by dials.index. However, in this case as the change of basis operator to
# the chosen setting is the identity operator (a,b,c) this step is not needed.
# We run it anyway to demonstrate its use:

dials.reindex indexed.pickle change_of_basis_op=a,b,c

# During indexing we saw the presence of outliers that we would like to exclude
# from refinement, and we also used a subset of reflections. Now we will repeat
# using all indexed reflections in the dataset and with outlier rejection switched on.

dials.refine bravais_setting_9.json reindexed_reflections.pickle \
outlier.algorithm=tukey use_all_reflections=true

# We have done the best we can with a static model for the experiment.
# However, a better model for the crystal might allow small misset rotations
# to occur over the course of the scan. There are usually even small changes
# to the cell dimensions (typically resulting in a net increase in cell volume)
# caused by exposure to radiation during data collection. To account for both of
# these effects we can extend our parameterisation to obtain a
# smoothed ‘scan-varying’ model for both the crystal orientation and unit cell.
# To do this, we run a further refinement job starting from the output of the
# previous job:

dials.refine refined_experiments.json refined.pickle \
outlier.algorithm=tukey use_all_reflections=true  \
scan_varying=true output.experiments=sv_refined_experiments.json

# Integration:
# Next step reads the indexed reflections to determine strong reflections for profile
# fitting and integrates the data in refined_experiments.json, using the default
# background determination with no outlier rejection and XDS-style 3D profile
# fitting. These commands are most likely to change and can be viewed by running

dials.integrate sv_refined_experiments.json refined.pickle \
outlier.algorithm=null nproc=4

# finally export the integrated measurements in an MTZ file - this should be
# properly formatted for immediate use in pointless / aimless

dials.export integrated.pickle refined_experiments.json hklout=integrated.mtz

# and as if to make a point, here is all we need to do is to sort the data with
# pointless and then scale the data (ignoring anomalous differences) to 1.3A,
# and run ctruncate for the intensity analysis...

pointless hklin integrated.mtz hklout sorted.mtz 2>&1 >pointless.log || exit 1

aimless hklin sorted.mtz hklout scaled.mtz 2>&1 >aimless.log << eof || exit 1
resolution 1.3
anomalous off
eof

ctruncate -hklin scaled.mtz -hklout truncated.mtz \
-colin '/*/*/[IMEAN,SIGIMEAN]' 2>&1 >ctruncate.log || exit 1

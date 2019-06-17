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

dials.find_spots min_spot_size=3 datablock.expt nproc=$nproc

# index these found spots - in this case the crystal is thaumatin which is known
# to be tetragonal, so impose appropriate lattice constraints in the indexing
# and refinement. can also specify unit cell here and also apply different
# indexing algorithms

dials.index strong.refl datablock.expt

# Refining: If you do want to use a time varying model,
# you will need to rerun the refinement with this new model as

dials.refine indexed.expt indexed.refl scan_varying=true

# Integration:
# Next step reads the indexed reflections to determine strong reflections for profile
# fitting and integrates the data in refined.expt, using the default
# background determination with no outlier rejection and XDS-style 3D profile
# fitting. These commands are most likely to change and can be viewed by running

dials.integrate outlier.algorithm=null refined.expt indexed.refl

# finally export the integrated measurements in an MTZ file - this should be
# properly formatted for immediate use in pointless / aimless

dials.export integrated.refl integrated.expt mtz.hklout=integrated.mtz

# and as if to make a point, here is all we need to do is to sort the data with
# pointless and then scale the data (ignoring anomalous differences) to 1.3A,
# and run ctruncate for the intensity analysis...

pointless hklin integrated.mtz hklout sorted.mtz > pointless.log

aimless hklin sorted.mtz hklout scaled.mtz > aimless.log << EOF
resolution 1.3
anomalous off
EOF

ctruncate -hklin scaled.mtz -hklout truncated.mtz \
-colin '/*/*/[IMEAN,SIGIMEAN]' > ctruncate.log

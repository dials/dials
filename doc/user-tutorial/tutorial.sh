############################################################################
# DIALS tutorial script - associated with data doi:10.5281/zenodo.10271
############################################################################

# first import the data - this reads the image headers and creates a file which
# describes the sweeps found therein. here there should only be one sweep found
# data should be passed on the command-line e.g. ./tutorial.sh /path/to/data
# or linked to the directory where the script is run as "data"

if test -z $1; then
  data_directory=./data
else
  data_directory=$1
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

dials.find_spots min_spot_size=3 datablock.json nproc=$nproc

# index these found spots - in this case the crystal is thaumatin which is known
# to be tetragonal, so impose appropriate lattice constraints in the indexing
# and refinement. can also specify unit cell here and also apply different
# indexing algorithms

dials.index datablock.json strong.pickle \
  refinement.reflections.use_all_reflections=true

# refine the indexing solution in all Bravais settings consistent with the
# indexed unit cell. In this example we would continue processing using the
# bravais_setting_9.json, i.e. solution number 9

dials.refine_bravais_settings experiments.json indexed.pickle

# If necessary (i.e. if the cb_op for the chosen solution is not a,b,c)
# run dials.reindex on the indexed.pickle file to match the chosen solution

dials.reindex indexed.pickle change_of_basis_op=a,b,c

# run the refinement of the data - the indexing includes refinement so the result
# will be refined already with a static model - however here we would like to use
# a scan varying model so re-run the refinement (you should find that the R.M.S.
# deviations are a little lower following the scan varying refinement). Here we
# use an expert option "bin_size_fraction=0.0" to ensure the refinement runs
# to RMSD convergence rather than terminating early with a good-enough RMSD.

dials.refine bravais_setting_9.json reindexed_reflections.pickle \
  refinement.parameterisation.crystal.scan_varying=true \
  refinement.reflections.use_all_reflections=true \
  refinement.target.bin_size_fraction=0.0

# now run the integration - complex choices of algorithms are shown here in
# terms of the peak measurement (fitrs) and background determination
# methods. pass reference reflections from indexing in to determine the
# profile parameters...

dials.integrate outlier.algorithm=null profile.fitting=True \
  input.experiments=refined_experiments.json \
  input.reflections=reindexed_reflections.pickle \
  nproc=$nproc

# finally export the integrated measurements in an MTZ file - this should be
# properly formatted for immediate use in pointless / aimless

dials.export_mtz integrated.pickle refined_experiments.json hklout=integrated.mtz

# and as if to make a point, here is all we need to do is to sort the data with
# pointless and then scale the data (ignoring anomalous differences) to 1.3A,
# and run ctruncate for the intensity analysis...

pointless hklin integrated.mtz hklout sorted.mtz > pointless.log

aimless hklin sorted.mtz hklout scaled.mtz > aimless.log << eof
resolution 1.3
anomalous off
eof

ctruncate -hklin scaled.mtz -hklout truncated.mtz \
-colin '/*/*/[IMEAN,SIGIMEAN]' > ctruncate.log

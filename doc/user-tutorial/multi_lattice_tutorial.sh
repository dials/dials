############################################################################
# DIALS tutorial script - associated with data doi:10.5281/zenodo.10820
############################################################################

# first import the data - this reads the image headers and creates a file which
# describes the sweeps found therein. here there should only be one sweep found
# data should be passed on the command-line e.g. ./tutorial.sh /path/to/data
# or linked to the directory where the script is run as "data"
# https://zenodo.org/record/10820/files/semisynthetic_multilattice_data_2.tar.bz2


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

if test ! -f ${data_directory}/2/ag/trp_ag_*.cbf; then
  echo "Data not found in directory: ${data_directory} -"
  echo "please download from doi:10.5281/zenodo.10820 and pass on command-line"
  echo "or create softlink from data directory to ./data"
  exit
fi

dials.import ${data_directory}/2/ag/trp_ag_*.cbf

# find spots on the images - this will take a while! also accept small spots
# in this case as the data are good

dials.find_spots min_spot_size=3 datablock.expt nproc=$nproc

# index these found spots, searching for multiple lattices

dials.index datablock.expt strong.refl \
  max_lattices=2

# refine each indexing solution (separately) in all Bravais settings consistent
# with the indexed unit cell. In this example we would continue processing
# using the bravais_setting_5.expt, i.e. solution number 5

dials.refine_bravais_settings indexed.expt indexed.refl crystal_id=0

dials.refine_bravais_settings indexed.expt indexed.refl crystal_id=1

# now re-run the indexing, this time imposing the lattice constraints for the
# chosen Bravais setting, in this case number 5, i.e. oP, or point group P222

dials.index datablock.expt strong.refl \
  max_lattices=2 \
  space_group=P222

dials.refine indexed.expt indexed.refl \
  scan_varying=True \
  outlier.algorithm=tukey

# now run the integration - complex choices of algorithms are shown here in
# terms of the peak measurement (fitrs) and background determination
# methods. pass reference reflections from indexing in to determine the
# profile parameters...

dials.integrate outlier.algorithm=null profile.fitting=True \
  input.experiments=refined.expt \
  input.reflections=indexed.refl \
  nproc=$nproc

# currently dials.export only supports one experiment at a time
# therefore split the refined.expt and integrated.refl into
# separate files

dials.split_experiments refined.expt integrated.refl \
  experiments_prefix=refined reflections_prefix=integrated

# finally export the integrated measurements in an MTZ file - this should be
# properly formatted for immediate use in pointless / aimless

dials.export integrated_0.refl refined_0.expt mtz.hklout=integrated_0.mtz
dials.export integrated_1.refl refined_1.expt mtz.hklout=integrated_1.mtz

rebatch hklin integrated_0.mtz hklout rebatch_0.mtz > rebatch_0.log << EOF
batch add 0
EOF

rebatch hklin integrated_1.mtz hklout rebatch_1.mtz > rebatch_1.log << EOF
batch add 200
EOF

# and as if to make a point, here is all we need to do is to sort the data with
# pointless and then scale the data (ignoring anomalous differences) to 1.3A,
# and run ctruncate for the intensity analysis...

pointless hklin rebatch_0.mtz rebatch_1.mtz hklout sorted.mtz > pointless.log

aimless hklin sorted.mtz hklout scaled.mtz > aimless.log << eof
#resolution 1.3
anomalous off
eof

ctruncate -hklin scaled.mtz -hklout truncated.mtz \
-colin '/*/*/[IMEAN,SIGIMEAN]' > ctruncate.log

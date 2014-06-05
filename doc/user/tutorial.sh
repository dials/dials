############################################################################
# DIALS tutorial script - associated with data doi:10.5281/zenodo.10271
############################################################################

# first import the data - this reads the image headers and creates a file which
# describes the sweeps found therein. here there should only be one sweep found

dials.import ~/data/th_8_2_0*.cbf

# find spots on the images - this will take a while! also accept small spots
# in this case as the data are good (N.B. to explore DIALS options you can
# run dials.parameters - this shows you the entire parameter set which is
# available for DIALS)

dials.find_spots min_spot_size=3 datablock.json

# index these found spots - in this case the crystal is thaumatin which is known
# to be tetragonal, so impose appropriate lattice constraints in the indexing
# and refinement. can also specify unit cell here and also apply different
# indexing algorithms

dials.index space_group=P4 datablock.json strong.pickle

# run the refinement of the data - the indexing includes refinement so the result
# will be refined already with a static model - however here we would like to use
# a scan varying model so re-run the refinement (you should find that the R.M.S.
# deviations are a little lower following the scan varying refinement)

dials.refine experiments.json indexed.pickle scan_varying=true

# now run the integration - at the moment need to specify the number of blocks
# to process the data in. if this is not set you may find that the computer
# struggles... 540 images makes good sense as 6 blocks of 9 images / 13.5 degrees.
# complex choices of algorithms are shown here in terms of the peak measurement
# and background determination methods.

# OLD command line

dials.integrate n_blocks=6 outlier.algorithm=null \
  intensity.algorithm=fitrs refined_experiments.json indexed.pickle

# NEW command-line

dials.integrate block_size=13.5 outlier.algorithm=null \
  intensity.algorithm=fitrs refined_experiments.json -r indexed.pickle


# finally export the integrated measurements in an MTZ file - this should be
# properly formatted for immediate use in pointless / aimless

dials.export_mtz integrated.pickle refined_experiments.json integrated.mtz

# and as if to make a point, here is all we need to do is to sort the data with
# pointless and then scale the data (ignoring anomalous differences) to 1.3A

pointless hklin integrated.mtz hklout sorted.mtz

aimless hklin sorted.mtz hklout scaled.mtz << eof
resolution 1.3
anomalous off
eof

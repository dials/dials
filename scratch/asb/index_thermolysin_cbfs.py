from __future__ import division
from libtbx import easy_run
import sys, os, shutil

results = []

for path in sys.argv[1:]:
  basename  = os.path.splitext(os.path.basename(path))[0]
  datablock = basename + "_datablock.json"
  spots     = basename + "_strong.pickle"

  easy_run.call("dials.import %s --output %s"%(path,datablock))
  easy_run.call("dials.find_spots %s threshold.xds.sigma_strong=15 min_spot_size=2 -o %s"%(datablock, spots))
  easy_run.call("dials.index %s %s method=fft1d beam.fix=all detector.fix=all known_symmetry.unit_cell=93,93,130,90,90,120 known_symmetry.space_group=P6122"%(spots, datablock))

  if os.path.exists("indexed.pickle"):
    assert os.path.exists("experiments.json")
    print basename, "indexed succesfully"

    indexed = basename + "_indexed.pickle"
    experiments = basename + "_experiments.json"

    results.append((indexed,experiments))

    shutil.move("indexed.pickle", indexed)
    shutil.move("experiments.json", experiments)
  else:
    print basename, "failed to index"

print "Writing phil input to metrology refiner..."

f = open("all.phil", 'w')
for (indexed, experiments) in results:
  f.write("input {\n")
  f.write("  experiments = %s\n"%experiments)
  f.write("  reflections = %s\n"%indexed)
  f.write("}\n")
f.close()

print "Done"

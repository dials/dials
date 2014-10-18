from __future__ import division
import os, sys

paths = sys.argv[1:]

for path in paths:
  for filename in os.listdir(path):
    if "integrated" in filename:
      print "input {"
      print "  experiments =", os.path.join(path, filename.rstrip("_integrated.pickle") + "_refined_experiments.json")
      print "  reflections =", os.path.join(path, filename)
      print "}"


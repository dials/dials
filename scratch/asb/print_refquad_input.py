from __future__ import division
import os, sys

path = sys.argv[1]

for filename in os.listdir(path):
  if "integrated" in filename:
    print "input {"
    print "  experiments =", filename.rstrip("_integrated.pickle") + "_refined_experiments.json"
    print "  reflections =", filename
    print "}"


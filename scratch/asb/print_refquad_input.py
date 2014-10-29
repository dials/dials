from __future__ import division
import os, sys

paths = sys.argv[1:]

for path in paths:
  for filename in os.listdir(path):
    if "indexed" in filename:
      exp_path = os.path.join(path, filename.rstrip("_indexed.pickle") + "_refined_experiments.json")
      if not os.path.exists(exp_path): continue
      print "input {"
      print "  experiments =", exp_path
      print "  reflections =", os.path.join(path, filename)
      print "}"


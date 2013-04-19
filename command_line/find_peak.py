from __future__ import division
from dxtbx.sweep import SweepFactory
from dials.algorithms.peak_finding.spot_finder_lui import SpotFinderLui
#from dials.scratch.luiso_s.to_dials_reg.func_fnd_pk import *

import numpy
def fnd_pk():
    import sys
    from dxtbx.format.Registry import Registry
    print len(sys.argv[1:]), "images given"
    #tm = 2
    #print "3D peak find"


    filenames = sys.argv[1:]
    sweep = SweepFactory.sweep(filenames)
    find_spots = SpotFinderLui()
    spot = find_spots(sweep)

if __name__ == '__main__':
    fnd_pk()


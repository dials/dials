from __future__ import division

from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d, write_2d, test_compare_2d
from matplotlib import pyplot as plt

# from dials.scratch.luiso_s import add_2d

from dials.algorithms.integration import add_2d

nrow = 10
ncol = 10

sumation = flex.double(flex.grid(21, 21))
descr = flex.double(flex.grid(1, 3))
descr[0, 0] = .5
descr[0, 1] = 5
descr[0, 2] = 5

ref2d = model_2d(5, 5, 2, 1, 0.2, 55, 0.5)

sumation = add_2d(descr, ref2d, sumation)
#sumation = test_compare_2d(descr, ref2d, sumation)

write_2d(sumation)


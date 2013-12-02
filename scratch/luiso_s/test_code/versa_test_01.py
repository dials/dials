from __future__ import division

from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d, write_2d, test_compare_2d
from matplotlib import pyplot as plt

from dials.scratch.luiso_s import add_2d

#from dials.algorithms.integration import add_2d

nrow = 10
ncol = 10

sumation = flex.double(flex.grid(9, 9),1)
descr = flex.double(flex.grid(1, 3))
descr[0, 0] = 2.5
descr[0, 1] = 2.5
descr[0, 2] = 1

#double centr_col = descriptor(0,0);
#double centr_row = descriptor(0,1);
#double scale = descriptor(0,2);


ref2d = model_2d(5, 5, 2, 1, 0.2, 55, 0.5)
write_2d(ref2d)
write_2d(sumation)
print "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
sumation01 = add_2d(descr, ref2d, sumation)

write_2d(ref2d)
write_2d(sumation01)
write_2d(sumation)


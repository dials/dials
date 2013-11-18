from __future__ import division
from scitbx.array_family import flex
from dials.scratch.luiso_s import model_2d, measure_2d_angl
import numpy
from matplotlib import pyplot as plt
ncol = 36
nrow = 36

for int_ref_ang in range(20):

#int_ref_ang = 10
  print "_________________________________________"
  ref_ang = float(int_ref_ang) / 20.0
  print "ref_ang =", ref_ang
  ref2d = model_2d(nrow, ncol, 3, 7, ref_ang, 25, 0.5)

  dat2d_ref = ref2d.as_numpy_array()
  mask2d = numpy.copy(dat2d_ref)
  for col in range(ncol):
    for row in range(nrow):
      if mask2d[row, col] > 2:
        mask2d[row, col] = 1
      else:
        mask2d[row, col] = 0
  msk_felx = flex.int(mask2d)
  #for pr_line in dat2d_ref:
  #    print pr_line
  print "_________________________________________"
  for pr_line in mask2d:
    print pr_line

  ang = measure_2d_angl(ref2d, msk_felx, 18, 18)

  print "ang =", ang
  print "amg/pi =", ang / 3.14159265358


  new_ref2d = model_2d(30, 30, 2, 5, 3.14159265358 - ang / 3.14159265358, 25, 0.5)
  dat2d_paint = new_ref2d.as_numpy_array()
  plt.imshow(dat2d_paint , interpolation = "nearest")
  plt.show()

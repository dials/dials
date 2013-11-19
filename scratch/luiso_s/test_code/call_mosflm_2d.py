from __future__ import division
from dials.scratch.luiso_s.test_code.mosflm_2D import \
 fit_profile_2d, calc_background_n_make_2d_profile

from dials.model.data import Reflection, ReflectionList
def mosflm_caller(rlist, xmax, ymax, n_div):
  ncol = n_div
  nrow = n_div
  arr_rlist = []

  for col in range(ncol):
    b = []
    for row in range(nrow):
      b.append(ReflectionList())

    arr_rlist.append(b)

  ncnt = 0
  lst_pos = []
  for r in rlist:
    if r.is_valid():
      #print r.is_valid(), r.bounding_box, r.shoebox.all()
      x, y = r.image_coord_px
      col = int(float(x) / float(xmax) * n_div)
      row = int(float(y) / float(ymax) * n_div)
      arr_rlist[row][col].append(r)
      ncnt += 1
      pos = [row, col, len(arr_rlist[row][col]) - 1]
      lst_pos.append(pos)

  for col in range(ncol):
    for row in range(nrow):
      profile, tr_hold = calc_background_n_make_2d_profile(arr_rlist[row][col])

      from matplotlib import pyplot as plt
      data2d = profile.as_numpy_array()
      plt.imshow(data2d, interpolation = "nearest", cmap = plt.gray())
      plt.show()

      arr_rlist[row][col] = fit_profile_2d(arr_rlist[row][col], profile, tr_hold)

  new_rlist = ReflectionList()
  for numpos in lst_pos:
    row = numpos[0]
    col = numpos[1]
    deep = numpos[2]
    new_rlist.append(arr_rlist[row][col][deep])
  #print new_rlist[1]
  return new_rlist

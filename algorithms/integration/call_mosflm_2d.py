#
# mosflm_like.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
from __future__ import division
from dials.algorithms.integration.mosflm_2D_profile import \
 fit_profile_2d, make_2d_profile

from dials.model.data import Reflection, ReflectionList

from dials.array_family import flex

#def mosflm_caller(rlist, xmax, ymax, n_div):
def mosflm_caller(ref_table_in, xmax, ymax, n_div):
  ncol = n_div
  nrow = n_div
  arr_rlist = []
  arr_proff = []
  print "Performing profile fitting  ...."
  for col in range(ncol):
    tmp_empty_ref_data = []
    tmp_empty_prof = []
    for row in range(nrow):
      tmp_empty_ref_data.append([])
      tmp_empty_prof.append([])
    arr_rlist.append(tmp_empty_ref_data)
    arr_proff.append(tmp_empty_prof)

  col_xyzcal = ref_table_in['xyzcal.px']

  #for r in rlist:
  for t_row in range(len(ref_table_in)):
    #in the future consider searcing for is_valid logical
    #if r.is_valid():

      # consider replasing xyzcal with centroid pos
      x, y = col_xyzcal[t_row][0:2]# r.image_coord_px

      col = int(float(x) / float(xmax) * n_div)
      row = int(float(y) / float(ymax) * n_div)
      log_print = '''
      print "x,y =", x, y
      print "col, row =", col, row
      print "___________________________end table row #", t_row
      #'''
      arr_rlist[row][col].append([t_row])
  print "Building profiles  ...."

  for col in range(ncol):
    for row in range(nrow):
      #profile, tr_hold = make_2d_profile(arr_rlist[row][col])
      profile, tr_hold = make_2d_profile(arr_rlist[row][col], ref_table_in)
      arr_proff[row][col] = [profile, tr_hold]

  print "Building profiles          ....       Done"
  for col in range(ncol):
    for row in range(nrow):
      ref_table_in = fit_profile_2d(arr_rlist[row][col], ref_table_in
                                    , arr_proff, row, col,  xmax, ymax)
      #arr_rlist[row][col] = fit_profile_2d(arr_rlist[row][col], ref_table_in
      #                      arr_proff, row, col,  xmax, ymax)
  print "profile fitting            ....       Done"

  new_ref_table = flex.reflection_table()
  return new_ref_table

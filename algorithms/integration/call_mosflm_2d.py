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

  for col in range(ncol):
    tmp_empty_ref_data = []
    tmp_empty_prof = []
    for row in range(nrow):
      tmp_empty_ref_data.append([])
      tmp_empty_prof.append([])
    arr_rlist.append(tmp_empty_ref_data)
    arr_proff.append(tmp_empty_prof)

  print "empty 3D list done"
  print "arr_proff =", arr_proff
  print "arr_rlist =", arr_rlist


  col_xyzcal = ref_table_in['xyzcal.px']

  #ncnt = 0
  #lst_pos = []

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
      #ncnt += 1
      #pos = [row, col, len(arr_rlist[row][col]) - 1]
      #lst_pos.append(pos)

  for col in range(ncol):
    for row in range(nrow):
      #profile, tr_hold = make_2d_profile(arr_rlist[row][col])
      profile, tr_hold = make_2d_profile(arr_rlist[row][col], ref_table_in)
      arr_proff[row][col] = [profile, tr_hold]

  print "performing profile fitting  ...."
  for col in range(ncol):
    for row in range(nrow):
      ref_table_in = fit_profile_2d(arr_rlist[row][col], ref_table_in
                                    , arr_proff, row, col,  xmax, ymax)
      #arr_rlist[row][col] = fit_profile_2d(arr_rlist[row][col], ref_table_in
      #                      arr_proff, row, col,  xmax, ymax)

  to_be_fixed_later = '''
  new_rlist = ReflectionList()
  for numpos in lst_pos:
    row = numpos[0]
    col = numpos[1]
    deep = numpos[2]
    new_rlist.append(arr_rlist[row][col][deep])
  print "profile fitting  .... done"
  #'''
  new_ref_table = flex.reflection_table()
  return new_ref_table

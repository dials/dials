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

  col_xyzcal = ref_table_in['xyzcal.px']

  for t_row in range(len(ref_table_in)):
    #in the future consider searcing for is_valid logical
    #if r.is_valid():

      # consider replasing xyzcal with centroid pos
      x, y = col_xyzcal[t_row][0:2]# r.image_coord_px

      col = int(float(x) / float(xmax) * ncol)
      row = int(float(y) / float(ymax) * nrow)
      log_print_n_debugg_way = '''
      print "x,y =", x, y
      print "col, row =", col, row
      print "___________________________end table row #", t_row

      try:
        arr_rlist[row][col].append([t_row])
      except Exception:
        from dials.util.command_line import interactive_console; interactive_console()
        break
      #'''
      arr_rlist[row][col].append([t_row])

  tbl_siz = ncol * nrow


  from dials.util.command_line import ProgressBar
  p_bar = ProgressBar(title = 'Performing profile building')
  tbl_prgr = 0
  for col in range(ncol):
    for row in range(nrow):

      p_bar.update(tbl_prgr * 100.0 / tbl_siz)
      tbl_prgr += 1

      profile, tr_hold = make_2d_profile(arr_rlist[row][col], ref_table_in)
      arr_proff[row][col] = [profile, tr_hold]
  p_bar.finished('Done Building profiles')

  p_bar = ProgressBar(title = 'Performing profile fitting')
  tbl_prgr = 0

  for col in range(ncol):
    for row in range(nrow):

      p_bar.update(tbl_prgr * 100.0 / tbl_siz)
      tbl_prgr += 1

      ref_table_in = fit_profile_2d(arr_rlist[row][col], ref_table_in
                                    , arr_proff, row, col,  xmax, ymax)
  p_bar.finished('Done profiles fitting')


  #new_ref_table = flex.reflection_table()
  #return new_ref_table
  return ref_table_in

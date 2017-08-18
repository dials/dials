from __future__ import absolute_import, division
import libtbx.phil
# LIBTBX_SET_DISPATCHER_NAME dev.dials.rfactor

scope = libtbx.phil.parse("""
  hklin = None
    .type = path
  d_min = 0
    .type = float
  fo_col = 'F'
    .type = str
  fc_col = 'FC_ALL'
    .type = str
  free_col = 'FreeR_flag'
    .type = str
""")

def find_column(miller_array_dict, name):
  for k in miller_array_dict:
    if k[-1] == name:
      return miller_array_dict[k]
  raise RuntimeError, 'Could not find column %s' % name

def main(args):
  work = scope
  for arg in args:
    work = work.fetch(libtbx.phil.parse(arg))
  params = work.extract()

  from iotbx import mtz

  m = mtz.object(params.hklin)
  mad = m.as_miller_arrays_dict()

  fo = find_column(mad, params.fo_col)
  fc = find_column(mad, params.fc_col)
  free = find_column(mad, params.free_col)

  data_d_min = fc.d_min()

  # filter on resolution - only need to do one column as match later
  if params.d_min:
    fc = fc.resolution_filter(d_min=params.d_min)

  # find common reflections i.e. in Fo and Fc set (can assume all Free flags)
  fo, fc = fo.common_sets(fc)
  ff, fc = free.common_sets(fc)

  # split off work and free set
  free_set = (ff.data() == 0)
  fow = fo.select(~free_set)
  fcw = fc.select(~free_set)
  fof = fo.select(free_set)
  fcf = fc.select(free_set)
  R = fow.r1_factor(fcw)
  Rf = fof.r1_factor(fcf)

  nw = free_set.count(False)
  nf = free_set.count(True)

  print 'resol  Rwork Rfree #work #free'
  print '%.3f %.4f %.4f %d %d' % (data_d_min, R, Rf, nw, nf)

if __name__ == '__main__':
  import sys
  main(sys.argv[1:])

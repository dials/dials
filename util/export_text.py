from __future__ import absolute_import, division
def export_text(integrated_data):
  '''Export contents of a dials reflection table as text.'''

  mi = integrated_data['miller_index']
  h, k, l = zip(*integrated_data['miller_index'])

  # FIXME Currently outputting either summation or profile fitting. Should do
  # both?
  if 'intensity.prf' in integrated_data:
    i = integrated_data['intensity.prf.value']
    v = integrated_data['intensity.prf.variance']
  else:
    i = integrated_data['intensity.sum.value']
    v = integrated_data['intensity.sum.variance']
  lp = integrated_data['lp']
  i *= lp
  v *= lp

  for _h, _k, _l, _i, _v in zip(h, k, l, i, v):
    print '%4d %4d %4d %f %f' % (_h, _k, _l, _i, _v)

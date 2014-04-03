def export_text(integrated_data):
  '''Export contents of a dials reflection table as text.'''

  mi = integrated_data['miller_index']
  h, k, l = zip(*integrated_data['miller_index'])

  i = integrated_data['intensity.cor.value']
  v = integrated_data['intensity.cor.variance']

  for _h, _k, _l, _i, _v in zip(h, k, l, i, v):
    print '%4d %4d %4d %f %f' % (_h, _k, _l, _i, _v)

  return

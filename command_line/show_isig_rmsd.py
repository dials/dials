from __future__ import division

def show_isig_rmsd(integrated_data):
  # check I have the columns I want

  assert('miller_index' in integrated_data)
  assert('xyzcal.px' in integrated_data)
  assert('xyzobs.px.value' in integrated_data)
  assert(('intensity.prf.value' in integrated_data) or \
         ('intensity.sum.value' in integrated_data))

  # strip out bad measurements

  if 'intensity.prf.value' in integrated_data:
    selection = integrated_data['intensity.prf.variance'] <= 0.0
    if selection.count(True) > 0:
      integrated_data.del_selected(selection)
  else:
    selection = integrated_data['intensity.sum.variance'] <= 0.0
    if selection.count(True) > 0:
      integrated_data.del_selected(selection)

  x_px, y_px, z_px = integrated_data['xyzcal.px'].parts()
  x_opx, y_opx, z_opx = integrated_data['xyzobs.px.value'].parts()

  if 'intensity.prf.value' in integrated_data:
    I = integrated_data['intensity.prf.value']
    V = integrated_data['intensity.prf.variance']
  else:
    I = integrated_data['intensity.sum.value']
    V = integrated_data['intensity.sum.variance']

  batch = flex.floor(z_px).iround()

  for b in range(min(batch), max(batch) + 1):
    sel = batch == b
    i = I.select(sel)
    v = V.select(sel)
    x_px_sel = x_px.select(sel)
    y_px_sel = y_px.select(sel)
    x_opx_sel = x_opx.select(sel)
    y_opx_sel = y_opx.select(sel)
    i_s = i / flex.sqrt(v)
    rmsd = flex.sqrt(flex.pow2(x_px_sel - x_opx_sel) +
                     flex.pow2(y_px_sel - y_opx_sel))
    n = len(i)
    print '%d %d %.3f %3f' % (b, n, flex.sum(i_s) / n, flex.sum(rmsd) / n)

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 2:
    raise RuntimeError, '%s integrated.pickle'

  import cPickle as pickle
  from dials.array_family import flex

  integrated_data = pickle.load(open(sys.argv[1], 'rb'))
  show_isig_rmsd(integrated_data)

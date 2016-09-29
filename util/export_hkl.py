from __future__ import division

from export_mtz import sum_partial_reflections
from export_mtz import scale_partial_reflections

def export_hkl(integrated_data, experiment_list, hklout, run=0,
               summation=False, include_partials=False, keep_partials=False):
  '''Export data from integrated_data corresponding to experiment_list to a
  HKL file for input to SADABS. FIXME probably need to make a .p4p file as
  well...'''

  from logging import info
  from dials.array_family import flex

  # for the moment assume (and assert) that we will convert data from exactly
  # one lattice...

  assert(len(experiment_list) == 1)
  # select reflections that are assigned to an experiment (i.e. non-negative id)

  integrated_data = integrated_data.select(integrated_data['id'] >= 0)
  assert max(integrated_data['id']) == 0

  if not summation:
    assert('intensity.prf.value' in integrated_data)

  # strip out negative variance reflections: these should not really be there
  # FIXME Doing select on summation results. Should do on profile result if
  # present? Yes

  if 'intensity.prf.variance' in integrated_data:
    selection = integrated_data.get_flags(
      integrated_data.flags.integrated,
      all=True)
  else:
    selection = integrated_data.get_flags(
      integrated_data.flags.integrated_sum)
  integrated_data = integrated_data.select(selection)

  selection = integrated_data['intensity.sum.variance'] <= 0
  if selection.count(True) > 0:
    integrated_data.del_selected(selection)
    info('Removing %d reflections with negative variance' % \
          selection.count(True))

  if 'intensity.prf.variance' in integrated_data:
    selection = integrated_data['intensity.prf.variance'] <= 0
    if selection.count(True) > 0:
      integrated_data.del_selected(selection)
      info('Removing %d profile reflections with negative variance' % \
            selection.count(True))

  # FIXME in here work on including partial reflections => at this stage best
  # to split off the partial refections into a different selection & handle
  # gracefully... better to work on a short list as will need to "pop" them &
  # find matching parts to combine.

  if include_partials:
    integrated_data = sum_partial_reflections(integrated_data)
    integrated_data = scale_partial_reflections(integrated_data)

  if 'partiality' in integrated_data:
    selection = integrated_data['partiality'] < 0.99
    if selection.count(True) > 0 and not keep_partials:
      integrated_data.del_selected(selection)
      info('Removing %d incomplete reflections' % \
        selection.count(True))

  experiment = experiment_list[0]

  # map to ASU first
  from cctbx.miller import map_to_asu
  map_to_asu(experiment.crystal.get_space_group().type(), False,
      integrated_data['miller_index'])

  # sort data before output
  nref = len(integrated_data['miller_index'])
  indices = flex.size_t_range(nref)
  perm = sorted(indices, key=lambda k: integrated_data['miller_index'][k])
  integrated_data = integrated_data.select(flex.size_t(perm))

  # tricky bit here will be calculating direction cosines of the incoming
  # beam and the outgoing beam w.r.t. the crystal unit cell axes at the time
  # so need to compute R U B for each reflection then - cache these in RUBs

  from scitbx import matrix

  assert (not experiment.goniometer is None)

  setting_rotation = matrix.sqr(experiment.goniometer.get_setting_rotation())
  axis = setting_rotation * matrix.col(
    experiment.goniometer.get_rotation_axis())

  beam = matrix.col(experiment.beam.get_direction())
  s0 = matrix.col(experiment.beam.get_s0())

  F = matrix.sqr(experiment.goniometer.get_fixed_rotation())
  S = matrix.sqr(experiment.goniometer.get_setting_rotation())
  unit_cell = experiment.crystal.get_unit_cell()

  from iotbx import mtz

  from scitbx.array_family import flex
  from math import floor, sqrt

  assert(not experiment.scan is None)
  image_range = experiment.scan.get_image_range()

  # somewhere to stash R*U*B for each frame where R is full goniometer
  # matrix i.e. setting * scan * fixed matrix

  RUBs = { }

  for b in range(image_range[0], image_range[1] + 1):
    UB = experiment.crystal.get_A_at_scan_point(b - image_range[0])
    phi_start, phi_range = experiment.scan.get_image_oscillation(b)
    phi = phi_start + 0.5 * phi_range
    R = axis.axis_and_angle_as_r3_rotation_matrix(phi, deg=True)
    RUBs[b - image_range[0]] = S * R * F * UB

  from cctbx.array_family import flex as cflex # implicit import
  from cctbx.miller import map_to_asu_isym # implicit import

  # gather the required information for the reflection file

  nref = len(integrated_data['miller_index'])
  zdet = flex.double(integrated_data['xyzcal.px'].parts()[2])

  # only using this to index into RUBs above
  iframe = flex.floor(zdet).iround()

  miller_index = integrated_data['miller_index']

  I = None
  sigI = None

  # export including scale factors

  if 'lp' in integrated_data:
    lp = integrated_data['lp']
  else:
    lp = flex.double(nref, 1.0)
  if 'dqe' in integrated_data:
    dqe = integrated_data['dqe']
  else:
    dqe = flex.double(nref, 1.0)
  scl = lp / dqe

  if summation:
    I = integrated_data['intensity.sum.value'] * scl
    V = integrated_data['intensity.sum.variance'] * scl * scl
    assert V.all_gt(0)
    sigI = flex.sqrt(V)
  else:
    I = integrated_data['intensity.prf.value'] * scl
    V = integrated_data['intensity.prf.variance'] * scl * scl
    assert V.all_gt(0)
    sigI = flex.sqrt(V)

  # figure out scaling to make sure data fit into format
  # 3I4, 2F8.2, I4, 6F8.5

  Imax = flex.max(I)

  info('Maximum intensity in file: %8.2f' % Imax)

  if Imax > 99999.0:
    scale = 99999.0 / Imax
    I = I * scale
    sigI = sigI * scale

  # detector scaling info
  assert(len(experiment.detector) == 1)
  panel = experiment.detector[0]
  scl_x = 512.0 / panel.get_image_size()[0]
  scl_y = 512.0 / panel.get_image_size()[1]

  fout = open(hklout, 'w')
  for j in range(nref):
    h, k, l = miller_index[j]
    x_px, y_px, z_px = integrated_data['xyzcal.px'][j]
    lp = scl[j]
    istol = int(round(10000 * unit_cell.stol((h, k, l))))
    RUB = RUBs[iframe[j]]
    x = RUB * (h, k, l)
    s = (s0 + x).normalize()
    astar = (RUB * (1, 0, 0)).normalize()
    bstar = (RUB * (0, 1, 0)).normalize()
    cstar = (RUB * (0, 0, 1)).normalize()
    ix = beam.dot(astar)
    iy = beam.dot(bstar)
    iz = beam.dot(cstar)
    dx = s.dot(astar)
    dy = s.dot(bstar)
    dz = s.dot(cstar)
    x = x_px * scl_x
    y = y_px * scl_y
    fout.write('%4d%4d%4d%8.2f%8.2f%4d%8.5f%8.5f%8.5f%8.5f%8.5f%8.5f' % \
               (h, k, l, I[j], sigI[j], run, ix, dx, iy, dy, iz, dz))
    fout.write('%7.2f%7.2f%8.2f%7.3f%5d\n' % (x, y, z_px, lp, istol))
  fout.close()
  info('Output %d reflections to %s' % (nref, hklout))
  return

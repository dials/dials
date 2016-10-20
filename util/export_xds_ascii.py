from __future__ import division

import logging
logger = logging.getLogger(__name__)

from export_mtz import sum_partial_reflections
from export_mtz import scale_partial_reflections

def export_xds_ascii(integrated_data, experiment_list, hklout, summation=False,
                     include_partials=False, keep_partials=False, var_model=(1,0)):
  '''Export data from integrated_data corresponding to experiment_list to
  an XDS_ASCII.HKL formatted text file.'''

  from dials.array_family import flex
  import math

  # for the moment assume (and assert) that we will convert data from exactly
  # one lattice...

  assert(len(experiment_list) == 1)
  # select reflections that are assigned to an experiment (i.e. non-negative id)

  integrated_data = integrated_data.select(integrated_data['id'] >= 0)
  assert max(integrated_data['id']) == 0

  if not summation:
    assert('intensity.prf.value' in integrated_data)

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
    logger.info('Removing %d reflections with negative variance' % \
          selection.count(True))

  if 'intensity.prf.variance' in integrated_data:
    selection = integrated_data['intensity.prf.variance'] <= 0
    if selection.count(True) > 0:
      integrated_data.del_selected(selection)
      logger.info('Removing %d profile reflections with negative variance' % \
            selection.count(True))

  if include_partials:
    integrated_data = sum_partial_reflections(integrated_data)
    integrated_data = scale_partial_reflections(integrated_data)

  if 'partiality' in integrated_data:
    selection = integrated_data['partiality'] < 0.99
    if selection.count(True) > 0 and not keep_partials:
      integrated_data.del_selected(selection)
      logger.info('Removing %d incomplete reflections' % \
        selection.count(True))

  experiment = experiment_list[0]

  # sort data before output
  nref = len(integrated_data['miller_index'])
  indices = flex.size_t_range(nref)

  import copy
  unique = copy.deepcopy(integrated_data['miller_index'])
  from cctbx.miller import map_to_asu
  map_to_asu(experiment.crystal.get_space_group().type(), False, unique)

  perm = sorted(indices, key=lambda k: unique[k])
  integrated_data = integrated_data.select(flex.size_t(perm))

  from scitbx import matrix
  from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame

  assert (not experiment.goniometer is None)

  unit_cell = experiment.crystal.get_unit_cell()

  from scitbx.array_family import flex
  from math import floor, sqrt

  assert(not experiment.scan is None)
  image_range = experiment.scan.get_image_range()
  phi_start, phi_range = experiment.scan.get_image_oscillation(image_range[0])

  # gather the required information for the reflection file

  nref = len(integrated_data['miller_index'])
  zdet = flex.double(integrated_data['xyzcal.px'].parts()[2])

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
  Iscl = lp / dqe
  Vscl = lp * lp / dqe

  # profile correlation
  if 'profile.correlation' in integrated_data:
    prof_corr = 100.0 * integrated_data['profile.correlation']
  else:
    prof_corr = flex.double(nref, 100.0)

  # partiality
  if 'partiality' in integrated_data:
    partiality = 100 * integrated_data['partiality']
  else:
    prof_corr = flex.double(nref, 100.0)

  if summation:
    I = integrated_data['intensity.sum.value'] * Iscl
    V = integrated_data['intensity.sum.variance'] * Vscl
    assert V.all_gt(0)
    V = var_model[0] * (V + var_model[1] * I * I)
    sigI = flex.sqrt(V)
  else:
    I = integrated_data['intensity.prf.value'] * Iscl
    V = integrated_data['intensity.prf.variance'] * Vscl
    assert V.all_gt(0)
    V = var_model[0] * (V + var_model[1] * I * I)
    sigI = flex.sqrt(V)

  fout = open(hklout, 'w')

  # first write the header - in the "standard" coordinate frame...

  panel = experiment.detector[0]
  fast = panel.get_fast_axis()
  slow = panel.get_slow_axis()
  Rd = align_reference_frame(fast, (1,0,0), slow, (0,1,0))
  print 'Coordinate change:'
  print '%5.2f %5.2f %5.2f\n%5.2f %5.2f %5.2f\n%5.2f %5.2f %5.2f\n' % Rd.elems

  fast = Rd * fast
  slow = Rd * slow

  qx, qy = panel.get_pixel_size()
  nx, ny = panel.get_image_size()
  distance = matrix.col(Rd * panel.get_origin()).dot(
      matrix.col(Rd * panel.get_normal()))
  org = Rd * (matrix.col(panel.get_origin()) - distance * matrix.col(
      panel.get_normal()))
  orgx = - org.dot(fast) / qx
  orgy = - org.dot(slow) / qy

  UB = Rd * matrix.sqr(experiment.crystal.get_A())
  real_space_ABC = UB.inverse().elems

  axis = Rd * experiment.goniometer.get_rotation_axis()
  beam = Rd * experiment.beam.get_s0()
  cell_fmt = '%9.3f %9.3f %9.3f %7.3f %7.3f %7.3f'
  axis_fmt = '%9.3f %9.3f %9.3f'

  fout.write('\n'.join([
    '!FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL\'S_LAW=TRUE',
    '!Generated by dials.export',
    '!DATA_RANGE= %d %d' % image_range,
    '!ROTATION_AXIS= %9.6f %9.6f %9.6f' % axis.elems,
    '!OSCILLATION_RANGE= %f' % phi_range,
    '!STARTING_ANGLE= %f' % phi_start,
    '!STARTING_FRAME= %d' % image_range[0],
    '!SPACE_GROUP_NUMBER= %d' % experiment.crystal.get_space_group().type().number(),
    '!UNIT_CELL_CONSTANTS= %s' % (cell_fmt % unit_cell.parameters()),
    '!UNIT_CELL_A-AXIS= %s' % (axis_fmt % real_space_ABC[0:3]),
    '!UNIT_CELL_B-AXIS= %s' % (axis_fmt % real_space_ABC[3:6]),
    '!UNIT_CELL_C-AXIS= %s' % (axis_fmt % real_space_ABC[6:9]),
    '!X-RAY_WAVELENGTH= %f' % experiment.beam.get_wavelength(),
    '!INCIDENT_BEAM_DIRECTION= %f %f %f' % beam.elems,
    '!NX= %d NY= %d QX= %f QY= %f' % (nx, ny, qx, qy),
    '!ORGX= %9.2f ORGY= %9.2f' % (orgx, orgy),
    '!DETECTOR_DISTANCE= %8.3f' % distance,
    '!DIRECTION_OF_DETECTOR_X-AXIS= %9.5f %9.5f %9.5f' % fast.elems,
    '!DIRECTION_OF_DETECTOR_Y-AXIS= %9.5f %9.5f %9.5f' % slow.elems,
    '!VARIANCE_MODEL= %7.3e %7.3e' % var_model,
    '!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=12',
    '!ITEM_H=1',
    '!ITEM_K=2',
    '!ITEM_L=3',
    '!ITEM_IOBS=4',
    '!ITEM_SIGMA(IOBS)=5',
    '!ITEM_XD=6',
    '!ITEM_YD=7',
    '!ITEM_ZD=8',
    '!ITEM_RLP=9',
    '!ITEM_PEAK=10',
    '!ITEM_CORR=11',
    '!ITEM_PSI=12',
    '!END_OF_HEADER',
    '']))

  # then write the data records

  s0 = Rd * matrix.col(experiment.beam.get_s0())

  for j in range(nref):
    x, y, z = integrated_data['xyzcal.px'][j]
    phi = phi_start + z * phi_range
    h, k, l = miller_index[j]
    X = (UB * (h, k, l)).rotate(axis, phi, deg=True)
    s = s0 + X
    g = s.cross(s0).normalize()
    f = (s - s0).normalize()

    # find component of beam perpendicular to f, e
    e = - (s + s0).normalize()
    if h == k and k == l:
      u = (h, -h, 0)
    else:
      u = (k - l, l - h, h - k)
    q = (matrix.col(u).transpose() * UB.inverse()).normalize(
        ).transpose().rotate(axis, phi, deg=True)

    psi = q.angle(g, deg=True)
    if q.dot(e) < 0:
      psi *= -1

    fout.write('%d %d %d %f %f %f %f %f %f %.1f %.1f %f\n' %
               (h, k, l, I[j], sigI[j], x, y, z, lp[j],
                partiality[j], prof_corr[j], psi))

  fout.write('!END_OF_DATA\n')
  fout.close()
  logger.info('Output %d reflections to %s' % (nref, hklout))
  return

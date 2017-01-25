from __future__ import absolute_import, division

# LIBTBX_SET_DISPATCHER_NAME dev.dials.compare_xds_dials

def pull_reference(integrate_hkl, d_min = 0.0):
  '''Generate reference data set from integrate.hkl, check out the calculated
  x, y and z centroids as well as the Miller indices as coordinates in some
  high dimensional space. Only consider measurements with meaningful
  centroids...'''

  uc = integrate_hkl_to_unit_cell(integrate_hkl)

  hkl = []
  i = []
  sigi = []
  xyz = []
  lp = []

  for record in open(integrate_hkl):
    if record.startswith('!'):
      continue

    f_tokens = map(float, record.split())

    if f_tokens[12:15] == [0.0, 0.0, 0.0]:
      continue

    _hkl = tuple(map(int, f_tokens[0:3]))
    if uc.d(_hkl) < d_min:
      continue

    hkl.append(_hkl)

    # need to scale by PEAK stored value to get comparison value

    peak = 0.01 * f_tokens[9]
    i.append(f_tokens[3] * peak)
    sigi.append(f_tokens[4] * peak)
    xyz.append(tuple(f_tokens[5:8]))
    lp.append(f_tokens[8])

  print 'Reference: %d observations' % len(hkl)
  return hkl, i, sigi, xyz, lp

def get_dials_matrix(experiments_json):
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  experiments = ExperimentListFactory.from_json_file(experiments_json)
  return experiments[0].crystal.get_A()

def get_dials_coordinate_frame(experiments_json):
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  experiments = ExperimentListFactory.from_json_file(experiments_json)
  return experiments[0].beam.get_direction(), \
    experiments[0].goniometer.get_rotation_axis()

def get_xds_coordinate_frame(integrate_hkl):
  from scitbx import matrix
  beam = axis = None

  for record in open(integrate_hkl):
    if not record.startswith('!'):
      break
    if record.startswith('!INCIDENT_BEAM_DIRECTION='):
      beam = matrix.col(map(float, record.split()[-3:])).normalize()
    if record.startswith('!ROTATION_AXIS='):
      axis = matrix.col(map(float, record.split()[-3:])).normalize()

  if not beam or not axis:
    raise RuntimeError, 'coordinate frame information not found'

  return beam, axis

def integrate_hkl_to_A_matrix(integrate_hkl):
  a = b = c = None

  for record in open(integrate_hkl):
    if not record.startswith('!'):
      break
    if record.startswith('!UNIT_CELL_A-AXIS='):
      a = tuple(map(float, record.split()[-3:]))
    if record.startswith('!UNIT_CELL_B-AXIS='):
      b = tuple(map(float, record.split()[-3:]))
    if record.startswith('!UNIT_CELL_C-AXIS='):
      c = tuple(map(float, record.split()[-3:]))

  if not a or not b or not c:
    raise RuntimeError, 'unit cell vectors not found'

  from scitbx import matrix
  return matrix.sqr(a + b + c).inverse()

def integrate_hkl_to_unit_cell(integrate_hkl):
  '''Generate a cctbx unit_cell from an integrate_hkl file.'''

  from cctbx.uctbx import unit_cell

  for record in open(integrate_hkl):
    if not record.startswith('!'):
      break
    if record.startswith('!UNIT_CELL_CONSTANTS='):
      return unit_cell(tuple(map(float, record.split()[-6:])))

  raise RuntimeError, 'unit cell not found'

def pull_calculated(integrate_pkl):
  from dials.array_family import flex # import dependency
  import cPickle as pickle
  import math

  r_list = pickle.load(open(integrate_pkl, 'rb'))

  strong_reflections = []

  for r in r_list:
    if r['intensity.sum.value'] ** 2 < r['intensity.sum.variance']:
      continue
    if r['intensity.sum.value'] <= 0.0:
      continue
    strong_reflections.append(r)

  del(r_list)

  hkl = []
  i = []
  sigi = []
  xyz = []
  lp = []

  for r in strong_reflections:
    hkl.append(r['miller_index'])
    i.append(r['intensity.cor.value'])
    sigi.append(math.sqrt(r['intensity.cor.variance']))
    lp.append(r['intensity.cor.value'] / r['intensity.sum.value'])
    x, y, z = r['xyzcal.px']
    xyz.append((x, y, z))

  print 'Computed: %d observations' % len(hkl)
  return hkl, i, sigi, xyz, lp

def meansd(values):
  import math

  assert(len(values) > 3)

  mean = sum(values) / len(values)
  var = sum([(v - mean) * (v - mean) for v in values]) / (len(values) - 1)

  return mean, math.sqrt(var)

def cc(a, b):

  assert(len(a) == len(b))

  ma, sa = meansd(a)
  mb, sb = meansd(b)

  r = (1 / (len(a) - 1)) * sum([((a[j] - ma) / sa) * ((b[j] - mb) / sb)
                                for j in range(len(a))])

  return r

def R(calc, obs, scale = None):

  import math

  assert(len(calc) == len(obs))

  if not scale:
    scale = sum(obs) / sum(calc)

  return sum([math.fabs(math.fabs(o) - math.fabs(scale * c)) \
              for c, o in zip(calc, obs)]) / \
              sum([math.fabs(o) for o in obs]), scale

def meansd(values):
  import math
  mean = sum(values) / len(values)
  var = sum([(v - mean) ** 2 for v in values]) / len(values)
  return mean, math.sqrt(var)

def compare_chunks(integrate_hkl, integrate_pkl, experiments_json, d_min = 0.0):

  from cctbx.array_family import flex
  from annlib_ext import AnnAdaptor as ann_adaptor

  rdx = derive_reindex_matrix(experiments_json, integrate_hkl)

  print 'Reindex matrix:\n%d %d %d\n%d %d %d\n%d %d %d' % (rdx.elems)

  uc = integrate_hkl_to_unit_cell(integrate_hkl)

  xhkl, xi, xsigi, xxyz, xlp = pull_reference(integrate_hkl, d_min=d_min)
  dhkl, di, dsigi, dxyz, dlp = pull_calculated(integrate_pkl)

  reference = flex.double()
  query = flex.double()

  for xyz in xxyz:
    reference.append(xyz[0])
    reference.append(xyz[1])
    reference.append(xyz[2])

  for xyz in dxyz:
    query.append(xyz[0])
    query.append(xyz[1])
    query.append(xyz[2])

  # perform the match
  ann = ann_adaptor(data = reference, dim = 3, k = 1)
  ann.query(query)

  XDS = []
  DIALS = []
  HKL = []

  # perform the analysis
  for j, hkl in enumerate(dhkl):
    c = ann.nn[j]
    if hkl == tuple(rdx * xhkl[c]):
      XDS.append(xi[c])
      DIALS.append(di[j])
      HKL.append(hkl)

  # now compute resolution for every reflection - or at least each unique
  # Miller index...

  unique = set(HKL)

  resolutions = { }

  for hkl in unique:
    resolutions[hkl] = uc.d(hkl)

  # then resort the list in terms of resolution, then reverse it

  sort_me = []
  for hkl, xds, dials in zip(HKL, XDS, DIALS):
    sort_me.append((resolutions[hkl], xds, dials))

  sort_me.sort()
  sort_me.reverse()

  resolutions = [sm[0] for sm in sort_me]
  XDS = [sm[1] for sm in sort_me]
  DIALS = [sm[2] for sm in sort_me]

  # then extract the original observation structure

  print 'Paired %d observations' % len(XDS)

  scale = sum(XDS) / sum(DIALS)

  chunks = [(i, i + 1000) for i in range(0, len(XDS), 1000)]

  ccs = []
  rs = []
  ss = []

  for chunk in chunks:
    xds = XDS[chunk[0]:chunk[1]]
    dials = DIALS[chunk[0]:chunk[1]]
    resols = resolutions[chunk[0]:chunk[1]]

    if len(xds) < 100:
      break

    c = cc(dials, xds)
    r, s = R(dials, xds)
    print '%7d %4d %.3f %.3f %.3f %.3f %.3f' % \
      (chunk[0], len(xds), min(resols), max(resols), c, r, s)
    ccs.append(c)
    rs.append(r)
    ss.append(s)

  chunks = [j for j in range(len(chunks))]

  # kludge - if we fall off

  chunks = chunks[:len(rs)]

  from matplotlib import pyplot
  pyplot.xlabel('Chunk')
  pyplot.ylabel('Statistic')
  pyplot.title('Statistics for 1000 reflection-pair chunks')
  pyplot.plot(chunks, ccs, label = 'CC')
  pyplot.plot(chunks, rs, label = 'R')
  pyplot.plot(chunks, ss, label = 'K')
  pyplot.legend()
  pyplot.savefig('plot-vs-xds.png')
  pyplot.close()

  return

def derive_reindex_matrix(experiments_json, integrate_hkl):
  '''Derive a reindexing matrix to go from the orientation matrix used
  for XDS integration to the one used for DIALS integration.'''

  dA = get_dials_matrix(experiments_json)
  dbeam, daxis = get_dials_coordinate_frame(experiments_json)
  xbeam, xaxis = get_xds_coordinate_frame(integrate_hkl)

  # want to align XDS -s0 vector...
  from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
  R = align_reference_frame(- xbeam, dbeam, xaxis, daxis)
  xA = R * integrate_hkl_to_A_matrix(integrate_hkl)

  # assert that this should just be a simple integer rotation matrix
  # i.e. reassignment of a, b, c so...

  from scitbx import matrix
  return matrix.sqr(map(int, map(round, (dA.inverse() * xA).elems)))

if __name__ == '__main__':
  import sys
  if len(sys.argv) < 4:
    raise RuntimeError, \
      '%s INTEGRATE.HKL integrate.pickle experiments.json [dmin]' % \
      sys.argv[0]

  if len(sys.argv) == 4:
    compare_chunks(sys.argv[1], sys.argv[2], sys.argv[3])
  else:
    compare_chunks(sys.argv[1], sys.argv[2], sys.argv[3],
                   d_min = float(sys.argv[4]))

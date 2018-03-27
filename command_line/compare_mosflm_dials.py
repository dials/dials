from __future__ import absolute_import, division
from __future__ import print_function

# LIBTBX_SET_DISPATCHER_NAME dev.dials.compare_mosflm_dials

def pull_reference(integrate_mtz):
  '''Generate reference data set from integrate.hkl, check out the calculated
  x, y and z centroids as well as the Miller indices as coordinates in some
  high dimensional space. Only consider measurements with meaningful
  centroids...'''

  # prepare input data as
  # sortmtz hklin thau_2_001.mtz hklout sorted.mtz << eof
  # H K L M/ISYM BATCH
  # eof
  #
  # scala hklin sorted.mtz hklout summed.mtz << eof
  # run 1 all
  # scales constant
  # output unmerged
  # sdcorrection noadjust norefine both 1 0 0
  # eof

  from iotbx import mtz

  m = mtz.object(integrate_mtz)

  hkl = m.extract_original_index_miller_indices()

  b0 = m.batches()[0]

  for c in m.columns():
    if c.label() == 'I':
      i_lp = c.extract_valid_values()
    elif c.label() == 'SIGI':
      sigi = c.extract_valid_values()
    elif c.label() == 'XDET':
      xdet = c.extract_valid_values()
    elif c.label() == 'YDET':
      ydet = c.extract_valid_values()
    elif c.label() == 'ROT':
      rot = c.extract_valid_values()
    elif c.label() == 'LP':
      lp = c.extract_valid_values()

  # extract x, y, z positions for this in image address space

  xyz = []

  dx = b0.detlm()[1]
  dz = b0.phiend() - b0.phistt()
  z0 = b0.phistt()

  for x, y, r in zip(xdet, ydet, rot):
    _x = y
    _y = x
    _z = (r - z0) / dz

    xyz.append((_x, _y, _z))

  i = i_lp / lp

  print('Reference: %d observations' % len(hkl))
  return hkl, i, sigi, xyz

def integrate_mtz_to_unit_cell(integrate_mtz):
  '''Generate a cctbx unit_cell from an integrate_mtz file.'''

  from iotbx import mtz
  m = mtz.object(integrate_mtz)
  for c in m.crystals():
    return c.unit_cell()

  raise RuntimeError('unit cell not found')

def pull_calculated(integrate_pkl):
  from dials.array_family import flex # import dependency
  import cPickle as pickle
  #import cPickle as pickle
  import math

  #table = pickle.load(open('integrated.pickle', 'rb'))

  r_list = pickle.load(open(integrate_pkl, 'rb'))

  strong_reflections = []

  for r in r_list:
    if r.intensity > math.sqrt(r.intensity_variance):
      strong_reflections.append(r)

  del(r_list)

  hkl = []
  i = []
  sigi = []
  xyz = []

  for r in strong_reflections:
    if not r.is_valid():
      continue
    hkl.append(r.miller_index)
    i.append(r.intensity)
    sigi.append(math.sqrt(r.intensity_variance))
    x, y = r.image_coord_px
    z = r.frame_number
    xyz.append((x, y, z))

  print('Computed: %d observations' % len(hkl))
  return hkl, i, sigi, xyz

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

def compare_chunks(integrate_mtz, integrate_pkl, crystal_json, sweep_json):

  from cctbx.array_family import flex
  from annlib_ext import AnnAdaptor as ann_adaptor

  uc = integrate_mtz_to_unit_cell(integrate_mtz)

  rdx = derive_reindex_matrix(crystal_json, sweep_json, integrate_mtz)

  print(rdx)

  xhkl, xi, xsigi, xxyz = pull_reference(integrate_mtz)
  dhkl, di, dsigi, dxyz = pull_calculated(integrate_pkl)

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

  MOS = []
  DIALS = []
  HKL = []

  # perform the analysis
  for j, hkl in enumerate(dhkl):
    c = ann.nn[j]
    if hkl == tuple(rdx * xhkl[c]):
      MOS.append(xi[c])
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
  for hkl, mos, dials in zip(HKL, MOS, DIALS):
    sort_me.append((resolutions[hkl], mos, dials))

  sort_me.sort()
  sort_me.reverse()

  resolutions = [sm[0] for sm in sort_me]
  MOS = [sm[1] for sm in sort_me]
  DIALS = [sm[2] for sm in sort_me]

  # then extract the original observation structure

  print('Paired %d observations' % len(MOS))

  scale = sum(MOS) / sum(DIALS)

  chunks = [(i, i + 1000) for i in range(0, len(MOS), 1000)]

  ccs = []
  rs = []
  ss = []

  for chunk in chunks:
    mos = MOS[chunk[0]:chunk[1]]
    dials = DIALS[chunk[0]:chunk[1]]
    resols = resolutions[chunk[0]:chunk[1]]

    if len(mos) < 100:
      break

    c = cc(dials, mos)
    r, s = R(dials, mos)
    uncomment_me = '''
    print '%7d %4d %.3f %.3f %.3f %.3f %.3f' % (chunk[0], len(mos),
                                                min(resols), max(resols),
                                                c, r, s)'''
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
  pyplot.show()
  #pyplot.legend()
  #pyplot.savefig('plot-vs-mosflm.png')
  #pyplot.close()

  return

def get_dials_matrix(crystal_json):
  from dials.model.serialize import load
  crystal = load.crystal(crystal_json)
  return crystal.get_A()

def get_dials_coordinate_frame(sweep_json):
  from dials.model.serialize import load
  sweep = load.sweep(sweep_json)
  return sweep.get_beam().get_direction(), \
    sweep.get_goniometer().get_rotation_axis()

def get_mosflm_coordinate_frame(integrate_mtz):
  from iotbx import mtz
  from scitbx import matrix
  m = mtz.object(integrate_mtz)
  b = m.batches()[0]
  return matrix.col(b.source()), matrix.col(b.e1())

def integrate_mtz_to_A_matrix(integrate_mtz):
  from iotbx import mtz
  from cctbx.uctbx import unit_cell
  from scitbx import matrix
  m = mtz.object(integrate_mtz)
  b = m.batches()[0]
  u = matrix.sqr(b.umat()).transpose()
  c = unit_cell(tuple(b.cell()))
  f = matrix.sqr(c.fractionalization_matrix()).transpose()

  return (u * f)

def derive_reindex_matrix(crystal_json, sweep_json, integrate_mtz):
  '''Derive a reindexing matrix to go from the orientation matrix used
  for MOSFLM integration to the one used for DIALS integration.'''

  dA = get_dials_matrix(crystal_json)
  dbeam, daxis = get_dials_coordinate_frame(sweep_json)
  mbeam, maxis = get_mosflm_coordinate_frame(integrate_mtz)

  from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
  R = align_reference_frame(mbeam, dbeam, maxis, daxis)
  mA = R * integrate_mtz_to_A_matrix(integrate_mtz)

  # assert that this should just be a simple integer rotation matrix
  # i.e. reassignment of a, b, c so...

  from scitbx import matrix
  return matrix.sqr(map(int, map(round, (dA.inverse() * mA).elems)))

if __name__ == '__main__':
  import sys

  compare_chunks(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

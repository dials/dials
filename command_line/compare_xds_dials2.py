from __future__ import absolute_import, division
from __future__ import print_function
#from scipy.interpolate import griddata

# LIBTBX_SET_DISPATCHER_NAME dev.dials.compare_xds_dials2

import matplotlib
matplotlib.use('WXAgg')

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

#    if f_tokens[3] / f_tokens[4] < 5 or f_tokens[3] / f_tokens[4] > 10:
#      continue

    hkl.append(_hkl)

    # need to scale by PEAK stored value to get comparison value

    peak = 0.01 * f_tokens[9]
    i.append(f_tokens[3] * peak / f_tokens[8])
    sigi.append(f_tokens[4] * peak / f_tokens[8])
    xyz.append(tuple(f_tokens[5:8]))
    lp.append(f_tokens[8])

  print('Reference: %d observations' % len(hkl))
  return hkl, i, sigi, xyz, lp

def get_dials_matrix(crystal_json):
  from dials.model.serialize import load
  crystal = load.crystal(crystal_json)
  return crystal.get_A()

def get_dials_coordinate_frame(sweep_json):
  from dials.model.serialize import load
  sweep = load.sweep(sweep_json)
  return sweep.get_beam().get_direction(), \
    sweep.get_goniometer().get_rotation_axis()

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
    raise RuntimeError('coordinate frame information not found')

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
    raise RuntimeError('unit cell vectors not found')

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

  raise RuntimeError('unit cell not found')

def pull_calculated(integrate_pkl):
  from dials.array_family import flex # import dependency
  import cPickle as pickle
  import math

  r_list = pickle.load(open(integrate_pkl, 'rb'))

  strong_reflections = []

  for r in r_list:
    if not r.is_valid():
      continue
    if r.intensity ** 2 < r.intensity_variance:
      continue
    if r.intensity <= 0.0:
      continue
    strong_reflections.append(r)

  del(r_list)

  hkl = []
  i = []
  sigi = []
  xyz = []
  lp = []

  for r in strong_reflections:
    hkl.append(r.miller_index)
#    i.append(r.corrected_intensity)
#    sigi.append(math.sqrt(r.corrected_intensity_variance))
    i.append(r.intensity)
    sigi.append(math.sqrt(r.intensity_variance))
    lp.append(r.corrected_intensity / r.intensity)
    x, y = r.image_coord_px
    z = r.frame_number
    xyz.append((x, y, z))

  print('Computed: %d observations' % len(hkl))
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

  var = sum([(scale * c - o) ** 2 for c, o in zip(calc, obs)]) / len(calc)

  return sum([math.fabs(math.fabs(o) - math.fabs(scale * c)) \
              for c, o in zip(calc, obs)]) / \
              sum([math.fabs(o) for o in obs]), scale, var

def meansd(values):
  import math
  mean = sum(values) / len(values)
  var = sum([(v - mean) ** 2 for v in values]) / len(values)
  return mean, math.sqrt(var)

def compare_chunks(integrate_hkl, integrate_pkl, crystal_json, sweep_json,
                   d_min = 0.0):

  from cctbx.array_family import flex
  from annlib_ext import AnnAdaptor as ann_adaptor
  from dials.model.serialize import load
  sweep = load.sweep(sweep_json)

  rdx = derive_reindex_matrix(crystal_json, sweep_json, integrate_hkl)

  print('Reindex matrix:\n%d %d %d\n%d %d %d\n%d %d %d' % (rdx.elems))

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
  XYZ = []
  SIGMA_XDS = []
  SIGMA_DIALS = []
  XLP = []
  DLP = []

  # perform the analysis
  for j, hkl in enumerate(dhkl):
    c = ann.nn[j]
    if hkl == tuple(rdx * xhkl[c]):
      XDS.append(xi[c])
      DIALS.append(di[j])
      HKL.append(hkl)
      XYZ.append(dxyz[j])
      SIGMA_XDS.append(xsigi[c])
      SIGMA_DIALS.append(dsigi[j])
      XLP.append(xlp[c])
      DLP.append(dlp[j])

  print("Found %d matches" % len(XDS))

  compare = CompareIntensity(sweep, uc, HKL, XYZ, XDS, DIALS, SIGMA_XDS, SIGMA_DIALS, XLP, DLP)
#  compare.plot_scale_factor_vs_resolution()
#  compare.plot_scale_factor_vs_frame_number()
  compare.plot_chunked_statistics_vs_resolution()
  compare.plot_chunked_statistics_vs_frame_number()
  compare.plot_chunked_statistics_vs_i_over_sigma()
  compare.plot_chunked_i_over_sigma_vs_frame_number()
  compare.plot_chunked_resolution_vs_frame_number()
  compare.plot_chunked_lp_vs_frame_number()
  compare.plot_scale_vs_x_y()
  compare.plot_scale_vs_i_over_sigma()

def derive_reindex_matrix(crystal_json, sweep_json, integrate_hkl):
  '''Derive a reindexing matrix to go from the orientation matrix used
  for XDS integration to the one used for DIALS integration.'''

  dA = get_dials_matrix(crystal_json)
  dbeam, daxis = get_dials_coordinate_frame(sweep_json)
  xbeam, xaxis = get_xds_coordinate_frame(integrate_hkl)

  # want to align XDS -s0 vector...
  from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
  R = align_reference_frame(- xbeam, dbeam, xaxis, daxis)
  xA = R * integrate_hkl_to_A_matrix(integrate_hkl)

  # assert that this should just be a simple integer rotation matrix
  # i.e. reassignment of a, b, c so...

  from scitbx import matrix
  return matrix.sqr(map(int, map(round, (dA.inverse() * xA).elems)))


class CompareIntensity(object):

  def __init__(self, sweep, uc, hkl, xyz, i_xds, i_dials, sigma_xds, sigma_dials, xlp, dlp):
    self.sweep = sweep
    self.hkl = hkl
    self.xyz = xyz
    self.i_xds = i_xds
    self.sigma_xds = sigma_xds
    self.i_dials = i_dials
    self.sigma_dials = sigma_dials
    self.d = [uc.d(h) for h in self.hkl]
    self.scale = [d / x for d, x in zip(self.i_dials, self.i_xds)]
    self.xlp = xlp
    self.dlp = dlp

  def plot_scale_factor_vs_resolution(self):
    print("plot_scale_factor_vs_resolution")
    from matplotlib import pyplot
    index = [i for i in range(len(self.scale)) if abs(self.scale[i]) < 10]
    res = [self.d[i] for i in index]
    scale = [self.scale[i] for i in index]
    pyplot.xlabel('resolution')
    pyplot.ylabel('I_dials / I_xds')
    pyplot.title('I_dials / I_xds vs resolution')
    pyplot.scatter(res, scale)
    pyplot.savefig('plot-scale-vs-res.png')
    pyplot.close()

  def plot_scale_factor_vs_frame_number(self):
    print("plot_scale_factor_vs_frame_number")
    from matplotlib import pyplot
    index = [i for i in range(len(self.scale)) if abs(self.scale[i]) < 10]
    frame = [self.xyz[i][2] for i in index]
    scale = [self.scale[i] for i in index]
    pyplot.xlabel('resolution')
    pyplot.ylabel('I_dials / I_xds')
    pyplot.title('I_dials / I_xds vs resolution')
    pyplot.scatter(frame, scale)
    pyplot.savefig('plot-scale-vs-frame.png')
    pyplot.close()

  def plot_scale_factor_vs_x_and_y(self):
    pass

  def plot_chunked_statistics_vs_resolution(self):
    print("plot_chunked_statistics_vs_resolution")
    # Sort by resolution
    index = sorted(range(len(self.d)), key=lambda i: self.d[i])
    index.reverse()
    i_xds = [self.i_xds[i] for i in index]
    i_dials = [self.i_dials[i] for i in index]
    d = [self.d[i] for i in index]

    # Get stats for chunks
    chunks = [(i, i + 1000) for i in range(0, len(self.hkl), 1000)]
    ccs = []
    rs = []
    ss = []
    vs = []
    for chunk in chunks:
        xds = i_xds[chunk[0]:chunk[1]]
        dials = i_dials[chunk[0]:chunk[1]]
        resols = d[chunk[0]:chunk[1]]
        if len(xds) < 100:
          break
        c = cc(dials, xds)
        r, s, v = R(dials, xds)
        import math
        print('%7d %4d %.3f %.3f %.3f %.3f %.3f %.3f' % \
          (chunk[0], len(xds), min(resols), max(resols), c, r, s, math.sqrt(v)))
        ccs.append(c)
        rs.append(r)
        ss.append(s)
        vs.append(math.sqrt(v))
    chunks = [j for j in range(len(chunks))]
    chunks = chunks[:len(rs)]

    from matplotlib import pyplot
    pyplot.xlabel('Chunk')
    pyplot.ylabel('Statistic')
    pyplot.title('Statistics for 1000 reflection-pair chunks')
    pyplot.plot(chunks, ccs, label = 'CC')
    pyplot.plot(chunks, rs, label = 'R')
    pyplot.plot(chunks, ss, label = 'K')
#    pyplot.plot(chunks, vs, label = 'stddev')
    pyplot.legend()
    pyplot.savefig('plot-statistics-vs-res.png')
    pyplot.close()

  def plot_chunked_statistics_vs_frame_number(self):
    print("plot_chunked_statistics_vs_frame_number")
    # Sort by resolution
    index = sorted(range(len(self.xyz)), key=lambda i: self.xyz[i][2])
    i_xds = [self.i_xds[i] for i in index]
    i_dials = [self.i_dials[i] for i in index]
    frame = [self.xyz[i][2] for i in index]

    # Get stats for chunks
    chunks = [0]
    for i, z in enumerate(frame):
      if z > (len(chunks) * 10):
        chunks.append(i)

    chunks = list(zip(chunks[:-1], chunks[1:]))
    ccs = []
    rs = []
    ss = []
    vs = []
    for chunk in chunks:
        xds = i_xds[chunk[0]:chunk[1]]
        dials = i_dials[chunk[0]:chunk[1]]
        frames = frame[chunk[0]:chunk[1]]
        if len(xds) < 10:
          break
        c = cc(dials, xds)
        r, s, v = R(dials, xds)
        import math
        print('%7d %4d %.3f %.3f %.3f %.3f %.3f %.3f' % \
          (chunk[0], len(xds), min(frames), max(frames), c, r, s, math.sqrt(v)))
        ccs.append(c)
        rs.append(r)
        ss.append(s)
        vs.append(math.sqrt(v))
    chunks = [j for j in range(len(chunks))]
    chunks = chunks[:len(rs)]

    from matplotlib import pyplot
    pyplot.xlabel('Chunk')
    pyplot.ylabel('Statistic')
    pyplot.title('Statistics for 1000 reflection-pair chunks')
    pyplot.plot(chunks, ccs, label = 'CC')
    pyplot.plot(chunks, rs, label = 'R')
    pyplot.plot(chunks, ss, label = 'K')
    #pyplot.plot(chunks, vs, label = 'stddev')
    pyplot.legend()
    pyplot.savefig('plot-statistics-vs-frame.png')
    pyplot.close()

  def plot_chunked_statistics_vs_i_over_sigma(self):
    print("plot_chunked_statistics_vs_frame_number")
    # Sort by resolution
    i_over_s = [i / s for i, s in zip(self.i_xds, self.sigma_xds)]
    index = list(reversed(sorted(range(len(i_over_s)), key=lambda i: i_over_s[i])))
    i_xds = [self.i_xds[i] for i in index]
    i_dials = [self.i_dials[i] for i in index]
    i_over_s = [i_over_s[i] for i in index]

    # Get stats for chunks
    chunks = [(i, i + 1000) for i in range(0, len(self.hkl), 1000)]
    ccs = []
    rs = []
    ss = []
    for chunk in chunks:
        xds = i_xds[chunk[0]:chunk[1]]
        dials = i_dials[chunk[0]:chunk[1]]
        ios = i_over_s[chunk[0]:chunk[1]]
        if len(xds) < 100:
          break
        c = cc(dials, xds)
        r, s, v = R(dials, xds)
        print('%7d %4d %.3f %.3f %.3f %.3f %.3f' % \
          (chunk[0], len(xds), min(ios), max(ios), c, r, s))
        ccs.append(c)
        rs.append(r)
        ss.append(s)
    chunks = [j for j in range(len(chunks))]
    chunks = chunks[:len(rs)]

    from matplotlib import pyplot
    pyplot.xlabel('Chunk')
    pyplot.ylabel('Statistic')
    pyplot.title('Statistics for 1000 reflection-pair chunks')
    pyplot.plot(chunks, ccs, label = 'CC')
    pyplot.plot(chunks, rs, label = 'R')
    pyplot.plot(chunks, ss, label = 'K')
    pyplot.legend()
    pyplot.savefig('plot-statistics-vs-i_over_s.png')
    pyplot.close()

  def plot_scale_vs_i_over_sigma(self):
    print("plot_scale_vs_i_over_sigma")
    # Sort by resolution
    i_over_s = [i / s for i, s in zip(self.i_xds, self.sigma_xds)]
    index = list(reversed(sorted(range(len(i_over_s)), key=lambda i: i_over_s[i])))
    i_xds = [self.i_xds[i] for i in index]
    i_dials = [self.i_dials[i] for i in index]
    i_over_s = [i_over_s[i] for i in index]

    x0 = min(i_over_s)
    x1 = max(i_over_s)
    n = int(len(i_over_s) / 1000)
    step = (x1 - x0) / n

    chunk = [x0 + i * step for i in range(n)]
    scale = [0 for i in range(n)]
    count = [0 for i in range(n)]
    for ios, i_x, i_d in zip(i_over_s, i_xds, i_dials):
      i = int((ios - x0) / step)
      if i >= 0 and i < len(chunk):
        scale[i] += i_x / i_d
        count[i] += 1
    for i in range(n):
      if count[i] > 0:
        scale[i] /= count[i]

    from matplotlib import pyplot
    pyplot.xlabel('I/sig(I)')
    pyplot.ylabel('Scale')
    pyplot.plot(chunk, scale)
    pyplot.savefig('plot-scale-vs-i_over_s.png')
    pyplot.close()

    pyplot.xlabel('I/sig(I)')
    pyplot.ylabel('Num')
    pyplot.plot(chunk, count)
    pyplot.savefig('plot-num-vs-i_over_s.png')
    pyplot.close()

  def plot_chunked_i_over_sigma_vs_frame_number(self):
    print("plot_chunked_i_over_sigma__vs_frame_number")
    # Sort by resolution
    index = sorted(range(len(self.xyz)), key=lambda i: self.xyz[i][2])
    i_over_s = [i / s for i, s in zip(self.i_xds, self.sigma_xds)]
    i_over_s = [i_over_s[i] for i in index]
    frame = [self.xyz[i][2] for i in index]

    # Get stats for chunks
    chunks = [0]
    for i, z in enumerate(frame):
      if z > (len(chunks) * 10):
        chunks.append(i)

    chunks = list(zip(chunks[:-1], chunks[1:]))
    mean_i_over_sigma = []
    for chunk in chunks:
        ios = i_over_s[chunk[0]:chunk[1]]
        frames = frame[chunk[0]:chunk[1]]
        if len(ios) < 10:
          break
        mios = sum(ios) / len(ios)
        print('%7d %4d %.3f %.3f %.3f' % \
          (chunk[0], len(ios), min(frames), max(frames), mios))
        mean_i_over_sigma.append(mios)
    chunks = [j for j in range(len(chunks))]
    chunks = chunks[:len(mean_i_over_sigma)]

    from matplotlib import pyplot
    pyplot.xlabel('Chunk')
    pyplot.ylabel('Statistic')
    pyplot.title('Statistics for 1000 reflection-pair chunks')
    pyplot.plot(chunks, mean_i_over_sigma)
    pyplot.savefig('plot-i-over-sigma-vs-frame.png')
    pyplot.close()

  def plot_chunked_resolution_vs_frame_number(self):
    print("plot_chunked_resolution__vs_frame_number")
    # Sort by resolution
    index = sorted(range(len(self.xyz)), key=lambda i: self.xyz[i][2])
    d = [self.d[i] for i in index]
    frame = [self.xyz[i][2] for i in index]

    # Get stats for chunks
    chunks = [0]
    for i, z in enumerate(frame):
      if z > (len(chunks) * 10):
        chunks.append(i)

    chunks = list(zip(chunks[:-1], chunks[1:]))
    mean_d = []
    for chunk in chunks:
        dd = d[chunk[0]:chunk[1]]
        frames = frame[chunk[0]:chunk[1]]
        if len(dd) < 10:
          break
        md = sum(dd) / len(dd)
        print('%7d %4d %.3f %.3f %.3f' % \
          (chunk[0], len(dd), min(frames), max(frames), md))
        mean_d.append(md)
    chunks = [j for j in range(len(chunks))]
    chunks = chunks[:len(mean_d)]

    from matplotlib import pyplot
    pyplot.xlabel('Chunk')
    pyplot.ylabel('Statistic')
    pyplot.title('Statistics for 1000 reflection-pair chunks')
    pyplot.plot(chunks, mean_d)
    pyplot.savefig('plot-resolution-vs-frame.png')
    pyplot.close()

  def plot_chunked_lp_vs_frame_number(self):
    print("plot_chunked_lp__vs_frame_number")
    # Sort by resolution
    index = sorted(range(len(self.xyz)), key=lambda i: self.xyz[i][2])
    xlp = [self.xlp[i] for i in index]
    dlp = [self.dlp[i] for i in index]
    frame = [self.xyz[i][2] for i in index]

    # Get stats for chunks
    chunks = [0]
    for i, z in enumerate(frame):
      if z > (len(chunks) * 10):
        chunks.append(i)

    chunks = list(zip(chunks[:-1], chunks[1:]))
    ss = []
    for chunk in chunks:
        XLP = xlp[chunk[0]:chunk[1]]
        DLP = dlp[chunk[0]:chunk[1]]
        frames = frame[chunk[0]:chunk[1]]
        if len(XLP) < 10:
          break
        r, s, v0 = R(DLP, XLP)
        print('%7d %4d %.3f %.3f %.3f' % \
          (chunk[0], len(xlp), min(frames), max(frames), s))
        ss.append(s)
    chunks = [j for j in range(len(chunks))]
    chunks = chunks[:len(ss)]

    from matplotlib import pyplot
    pyplot.xlabel('Chunk')
    pyplot.ylabel('Statistic')
    pyplot.title('scale LP for 1000 reflection-pair chunks')
    pyplot.plot(chunks, ss)
    pyplot.savefig('plot-lp-vs-frame.png')
    pyplot.close()

  def plot_scale_vs_x_y(self):
    from scitbx.array_family import flex
    from math import ceil
    print("Getting scale")
    points = [(int(xyz[1] / 8), int(xyz[0] / 8)) for xyz in self.xyz]
    scale = [x / d for x, d in zip(self.i_xds, self.i_dials)]

    print("Creating Grid")
    image_size = self.sweep.get_detector()[0].get_image_size()[::-1]
    image_size = (int(ceil(image_size[0] / 8)), int(ceil(image_size[1] / 8)))
    grid = flex.double(flex.grid(image_size))
    count = flex.int(flex.grid(image_size))
    for p, s in zip(points, scale):
      grid[p] += s
      count[p] += 1
    for i in range(len(grid)):
      if count[i] > 0:
        grid[i] /= count[i]

    #grid_points = [(j,i) for j in range(image_size[0]) for i in range(image_size[1])]

    #grid = griddata(points, scale, grid_points)
    #grid.shape = image_size
    from matplotlib import pyplot
    fig, ax = pyplot.subplots()
    pyplot.title('scale vs x/y')
    cax = pyplot.imshow(grid.as_numpy_array())
    cbar = fig.colorbar(cax)
    pyplot.savefig('plot-scale-vs-xy.png')
    pyplot.close()


if __name__ == '__main__':
  import sys
  if len(sys.argv) < 5:
    raise RuntimeError( \
      '%s INTEGRATE.HKL integrate.pickle crystal.json sweep.json [dmin]' % \
      sys.argv[0])

  if len(sys.argv) == 5:
    compare_chunks(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
  else:
    compare_chunks(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                   d_min = float(sys.argv[5]))

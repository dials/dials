
from __future__ import division
from dials.array_family import flex

class GridIndex(object):

  def __init__(self, experiment, cs, offset, step_size):
    from dials.algorithms.reflection_basis import CoordinateSystem
    self.beam = experiment.beam
    self.detector = experiment.detector
    self.goniometer = experiment.goniometer
    self.scan = experiment.scan
    self.cs = cs
    self.phi = self.cs.phi()
    self.s1 = self.cs.s1()
    self.offset = offset
    self.step_size = step_size

  def __call__(self, k, j, i):
    from scitbx import matrix
    s1 = matrix.col(self.detector[0].get_pixel_lab_coord((i, j)))
    s1 = s1.normalize() * matrix.col(self.beam.get_s0()).length()
    e1, e2 = self.cs.from_beam_vector(s1)
    phi = self.scan.get_angle_from_array_index(k, deg=False)
    e3 = self.cs.from_rotation_angle(phi)
    gk = self.offset + e3 / self.step_size[0]
    gj = self.offset + e2 / self.step_size[1]
    gi = self.offset + e1 / self.step_size[2]
    return (int(gk), int(gj), int(gi))


def test_transform(experiment, shoebox, s1, phi, grid_size, sigma_m, sigma_b,
                   n_sigma, ndiv=5):
  from dials.algorithms.reflection_basis.transform import MapFramesForward
  from dials.algorithms.reflection_basis import CoordinateSystem

  bbox = shoebox.bbox
  data = shoebox.data
  mask = shoebox.mask
  x0, x1, y0, y1, z0, z1 = bbox

  step_size =  (
    n_sigma * sigma_m / (grid_size + 0.5),
    n_sigma * sigma_b / (grid_size + 0.5),
    n_sigma * sigma_b / (grid_size + 0.5))

  grid_range = 2*grid_size+1
  offset = grid_size + 0.5
  grid = flex.double(flex.grid(grid_range, grid_range, grid_range))

  m2 = experiment.goniometer.get_rotation_axis()
  s0 = experiment.beam.get_s0()
  cs = CoordinateSystem(m2, s0, s1, phi)

  grid_index = GridIndex(experiment, cs, offset, step_size)

  #phi0, dphi = experiment.scan.get_oscillation(deg=False)
  #map_frames = MapFramesForward(phi0, dphi, sigma_m, n_sigma, grid_size)
  #zfraction = map_frames((z0, z1), phi, cs.zeta())

  #for j in range(zfraction.all()[0]):
    #print ' '.join(str(zfraction[j,i]) for i in range(zfraction.all()[1]))

  #for j in range(zfraction.all()[0]):
    #print "Fraction %d: %f" % (j, flex.sum(zfraction[j:j+1,:]))

  fraction = 1.0 / (ndiv * ndiv * ndiv)
  for k in range(data.all()[0]):
    for j in range(data.all()[1]):
      for i in range(data.all()[2]):
        if mask[k,j,i]:
          value = data[k,j,i] * fraction
          for kk in range(0, ndiv):
            for jj in range(0, ndiv):
              for ii in range(0, ndiv):
                kkk = z0 + k + 1.0 * (kk + 0.5) / ndiv
                jjj = y0 + j + 1.0 * (jj + 0.5) / ndiv
                iii = x0 + i + 1.0 * (ii + 0.5) / ndiv
                gk, gj, gi = grid_index(kkk,jjj,iii)
                if (gk >= 0 and gk < grid_range and
                    gj >= 0 and gj < grid_range and
                    gi >= 0 and gi < grid_range):
                  grid[gk,gj,gi] += value

  #fraction = 1.0 / (ndiv * ndiv)
  #for k in range(data.all()[0]):
    #for j in range(data.all()[1]):
      #for i in range(data.all()[2]):
        #if mask[k,j,i]:
          #value = data[k,j,i] * fraction
          #for jj in range(0, ndiv):
            #for ii in range(0, ndiv):
              #jjj = y0 + j + 1.0 * (jj + 0.5) / ndiv
              #iii = x0 + i + 1.0 * (ii + 0.5) / ndiv
              #gj, gi = grid_index(jjj,iii)
              #if (gj >= 0 and gj < grid_range and
                  #gi >= 0 and gi < grid_range):
                #for gk in range(0, 2*grid_size+1):
                  #grid[gk,gj,gi] += value * zfraction[k,gk]

  return grid

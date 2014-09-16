
def rs_coords(bbox, cs, beam, detector, scan):
  from dials.array_family import flex
  from scitbx import matrix
  x0, x1, y0, y1, z0, z1 = bbox
  xs = bbox[1]-bbox[0]
  ys = bbox[3]-bbox[2]
  zs = bbox[5]-bbox[4]
  coords = flex.vec3_double(flex.grid(zs+1, ys+1, xs+1))
  for k in range(z0, z1+1):
    phid = scan.get_angle_from_array_index(k)
    e3 = cs.from_rotation_angle(phid)
    for j in range(y0, y1+1):
      for i in range(x0, x1+1):
        s1d = matrix.col(detector[0].get_pixel_lab_coord((i, j))).normalize()
        s1d *= matrix.col(beam.get_s0()).length()
        e1, e2 = cs.from_beam_vector(s1d)
        zz = k - z0
        yy = j - y0
        xx = i - x0
        coords[zz,yy,xx] = (e1, e2, e3)
  return coords

def to_s1_and_phi(e1, e2, e3, cs):
  from scitbx import matrix
  phid = cs.to_rotation_angle_fast(e3)
  s1 = matrix.col(cs.s1())
  s1l = s1.length()
  e1_axis = matrix.col(cs.e1_axis())*s1l
  e2_axis = matrix.col(cs.e2_axis())*s1l
  s1d = (e1_axis * e1 + e2_axis * e2) + s1
  return s1d, phid

def px_coords(e, cs, detector, scan):
  xyz = flex.vec3_double(flex.grid(e.all()))
  for k in range(e.all()[0]):
    for j in range(e.all()[1]):
      for i in range(e.all()[2]):
        e1, e2, e3 = e[k,j,i]
        s1d, phid = to_s1_and_phi(e1, e2, e3, cs)
        # s1d = cs.to_beam_vector((e1, e2))
        # phid = cs.to_rotation_angle(e3)
        x, y = detector[0].get_ray_intersection_px(s1d)
        z = scan.get_array_index_from_angle(phid)
        xyz[k,j,i] = (x, y, z)
  return xyz


def drizzle_test(sbox, beam_vectors, angles, experiment):
  from dials.algorithms.reflection_basis import CoordinateSystem
  # from dials.algorithms.image.drizzle import drizzle

  beam = experiment.beam
  detector = experiment.detector
  goniometer = experiment.goniometer
  scan = experiment.scan

  s0 = beam.get_s0()
  m2 = goniometer.get_rotation_axis()

  s1 = beam_vectors[0]
  phi = angles[0]
  cs_target = CoordinateSystem(m2, s0, s1, phi)
  from matplotlib import pylab
  fig = pylab.figure()
  ax = fig.gca()
  colour = ['black', 'blue', 'red', 'green', 'yellow', 'cyan', 'magenta']
  output = flex.double(flex.grid(sbox[0].data.all()))
  for r in range(len(sbox)):
    s1 = beam_vectors[r]
    phi = angles[r]
    cs = CoordinateSystem(m2, s0, s1, phi)
    e = rs_coords(sbox[r].bbox, cs, beam, detector, scan)
    e1, e2, e3 = e.parts()
    p = px_coords(e, cs_target, detector, scan)
    x, y, z = p.parts()
    pylab.scatter(x, y, color=colour[r])
    input = sbox[r].data.all()
    x -= sbox[0].bbox[0]
    y -= sbox[0].bbox[2]
    z -= sbox[0].bbox[4]
    # coords = flex.vec3_double(x, y, z)
    # drizzle(output, input, coords)

  ax.set_xticks(range(sbox[0].bbox[0], sbox[0].bbox[1]+1))
  ax.set_yticks(range(sbox[0].bbox[2], sbox[0].bbox[3]+1))
  pylab.grid()
  pylab.show()




if __name__ == "__main__":

  import sys
  from dials.array_family import flex
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory

  exlist = ExperimentListFactory.from_json_file(sys.argv[1])
  rlist = flex.reflection_table.from_pickle(sys.argv[2])

  I = rlist['intensity.sum.value']
  indices = flex.size_t(list(reversed(
    sorted(range(len(I)), key=lambda x:I[x])))[0:7])
  rlist = rlist.select(indices)

  sweep = exlist[0].imageset

  sbox = flex.shoebox(rlist['panel'], rlist['bbox'])
  s1 = rlist['s1']
  phi = rlist['xyzcal.mm'].parts()[2]

  print "Read pixels"
  for i in range(len(sbox)):
    print i
    x0, x1, y0, y1, z0, z1 = sbox[i].bbox
    data = sweep.to_array((z0, z1, y0, y1, x0, x1))
    sbox[i].data = data.as_double()

  profile = drizzle_test(sbox, s1, phi, exlist[0])

from __future__ import division

def run(args):
  from scitbx.array_family import flex
  from scitbx import matrix
  from dials.util.command_line import Importer
  from dials.algorithms.reflection_basis import zeta_factor

  importer = Importer(args, check_format=False)
  assert importer.datablocks is not None
  assert len(importer.datablocks) == 1
  datablock = importer.datablocks[0]
  imagesets = datablock.extract_imagesets()
  assert len(imagesets) == 1
  imageset = imagesets[0]

  detector = imageset.get_detector()
  beam = imageset.get_beam()
  goniometer = imageset.get_goniometer()
  assert goniometer is not None
  assert len(detector) == 1
  panel = detector[0]

  lab_coords = flex.vec3_double(flex.grid(panel.get_image_size()))

  for i in range(panel.get_image_size()[0]):
    for j in range(panel.get_image_size()[1]):
      lab_coords[i,j] = panel.get_lab_coord(panel.pixel_to_millimeter((i,j)))

  axis = matrix.col(goniometer.get_rotation_axis())
  s0 = matrix.col(beam.get_s0())
  s1 = (lab_coords.as_1d()/lab_coords.as_1d().norms()) * s0.length()
  s1_cross_s0 = s1.cross(flex.vec3_double(s1.size(), s0.elems))
  p_volume = flex.abs(s1_cross_s0.dot(axis.elems))
  p_volume.reshape(flex.grid(panel.get_image_size()))
  zeta = flex.abs(zeta_factor(axis.elems, s0.elems, s1.as_1d()))
  zeta.reshape(flex.grid(panel.get_image_size()))

  from matplotlib import pyplot
  pyplot.figure()
  pyplot.title('parallelepiped volume')
  CS = pyplot.contour(p_volume.matrix_transpose().as_numpy_array(), 10)
  pyplot.clabel(CS, inline=1, fontsize=10, fmt="%6.3f")
  pyplot.axes().set_aspect('equal')
  pyplot.show()
  pyplot.title('zeta factor')
  CS = pyplot.contour(zeta.matrix_transpose().as_numpy_array(), 10)
  pyplot.clabel(CS, inline=1, fontsize=10, fmt="%6.3f")
  pyplot.axes().set_aspect('equal')
  pyplot.show()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

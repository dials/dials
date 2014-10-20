from dials.array_family import flex
from dials.viewer.next_gen_viewer.simple_3D_slice_viewer import show_3d
if(__name__ == "__main__"):
  size_xy = 6

  data_xyz_flex = flex.double(flex.grid(size_xy, size_xy, size_xy),15)
  data_xyz_flex[1, 2, 2] = 15
  data_xyz_flex[2, 2, 2] = 20
  data_xyz_flex[3, 2, 2] = 25
  data_xyz_flex[4, 2, 2] = 20
  data_xyz_flex[5, 2, 2] = 15

  for frm in range(size_xy):
    for row in range(size_xy):
      for col in range(size_xy):
        data_xyz_flex[frm, row, col] += (row * 2 + col * 2 + frm * 2)


  show_3d(data_xyz_flex)

from __future__ import absolute_import, division
from __future__ import print_function
def build_prediction_matrix(hkl, mhkl, phkl, hmkl, hpkl, hkml, hkpl, tst=False):
  '''Build a prediction matrix around reflection hkl, coomputing dx / dh etc.'''

  from scitbx import matrix

  refl = hkl, mhkl, phkl, hmkl, hpkl, hkml, hkpl

  x = [r.image_coord_px[0] for r in refl]
  y = [r.image_coord_px[1] for r in refl]
  z = [r.frame_number for r in refl]

  d = 0.5 * matrix.sqr((x[2] - x[1], x[4] - x[3], x[6] - x[5],
                        y[2] - y[1], y[4] - y[3], y[6] - y[5],
                        z[2] - z[1], z[4] - z[3], z[6] - z[5]))

  if not tst:
    return d

  # now assess how accurate this d matrix is

  for j, r in enumerate(refl[1:]):
    dh = matrix.col((r.miller_index[0] - hkl.miller_index[0],
                     r.miller_index[1] - hkl.miller_index[1],
                     r.miller_index[2] - hkl.miller_index[2]))
    dxyz = d * dh
    dx = x[0] + dxyz[0] - x[j + 1]
    dy = y[0] + dxyz[1] - y[j + 1]
    dz = z[0] + dxyz[2] - z[j + 1]

    print(j, dx, dy, dz)

  return d

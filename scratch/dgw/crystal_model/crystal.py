from __future__ import division
from scitbx import matrix
from cctbx.uctbx import unit_cell
from cctbx.sgtbx import space_group, space_group_symbols
from cctbx.crystal_orientation import crystal_orientation

class Crystal:
  """Simple model for the crystal lattice geometry and symmetry"""

  def __init__(self, real_space_a, real_space_b, real_space_c, sg = 1):

    self._sg = space_group(space_group_symbols(sg).hall())

    # setting matrix at initialisation
    A = matrix.sqr(real_space_a.elems +  real_space_b.elems + \
                   real_space_c.elems).inverse()

    # unit cell
    self.set_unit_cell(real_space_a, real_space_b, real_space_c)

    # reciprocal space orthogonalisation matrix (is the transpose of the
    # real space fractionalisation matrix, see http://goo.gl/H3p1s)
    self._update_B()

    # initial orientation matrix
    self._U = A * self._B.inverse()

  def set_unit_cell(self, real_space_a, real_space_b, real_space_c):
    cell = (real_space_a.length(),
            real_space_b.length(),
            real_space_c.length(),
            real_space_b.angle(real_space_c, deg = True),
            real_space_c.angle(real_space_a, deg = True),
            real_space_a.angle(real_space_b, deg = True))
    self._uc = unit_cell(cell)
    self._update_B()

  def _update_B(self):
    self._B = matrix.sqr(self._uc.fractionalization_matrix()).transpose()

  def set_U(self, U):

    # check U is a rotation matrix.
    assert(U.is_r3_rotation_matrix())
    self._U = U

  def get_U(self):
    return self._U

  def get_B(self):
    return self._B

  def set_B(self, B):

    # also set the unit cell
    co = crystal_orientation(B,True)
    self._uc = co.unit_cell()
    self._B = matrix.sqr(self._uc.fractionalization_matrix()).transpose()

  def get_unit_cell(self):
    return self._uc

  def get_space_group(self):
    return self._sg

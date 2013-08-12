from __future__ import division
import math
from cStringIO import StringIO
from libtbx.test_utils import approx_equal, show_diff
from scitbx import matrix
from cctbx import sgtbx, uctbx
from dials.model.experiment.crystal_model import Crystal

def exercise_crystal_model():

  real_space_a = matrix.col((10,0,0))
  real_space_b = matrix.col((0,11,0))
  real_space_c = matrix.col((0,0,12))
  model = Crystal(real_space_a=(10,0,0),
                  real_space_b=(0,11,0),
                  real_space_c=(0,0,12),
                  space_group_symbol="P 1")
  assert isinstance(model.get_unit_cell(), uctbx.unit_cell)
  assert model.get_unit_cell().parameters() == (
    10.0, 11.0, 12.0, 90.0, 90.0, 90.0)
  assert approx_equal(model.get_A(), (1/10, 0, 0, 0, 1/11, 0, 0, 0, 1/12))
  assert approx_equal(model.get_A().inverse(), (10, 0, 0, 0, 11, 0, 0, 0, 12))
  assert approx_equal(model.get_B(), model.get_A())
  assert approx_equal(model.get_U(), (1, 0, 0, 0, 1, 0, 0, 0, 1))
  assert approx_equal(model.get_real_space_vectors(),
                      (real_space_a, real_space_b, real_space_c))
  assert approx_equal(model.get_mosaicity(), 0)

  model2 = Crystal(real_space_a=(10,0,0),
                  real_space_b=(0,11,0),
                  real_space_c=(0,0,12),
                  space_group_symbol="P 1")
  assert model == model2
  model2.set_mosaicity(0.01)
  assert model != model2
  # rotate 45 degrees about x-axis
  R1 = matrix.sqr((1, 0, 0,
                   0, math.cos(math.pi/4), -math.sin(math.pi/4),
                   0, math.sin(math.pi/4), math.cos(math.pi/4)))
  # rotate 30 degrees about y-axis
  R2 = matrix.sqr((math.cos(math.pi/6), 0, math.sin(math.pi/6),
                   0, 1, 0,
                   -math.sin(math.pi/6), 0, math.cos(math.pi/6)))
  # rotate 60 degrees about z-axis
  R3 = matrix.sqr((math.cos(math.pi/3), -math.sin(math.pi/3), 0,
                   math.sin(math.pi/3), math.cos(math.pi/3), 0,
                   0, 0, 1))
  R = R1 * R2 * R3
  model.set_U(R)
  # B is unchanged
  assert approx_equal(model.get_B(), (1/10, 0, 0, 0, 1/11, 0, 0, 0, 1/12))
  assert approx_equal(model.get_U(), R)
  assert approx_equal(model.get_A(), model.get_U() * model.get_B())
  a_, b_, c_ = model.get_real_space_vectors()
  assert approx_equal(a_, R * real_space_a)
  assert approx_equal(b_, R * real_space_b)
  assert approx_equal(c_, R * real_space_c)
  s = StringIO()
  print >> s, model
  assert not show_diff(s.getvalue(), """\
Crystal:
    Unit cell: (10.000, 11.000, 12.000, 90.000, 90.000, 90.000)
    Space group: P 1
    U matrix:  {{ 0.4330, -0.7500,  0.5000},
                { 0.7891,  0.0474, -0.6124},
                { 0.4356,  0.6597,  0.6124}}
    B matrix:  {{ 0.1000,  0.0000,  0.0000},
                {-0.0000,  0.0909,  0.0000},
                {-0.0000,  0.0000,  0.0833}}
    A = UB:    {{ 0.0433, -0.0682,  0.0417},
                { 0.0789,  0.0043, -0.0510},
                { 0.0436,  0.0600,  0.0510}}

""")
  model.set_B((1/12, 0, 0, 0, 1/12, 0, 0, 0, 1/12))
  assert approx_equal(model.get_unit_cell().parameters(),
                      (12, 12, 12, 90, 90, 90))

  model3 = Crystal(
    real_space_a=(10,0,0),
    real_space_b=(0,11,0),
    real_space_c=(0,0,12),
    space_group=sgtbx.space_group_info("P 222").group(),
    mosaicity=0.01)
  assert model3.get_space_group().type().hall_symbol() == " P 2 2"
  assert model != model3


def run():
  exercise_crystal_model()

if __name__ == '__main__':
  from libtbx.utils import show_times_at_exit
  show_times_at_exit()
  run()

from __future__ import division
from scitbx import matrix
from dials.scratch.dgw.crystal_model import Crystal

if __name__ == '__main__':

  #### Imports

  import random

  #### Local functions

  def random_direction_close_to(vector):
    return vector.rotate_around_origin(matrix.col(
                (random.random(),
                 random.random(),
                 random.random())).normalize(),
                 random.gauss(0, 5.0),  deg = True)

  #### Unit tests

  # instantiation
  a = random.uniform(10,50) * random_direction_close_to(matrix.col((1, 0, 0)))
  b = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 1, 0)))
  c = random.uniform(10,50) * random_direction_close_to(matrix.col((0, 0, 1)))

  xl = Crystal(a, b, c, space_group_symbol="P 1")

  print "OK"

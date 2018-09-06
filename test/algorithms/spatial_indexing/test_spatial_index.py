from __future__ import absolute_import, division, print_function

import random

from dials.algorithms.spatial_indexing import make_spatial_index
import pytest

@pytest.fixture
def vec2_double():
  from scitbx.array_family import flex
  num = 100
  x0, x1, y0, y1 = 0, 100, 0, 100
  result = flex.vec2_double(num)
  for i in range(num):
    result[i] = (random.uniform(x0, x1), random.uniform(y0, y1))
  return result

@pytest.fixture
def vec3_double():
  from scitbx.array_family import flex
  num = 100
  x0, x1, y0, y1, z0, z1 = 0, 100, 0, 100, 0, 100
  result = flex.vec3_double(num)
  for i in range(num):
    result[i] = (random.uniform(x0, x1), random.uniform(y0, y1), random.uniform(z0, z1))
  return result

def inside2d(rg, a):
  inside = (rg[0] <= a[0] and rg[1] >= a[0] and
            rg[2] <= a[1] and rg[3] >= a[1])
  return inside

def inside3d(rg, a):
  inside = (rg[0] <= a[0] and rg[1] >= a[0] and
            rg[2] <= a[1] and rg[3] >= a[1] and
            rg[4] <= a[2] and rg[5] >= a[2])
  return inside

def test_vec2_double(vec2_double):
  index = make_spatial_index(vec2_double)
  x, y = zip(*vec2_double)
  x0, x1, y0, y1 = (int(min(x)), int(max(x)),
                    int(min(y)), int(max(y)))
  for i in range(10):
    xx0 = random.randint(x0, x1)
    xx1 = random.randint(x0, x1)
    yy0 = random.randint(y0, y1)
    yy1 = random.randint(y0, y1)
    rg = (min(xx0, xx1), max(xx0, xx1)+1,
          min(yy0, yy1), max(yy0, yy1)+1)
    idx = index.query_range(rg)
    for j in idx:
      assert inside2d(rg, vec2_double[j])
    for j in set(range(len(vec2_double))).difference(set(idx)):
      assert not inside2d(rg, vec2_double[j])

def test_vec3_double(vec3_double):
  index = make_spatial_index(vec3_double)
  x, y, z = zip(*vec3_double)
  x0, x1, y0, y1, z0, z1 = (int(min(x)), int(max(x)),
                            int(min(y)), int(max(y)),
                            int(min(z)), int(max(z)))
  for i in range(10):
    xx0 = random.randint(x0, x1)
    xx1 = random.randint(x0, x1)
    yy0 = random.randint(y0, y1)
    yy1 = random.randint(y0, y1)
    zz0 = random.randint(y0, y1)
    zz1 = random.randint(y0, y1)
    rg = (min(xx0, xx1), max(xx0, xx1)+1,
          min(yy0, yy1), max(yy0, yy1)+1,
          min(zz0, zz1), max(zz0, zz1)+1)
    idx = index.query_range(rg)
    for j in idx:
      assert inside3d(rg, vec3_double[j])
    for j in set(range(len(vec3_double))).difference(set(idx)):
      assert not inside3d(rg, vec3_double[j])

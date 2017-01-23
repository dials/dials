
from __future__ import division

class Test:

  def __init__(self):
    pass

  def run(self):
    from dials.algorithms.image.filter import manhatten_distance
    from scitbx.array_family import flex
    from math import sqrt

    data = flex.bool(flex.grid(100, 100), True)

    for j in range(100):
      for i in range(100):
        if sqrt((j - 50)**2 + (i - 50)**2) <= 10.5:
          data[j,i] = False

    distance = manhatten_distance(data)

    M = distance[1:-1,1:-1]
    D = data[1:-1,1:-1]

    selection = D.as_1d().select(M.as_1d() == 0)
    assert selection.all_eq(True)

    while True:
      N = data[0:-2,1:-1]
      S = data[2:,1:-1]
      E = data[1:-1,0:-2]
      W = data[1:-1,2:]
      selection = M.as_1d() == 1
      neighbours = (
        N.select(selection) |
        S.select(selection) |
        E.select(selection) |
        W.select(selection))
      assert (neighbours.all_eq(True))
      indices = flex.size_t(range(len(D))).select(selection)
      for i in indices:
        D[i] = True
      data[1:-1,1:-1] = D
      M -= 1
      if (M > 0).count(True) == 0:
        break

    # Test passed
    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

from dials.array_family import flex
from dials.algorithms.integration.profile import MLPoisson2

c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
b = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
s = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]

b = [bb/ sum(b) for bb in b]
s = [ss /sum(s) for ss in s]

result = MLPoisson2(flex.double(c), flex.double(b), flex.double(s), 1.0, 1.0,
                    1.0)


for i in range(10):
  result.step()
  print result.S1(), result.S2()

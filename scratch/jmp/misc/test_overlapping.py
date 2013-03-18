
from dials.model.data import Reflection, ReflectionList
from dials.algorithms.integration import find_overlapping_reflections
rl = ReflectionList()
r1 = Reflection()
r1.shoebox = (10, 20, 10, 20, 10, 20)
r2 = Reflection()
r2.shoebox = (15, 25, 15, 25, 15, 25)
r3 = Reflection()
r3.shoebox = (20, 30, 20, 30, 20, 30)
rl.append(r1)
rl.append(r2)
rl.append(r3)
overlapping = find_overlapping_reflections(rl)

for v in overlapping.vertices():
  print "Vertex: ", v, " => ", [a for a in overlapping.adjacent_vertices(v)] 


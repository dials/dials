


from dials.algorithms.integration import ReflectionExtractor
from dials.model.serialize import load
import os
path = '/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data'

sweep = load.sweep(os.path.join(path, 'sweep.json'))
crystal = load.crystal(os.path.join(path, 'crystal.json'))

extract = ReflectionExtractor(3)
refl = extract(sweep, crystal)


import copy

refl2 = copy.deepcopy(refl)


from dials.algorithms.reflection_basis import transform
from dials.util.command_line import Command

Command.start('Init transform')
trans = transform.Forward(sweep, crystal, 3, 4)
Command.end('Init transform')

Command.start('Transform')
for i in range(1):
  trans(refl)
Command.end('Transform')

Command.start('Init Transform')
spec = transform.TransformSpec(sweep, crystal, 3, 4)
Command.end('Init Transform')

Command.start('Transform')
for i in range(1):
  transform.forward_batch(spec, refl2)
Command.end('Transform')

assert(len(refl) == len(refl2))
for r1, r2 in zip(refl, refl2):
  assert(r1.is_valid() == r2.is_valid())
  if r1.is_valid():

    prof1 = r1.transformed_shoebox
    prof2 = r2.transformed_shoebox

    assert(len(prof1) == len(prof2))

    diff = abs(prof1 - prof2)
    assert(diff.all_lt(1e-7))

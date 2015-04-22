from cctbx import crystal, sgtbx, uctbx
from cctbx.sgtbx import lattice_symmetry
from cctbx.sgtbx.bravais_types import bravais_lattice

sg = sgtbx.space_group_info('C2').group()
uc = uctbx.unit_cell(
  parameters='102.965614947,79.2397736681,77.45579668,90.0,139.074859178,90.0')
cs = crystal.symmetry(unit_cell=uc, space_group=sg)
subgroups = lattice_symmetry.metric_subgroups(cs, max_delta=0.1)
for subgroup in subgroups.result_groups:
  bravais_t = bravais_lattice(
    group=subgroup['ref_subsym'].space_group())
  if bravais_t == 'mC':
    print subgroup['ref_subsym'].unit_cell()
    print subgroup['best_subsym'].as_reference_setting().unit_cell()

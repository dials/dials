from cctbx import crystal, sgtbx, uctbx
from cctbx.sgtbx import lattice_symmetry
from cctbx.sgtbx.bravais_types import bravais_lattice

sg = sgtbx.space_group_info('C2').group()
uc = uctbx.unit_cell(
  parameters='102.965614947,79.2397736681,77.45579668,90.0,139.074859178,90.0')
cs = crystal.symmetry(unit_cell=uc, space_group=sg)
print "Input cell:"
cs.show_summary()
print

print "Best cell:"
cs.best_cell().show_summary()
print

print "Reference settings:"
print "Best cell -> reference setting"
cs.best_cell().as_reference_setting().show_summary()
print "Input -> reference setting"
cs.as_reference_setting().show_summary()
print "Input -> primitive setting -> reference setting"
cs.primitive_setting().as_reference_setting().show_summary()
print "Best cell -> primitive setting -> reference setting"
print

print "Primitive settings:"
cs.best_cell().primitive_setting().as_reference_setting().show_summary()
print "Best cell -> primitive setting"
cs.best_cell().primitive_setting().show_summary()
print "Best cell -> reference setting -> primitive setting"
cs.best_cell().as_reference_setting().primitive_setting().show_summary()
print "Input -> reference setting -> primitive setting"
cs.as_reference_setting().primitive_setting().show_summary()
print "Input -> primitive setting"
cs.primitive_setting().show_summary()
print

print "Minimum/niggli cell:"
print "Input -> minimum cell"
cs.minimum_cell().show_summary()
print "Input -> niggli cell"
cs.niggli_cell().show_summary()
print "Best cell -> minimum cell"
cs.best_cell().minimum_cell().show_summary()
print "Best cell -> niggli cell"
cs.best_cell().niggli_cell().show_summary()
print "Reference setting -> minimum cell"
cs.as_reference_setting().minimum_cell().show_summary()
print "Reference setting -> niggli cell"
cs.as_reference_setting().niggli_cell().show_summary()
print

print "Reference settings (via minimum/niggli cell):"
print "Input -> minimum cell -> reference setting"
cs.minimum_cell().as_reference_setting().show_summary()
print "Input -> niggli cell -> reference setting"
cs.niggli_cell().as_reference_setting().show_summary()
print "Best cell -> minimum cell -> reference setting"
cs.best_cell().minimum_cell().as_reference_setting().show_summary()
print "Best cell -> niggli cell -> reference setting"
cs.best_cell().niggli_cell().as_reference_setting().show_summary()
print "Reference setting -> minimum cell -> reference setting"
cs.as_reference_setting().minimum_cell().as_reference_setting().show_summary()
print "Reference setting -> niggli cell -> reference setting"
cs.as_reference_setting().niggli_cell().as_reference_setting().show_summary()
print

subgroups = lattice_symmetry.metric_subgroups(cs, max_delta=0.1)
for subgroup in subgroups.result_groups:
  bravais_t = bravais_lattice(
    group=subgroup['ref_subsym'].space_group())
  if bravais_t == 'mC':
    print subgroup['ref_subsym'].unit_cell()
    print subgroup['best_subsym'].as_reference_setting().unit_cell()

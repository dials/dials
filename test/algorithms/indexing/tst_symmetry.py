from __future__ import absolute_import, division

def exercise_find_matching_symmetry():
  from cctbx import crystal, sgtbx, uctbx
  from dials.algorithms.indexing import symmetry

  def run_once(crystal_symmetry):
    cs = crystal_symmetry
    #print
    cs.show_summary()

    from cctbx.sgtbx import lattice_symmetry
    subgroups = lattice_symmetry.metric_subgroups(cs, max_delta=5)

    for op in ('x,y,z', 'z,x,y', 'y,z,x', '-x,z,y', 'y,x,-z', 'z,-y,x')[:]:
      cb_op = sgtbx.change_of_basis_op(op)

      uc_inp = cs.unit_cell().change_basis(cb_op)

      for ref_uc, ref_sg in [(cs.unit_cell(), cs.space_group()),
                             (None, cs.space_group())][:]:

        best_subgroup = symmetry.find_matching_symmetry(
          uc_inp,
          target_space_group=ref_sg)
        cb_op_inp_best = best_subgroup['cb_op_inp_best']

        assert uc_inp.change_basis(cb_op_inp_best).is_similar_to(
          cs.as_reference_setting().best_cell().unit_cell())

  # 1) known unit cell and space group
  #      -> cb op for test unit cell to best match for reference symmetry
  # 2) known unit cell
  #      -> cb op for test unit cell to best match for given unit cell
  # 3) known space group
  #      -> cb op for test unit cell to best match for given space group

  # Assumptions:
  # 1) test unit cell is minimum cell

  uc = uctbx.unit_cell("76, 115, 134, 90, 99.07, 90")
  sgi = sgtbx.space_group_info(symbol="I2")
  cs = crystal.symmetry(unit_cell=uc, space_group_info=sgi)
  uc_inp = uc.minimum_cell()
  cs_min = crystal.symmetry(
    unit_cell=cs.minimum_cell().unit_cell(), space_group=sgtbx.space_group())


  uc = uctbx.unit_cell("42,42,40,90,90,90")
  sgi = sgtbx.space_group_info(symbol="P41212")
  cs = crystal.symmetry(unit_cell=uc, space_group_info=sgi)
  cb_op = sgtbx.change_of_basis_op("c,a,b")
  uc_inp = uc.change_basis(cb_op)

  run_once(cs.change_basis(cb_op))

  from cctbx.sgtbx import bravais_types
  for symbol in bravais_types.acentric:
    sgi = sgtbx.space_group_info(symbol=symbol)
    uc = sgi.any_compatible_unit_cell(volume=1000)
    cs = crystal.symmetry(unit_cell=uc, space_group_info=sgi)
    cs = cs.niggli_cell().as_reference_setting().primitive_setting()
    run_once(cs)

def run():
  exercise_find_matching_symmetry()
  print "OK"

if __name__ == '__main__':
  run()

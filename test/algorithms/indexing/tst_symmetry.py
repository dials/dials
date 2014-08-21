from __future__ import division

def exercise_change_of_basis_op_to_best_cell():
  from cctbx import crystal, sgtbx, uctbx
  from dials.algorithms.indexing import symmetry

  uc = uctbx.unit_cell("42,42,40,90,90,90")
  sgi = sgtbx.space_group_info(symbol="P41212")
  cs = crystal.symmetry(unit_cell=uc, space_group_info=sgi)
  cb_op = sgtbx.change_of_basis_op("c,a,b")
  uc_inp = uc.change_basis(cb_op)

  cb_op_to_best_cell = symmetry.change_of_basis_op_to_best_cell(
    unit_cell=uc_inp,
    target_unit_cell=cs.unit_cell(),
    target_space_group=cs.space_group())
  assert str(cb_op_to_best_cell) == str(cb_op.inverse())
  assert cs.primitive_setting().unit_cell().is_similar_to(
    uc_inp.change_basis(cb_op_to_best_cell))
  print uc_inp
  print cb_op_to_best_cell

  cb_op_to_best_cell = symmetry.change_of_basis_op_to_best_cell(
    unit_cell=uc_inp,
    target_unit_cell=cs.unit_cell(),
    target_space_group=None)
  assert str(cb_op_to_best_cell) == str(cb_op.inverse())
  assert cs.primitive_setting().unit_cell().is_similar_to(
    uc_inp.change_basis(cb_op_to_best_cell))
  print uc_inp
  print cb_op_to_best_cell

  cb_op_to_best_cell = symmetry.change_of_basis_op_to_best_cell(
    unit_cell=uc_inp,
    target_unit_cell=None,
    target_space_group=cs.space_group())
  assert str(cb_op_to_best_cell) == str(cb_op.inverse())
  assert cs.primitive_setting().unit_cell().is_similar_to(
    uc_inp.change_basis(cb_op_to_best_cell))
  print uc_inp
  print cb_op_to_best_cell

  from cctbx.sgtbx import bravais_types
  for symbol in bravais_types.acentric:
    sgi = sgtbx.space_group_info(symbol=symbol)
    uc = sgi.any_compatible_unit_cell(volume=1000)
    cs = crystal.symmetry(unit_cell=uc, space_group_info=sgi)
    cs = cs.niggli_cell().as_reference_setting().primitive_setting()
    #print
    cs.show_summary()

    from cctbx.sgtbx import lattice_symmetry
    subgroups = lattice_symmetry.metric_subgroups(cs, max_delta=5)

    for op in ('x,y,z', 'z,x,y', 'y,z,x', '-x,z,y', 'y,x,-z', 'z,-y,x')[:]:
      cb_op = sgtbx.change_of_basis_op(op)

      uc_inp = cs.unit_cell().change_basis(cb_op)

      for ref_uc, ref_sg in [(cs.unit_cell(), cs.space_group()),
                             (cs.unit_cell(), None),
                             (None, cs.space_group())][:]:

        cb_op_to_best_cell = symmetry.change_of_basis_op_to_best_cell(
          unit_cell=uc_inp,
          target_unit_cell=ref_uc,
          target_space_group=ref_sg)

        if ref_uc is not None:
          assert cs.primitive_setting().unit_cell().is_similar_to(
            uc_inp.change_basis(cb_op_to_best_cell))
        assert cs.niggli_cell().unit_cell().is_similar_to(
          uc_inp.change_basis(cb_op_to_best_cell).niggli_cell())
        #print

def run():
  exercise_change_of_basis_op_to_best_cell()
  print "OK"

if __name__ == '__main__':
  run()

from __future__ import absolute_import, division, print_function

import pytest

from cctbx import crystal, sgtbx, uctbx
from cctbx.sgtbx import bravais_types
import scitbx.matrix
from dxtbx.model import Crystal
from dials.algorithms.indexing import symmetry


@pytest.mark.parametrize("space_group_symbol", bravais_types.acentric)
def test_SymmetryHandler(space_group_symbol):

    sgi = sgtbx.space_group_info(symbol=space_group_symbol)
    sg = sgi.group()
    cs = sgi.any_compatible_crystal_symmetry(volume=10000)
    uc = cs.unit_cell()

    handler = symmetry.SymmetryHandler(unit_cell=uc, space_group=sg)

    assert handler.target_symmetry_primitive.unit_cell().parameters() == pytest.approx(
        cs.best_cell().primitive_setting().unit_cell().parameters()
    )
    assert (
        handler.target_symmetry_primitive.space_group()
        == sg.build_derived_patterson_group().info().primitive_setting().group()
    )
    assert handler.target_symmetry_reference_setting.unit_cell().parameters() == pytest.approx(
        cs.best_cell().as_reference_setting().unit_cell().parameters()
    )
    assert (
        handler.target_symmetry_reference_setting.space_group()
        == sg.build_derived_patterson_group().info().reference_setting().group()
    )

    # test apply_symmetry on the primitive setting
    cs_primitive = cs.primitive_setting()
    B = scitbx.matrix.sqr(
        cs_primitive.unit_cell().fractionalization_matrix()
    ).transpose()
    crystal = Crystal(B, sgtbx.space_group())
    crystal_new, cb_op = handler.apply_symmetry(crystal)
    crystal_new.get_crystal_symmetry(assert_is_compatible_unit_cell=True)

    # test apply_symmetry on the minimum cell setting
    cs_min_cell = cs.minimum_cell()
    B = scitbx.matrix.sqr(
        cs_min_cell.unit_cell().fractionalization_matrix()
    ).transpose()
    crystal = Crystal(B, sgtbx.space_group())
    crystal_new, cb_op = handler.apply_symmetry(crystal)
    crystal_new.get_crystal_symmetry(assert_is_compatible_unit_cell=True)

    handler = symmetry.SymmetryHandler(space_group=sg)
    assert handler.target_symmetry_primitive.unit_cell() is None
    assert (
        handler.target_symmetry_primitive.space_group()
        == sg.build_derived_patterson_group().info().primitive_setting().group()
    )
    assert handler.target_symmetry_reference_setting.unit_cell() is None
    assert (
        handler.target_symmetry_reference_setting.space_group()
        == sg.build_derived_patterson_group().info().reference_setting().group()
    )

    # test apply_symmetry on the primitive setting
    cs_primitive = cs.primitive_setting()
    B = scitbx.matrix.sqr(
        cs_primitive.unit_cell().fractionalization_matrix()
    ).transpose()
    crystal = Crystal(B, sgtbx.space_group())
    crystal_new, cb_op = handler.apply_symmetry(crystal)
    crystal_new.get_crystal_symmetry(assert_is_compatible_unit_cell=True)

    handler = symmetry.SymmetryHandler(unit_cell=cs.minimum_cell().unit_cell())
    assert handler.target_symmetry_primitive.unit_cell().parameters() == pytest.approx(
        cs.minimum_cell().unit_cell().parameters()
    )
    assert handler.target_symmetry_primitive.space_group() == sgtbx.space_group()
    assert handler.target_symmetry_reference_setting.unit_cell().parameters() == pytest.approx(
        cs.minimum_cell().unit_cell().parameters()
    )
    assert (
        handler.target_symmetry_reference_setting.space_group() == sgtbx.space_group()
    )


crystal_symmetries = []
cs = crystal.symmetry(
    unit_cell=uctbx.unit_cell("76, 115, 134, 90, 99.07, 90"),
    space_group_info=sgtbx.space_group_info(symbol="I2"),
)
crystal_symmetries.append(
    crystal.symmetry(
        unit_cell=cs.minimum_cell().unit_cell(), space_group=sgtbx.space_group()
    )
)
cs = crystal.symmetry(
    unit_cell=uctbx.unit_cell("42,42,40,90,90,90"),
    space_group_info=sgtbx.space_group_info(symbol="P41212"),
)
crystal_symmetries.append(cs.change_basis(sgtbx.change_of_basis_op("c,a,b")))

for symbol in bravais_types.acentric:
    sgi = sgtbx.space_group_info(symbol=symbol)
    cs = crystal.symmetry(
        unit_cell=sgi.any_compatible_unit_cell(volume=1000), space_group_info=sgi
    )
    cs = cs.niggli_cell().as_reference_setting().primitive_setting()
    crystal_symmetries.append(cs)


@pytest.mark.parametrize("crystal_symmetry", crystal_symmetries)
def test_find_matching_symmetry(crystal_symmetry):

    cs = crystal_symmetry
    cs.show_summary()

    for op in ("x,y,z", "z,x,y", "y,z,x", "-x,z,y", "y,x,-z", "z,-y,x")[:]:
        cb_op = sgtbx.change_of_basis_op(op)

        uc_inp = cs.unit_cell().change_basis(cb_op)

        for ref_uc, ref_sg in [
            (cs.unit_cell(), cs.space_group()),
            (None, cs.space_group()),
        ][:]:

            best_subgroup = symmetry.find_matching_symmetry(
                uc_inp, target_space_group=ref_sg
            )
            cb_op_inp_best = best_subgroup["cb_op_inp_best"]

            assert uc_inp.change_basis(cb_op_inp_best).is_similar_to(
                cs.as_reference_setting().best_cell().unit_cell()
            )

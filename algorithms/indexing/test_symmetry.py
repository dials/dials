from __future__ import absolute_import, division, print_function

import pytest

from cctbx import sgtbx
from cctbx.sgtbx import bravais_types
import scitbx.matrix
from dxtbx.model import Crystal
from dials.algorithms.indexing.symmetry import SymmetryHandler


@pytest.mark.parametrize("space_group_symbol", bravais_types.acentric)
def test_SymmetryHandler(space_group_symbol):

    sgi = sgtbx.space_group_info(symbol=space_group_symbol)
    sg = sgi.group()
    cs = sgi.any_compatible_crystal_symmetry(volume=10000)
    uc = cs.unit_cell()

    handler = SymmetryHandler(unit_cell=uc, space_group=sg)

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
    cs_new = crystal_new.get_crystal_symmetry(assert_is_compatible_unit_cell=True)

    # test apply_symmetry on the minimum cell setting
    cs_min_cell = cs.minimum_cell()
    B = scitbx.matrix.sqr(
        cs_min_cell.unit_cell().fractionalization_matrix()
    ).transpose()
    crystal = Crystal(B, sgtbx.space_group())
    crystal_new, cb_op = handler.apply_symmetry(crystal)
    cs_new = crystal_new.get_crystal_symmetry(assert_is_compatible_unit_cell=True)

    handler = SymmetryHandler(space_group=sg)
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
    cs_new = crystal_new.get_crystal_symmetry(assert_is_compatible_unit_cell=True)

    handler = SymmetryHandler(unit_cell=cs.minimum_cell().unit_cell())
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

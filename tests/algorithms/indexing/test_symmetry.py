from __future__ import annotations

import pytest

import scitbx.matrix
from cctbx import crystal, sgtbx, uctbx
from cctbx.sgtbx import bravais_types
from dxtbx.model import Crystal

from dials.algorithms.indexing import symmetry


@pytest.mark.parametrize("space_group_symbol", bravais_types.acentric)
def test_SymmetryHandler(space_group_symbol):
    sgi = sgtbx.space_group_info(symbol=space_group_symbol)
    sg = sgi.group()
    cs = sgi.any_compatible_crystal_symmetry(volume=10000)
    uc = cs.unit_cell()

    handler = symmetry.SymmetryHandler(unit_cell=uc, space_group=sg)

    assert (
        handler.target_symmetry_primitive.space_group()
        == sg.build_derived_patterson_group().info().primitive_setting().group()
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

    handler = symmetry.SymmetryHandler(
        unit_cell=cs_min_cell.unit_cell(),
        space_group=sgtbx.space_group(),
    )
    assert handler.target_symmetry_primitive.unit_cell().volume() == pytest.approx(
        cs_min_cell.unit_cell().volume()
    )
    assert handler.target_symmetry_primitive.space_group() == sgtbx.space_group("P-1")
    assert (
        handler.target_symmetry_reference_setting.unit_cell().volume()
        == pytest.approx(cs_min_cell.unit_cell().volume())
    )
    assert handler.target_symmetry_reference_setting.space_group() == sgtbx.space_group(
        "P-1"
    )


# https://github.com/dials/dials/issues/1254
def test_SymmetryHandler_no_match():
    sgi = sgtbx.space_group_info(symbol="P422")
    cs = sgi.any_compatible_crystal_symmetry(volume=10000)
    B = scitbx.matrix.sqr(cs.unit_cell().fractionalization_matrix()).transpose()
    crystal = Crystal(B, sgtbx.space_group())

    handler = symmetry.SymmetryHandler(
        unit_cell=None, space_group=sgtbx.space_group_info("I23").group()
    )
    assert handler.apply_symmetry(crystal) == (None, None)


# https://github.com/dials/dials/issues/1217
@pytest.mark.parametrize(
    "crystal_symmetry",
    [
        crystal.symmetry(
            unit_cell=(
                44.66208171,
                53.12629403,
                62.53397661,
                64.86329707,
                78.27343894,
                90,
            ),
            space_group_symbol="C 1 2/m 1 (z,x+y,-2*x)",
        ),
        crystal.symmetry(
            unit_cell=(44.3761, 52.5042, 61.88555952, 115.1002877, 101.697107, 90),
            space_group_symbol="C 1 2/m 1 (-z,x+y,2*x)",
        ),
    ],
)
def test_symmetry_handler_c2_i2(crystal_symmetry):
    cs_ref = crystal_symmetry.as_reference_setting()
    cs_ref = cs_ref.change_basis(
        cs_ref.change_of_basis_op_to_best_cell(best_monoclinic_beta=False)
    )
    cs_best = cs_ref.best_cell()
    # best -> ref is different to cs_ref above
    cs_best_ref = cs_best.as_reference_setting()
    assert not cs_ref.is_similar_symmetry(cs_best_ref)

    B = scitbx.matrix.sqr(
        crystal_symmetry.unit_cell().fractionalization_matrix()
    ).transpose()
    cryst = Crystal(B, sgtbx.space_group())

    for cs in (crystal_symmetry, cs_ref, cs_best):
        print(cs)
        handler = symmetry.SymmetryHandler(space_group=cs.space_group())
        new_cryst, cb_op = handler.apply_symmetry(cryst)
        assert (
            new_cryst.change_basis(cb_op).get_crystal_symmetry().is_similar_symmetry(cs)
        )

    for cs in (crystal_symmetry, cs_ref, cs_best, cs_best_ref):
        print(cs)
        handler = symmetry.SymmetryHandler(
            unit_cell=cs.unit_cell(), space_group=cs.space_group()
        )
        new_cryst, cb_op = handler.apply_symmetry(cryst)
        assert (
            new_cryst.change_basis(cb_op).get_crystal_symmetry().is_similar_symmetry(cs)
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


# A monoclinic crystal whose beta is within a fraction of a degree of 90 has a
# metrically orthorhombic lattice, which offers a monoclinic subgroup for each of
# the three axes. The Le Page deltas of the three are all tiny and their ordering
# follows the noise in the refined angles, so choosing by delta alone put the
# twofold on a different axis for different crystals. With the unit cell known,
# the subgroup that reproduces it must win.
pseudo_orthorhombic_target = crystal.symmetry(
    unit_cell=(188.812, 123.851, 209.336, 90, 89.991, 90),
    space_group_symbol="P 1 21 1",
)
# Triclinic cells as indexing in P1 produces them: lengths in the reduced-cell
# order and each angle a little off 90, arranged so that each axis in turn is
# the one with the smallest deviation from monoclinic metric symmetry.
pseudo_orthorhombic_triclinic_cells = [
    (123.76, 188.79, 209.45, 90.05, 90.17, 90.22),
    (123.76, 188.79, 209.45, 90.22, 90.05, 90.17),
    (123.76, 188.79, 209.45, 90.17, 90.22, 90.05),
    (123.76, 188.79, 209.45, 89.80, 90.10, 89.95),
    (123.76, 188.79, 209.45, 90.10, 89.80, 90.12),
    (123.76, 188.79, 209.45, 89.90, 90.14, 89.75),
]


@pytest.mark.parametrize("unit_cell", pseudo_orthorhombic_triclinic_cells)
def test_apply_symmetry_pseudo_orthorhombic_known_cell(unit_cell):
    target = pseudo_orthorhombic_target
    B = scitbx.matrix.sqr(
        uctbx.unit_cell(unit_cell).fractionalization_matrix()
    ).transpose()
    cryst = Crystal(B, sgtbx.space_group())

    handler = symmetry.SymmetryHandler(
        unit_cell=target.unit_cell(), space_group=target.space_group()
    )
    new_cryst, cb_op = handler.apply_symmetry(cryst)
    assert new_cryst is not None
    final = new_cryst.change_basis(cb_op)
    # The twofold lands on the 124 A axis, i.e. the cell comes out in the known
    # setting rather than one of the two other monoclinic settings of the lattice.
    assert final.get_unit_cell().is_similar_to(
        target.unit_cell(), relative_length_tolerance=0.01, absolute_angle_tolerance=0.5
    )
    assert target.space_group().is_compatible_unit_cell(final.get_unit_cell())


@pytest.mark.parametrize("unit_cell", pseudo_orthorhombic_triclinic_cells)
def test_find_matching_symmetry_pseudo_orthorhombic_known_cell(unit_cell):
    target = pseudo_orthorhombic_target
    target_best = target.as_reference_setting().best_cell().unit_cell()
    uc = uctbx.unit_cell(unit_cell)

    best_subgroup = symmetry.find_matching_symmetry(
        uc, target.space_group(), target_unit_cell=target_best
    )
    assert (
        best_subgroup["best_subsym"]
        .unit_cell()
        .is_similar_to(
            target_best, relative_length_tolerance=0.01, absolute_angle_tolerance=0.5
        )
    )
    assert uc.change_basis(best_subgroup["cb_op_inp_best"]).is_similar_to(
        target_best, relative_length_tolerance=0.01, absolute_angle_tolerance=0.5
    )

from __future__ import absolute_import, division, print_function

from dials.util.command_line import OptionParser
from dials.algorithms.indexing.compare_orientation_matrices import (
    difference_rotation_matrix_axis_angle,
)
from libtbx.phil import parse

master_scope = parse(
    """
check_format = False
  .type = bool
"""
)


def reflections_to_keys(reflections):
    hkl = tuple(
        p.iround() for p in reflections["miller_index"].as_vec3_double().parts()
    )
    e = reflections["entering"]
    return list(zip(*hkl + (e,)))


def assess_error_model(experiments, input_reflections):

    # reindex reflections - ensure they are all on a common setting based
    # on crystal orientation matrix

    c0 = experiments[0].crystal
    reflections = []

    for j, expt in enumerate(experiments):
        refl = input_reflections[j]
        refl = refl.select(refl.get_flags(refl.flags.integrated))
        c = expt.crystal
        cb_op = difference_rotation_matrix_axis_angle(c, c0)[-1]
        if str(cb_op) == "a,b,c":
            reflections.append(refl)
            continue
        expt.crystal = c.change_basis(cb_op)
        refl["miller_index"] = cb_op.apply(refl["miller_index"])
        reflections.append(refl)

    r0 = reflections[0]
    for refl in reflections[1:]:
        m = r0.match(refl)
        r0 = r0.select(m[0])

    print(f"{r0.size()} reflections found to be common")

    for refl in reflections:
        m = r0.match(refl)

    matched = [refl.select(r0.match(refl)[1]) for refl in reflections]

    for m in matched:
        print(m["miller_index"][100])


def main():
    op = OptionParser(phil=master_scope)

    print(op)

    experiments = op.read_experiments(
        op.input_experiments(), check_format=op.params().check_format
    )
    reflections = op.read_reflections(op.input_reflections())
    assess_error_model(experiments, reflections)


if __name__ == "__main__":
    main()

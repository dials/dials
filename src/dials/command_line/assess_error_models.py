from __future__ import annotations

import math

from libtbx.phil import parse

from dials.algorithms.indexing.compare_orientation_matrices import (
    difference_rotation_matrix_axis_angle,
)
from dials.array_family import flex
from dials.util.command_line import OptionParser
from dials.util.filter_reflections import filter_reflection_table

master_scope = parse(
    """
check_format = False
  .type = bool
d_min = None
  .type = float
  .help = "High resolution limit for analysis"
d_max = None
  .type = float
  .help = "Low resolution limit for analysis"

"""
)


def reflections_to_keys(reflections):
    hkl = tuple(
        p.iround() for p in reflections["miller_index"].as_vec3_double().parts()
    )
    e = reflections["entering"]
    c = flex.floor(reflections["xyzcal.mm"].parts()[2] / (math.pi)).iround()
    return list(zip(*(hkl + (e, c))))


def assess_error_model(experiments, input_reflections):

    # reindex reflections - ensure they are all on a common setting based
    # on crystal orientation matrix

    c0 = experiments[0].crystal
    reflections = []

    for j, expt in enumerate(experiments):
        refl = filter_reflection_table(
            input_reflections[j],
            intensity_choice=["scale"],
            partiality_threshold=0.99,
            combine_partials=True,
        )
        c = expt.crystal
        cb_op = difference_rotation_matrix_axis_angle(c, c0)[-1]
        if str(cb_op) == "a,b,c":
            reflections.append(refl)
            continue
        expt.crystal = c.change_basis(cb_op)
        refl["miller_index"] = cb_op.apply(refl["miller_index"])
        reflections.append(refl)

    nn = len(reflections)
    keys = [reflections_to_keys(refl) for refl in reflections]

    # check keys unique in each list
    for k in keys:
        assert len(set(k)) == len(k)

    # find the set of common keys
    k0 = set(keys[0])
    for k in keys[1:]:
        k0 = k0.intersection(set(k))

    print(f"{len(k0)} reflections found to be common")

    # put all of the reflection lists into a consistent sort order - really
    # slowly...
    k0 = sorted(k0)
    for j, k in enumerate(keys):
        sel = flex.size_t([k.index(_) for _ in k0])
        reflections[j] = reflections[j].select(sel)

    for n, _ in enumerate(k0):
        obs = [reflections[j][n]["intensity.scale.value"] for j in range(nn)]
        mean = sum(obs) / len(obs)
        var = sum((o - mean) * (o - mean) for o in obs) * (1.0 / (nn - 1))
        vars = [reflections[j][n]["intensity.scale.variance"] for j in range(nn)]
        mvar = sum(vars) / len(vars)
        print(mean, mvar, var)


def main():
    op = OptionParser(phil=master_scope)

    print(op)

    params = op.params()

    experiments = op.read_experiments(
        op.input_experiments(), check_format=params.check_format
    )

    reflections = op.read_reflections(op.input_reflections())

    if params.d_min is not None:
        reflections = [refl.select(refl["d"] >= params.d_min) for refl in reflections]

    if params.d_max is not None:
        reflections = [refl.select(refl["d"] <= params.d_max) for refl in reflections]

    assess_error_model(experiments, reflections)


if __name__ == "__main__":
    main()

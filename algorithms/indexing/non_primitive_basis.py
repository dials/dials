from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
import scitbx.matrix
from rstbx.indexing_api import tools
from rstbx.dps_core.cell_assessment import unit_cell_too_small


import logging

logger = logging.getLogger(__name__)


def correct(experiments, reflections, assign_indices, threshold=0.9):
    assert len(experiments.crystals()) == 1
    while True:
        sel = reflections["miller_index"] != (0, 0, 0)
        if sel.count(True) == 0:
            break
        T = detect(reflections["miller_index"].select(sel), threshold=threshold)
        if T is None:
            break

        crystal_model = experiments.crystals()[0]
        direct_matrix = scitbx.matrix.sqr(crystal_model.get_A()).inverse()
        M = T.inverse().transpose()
        new_direct_matrix = M * direct_matrix

        crystal_model.set_A(new_direct_matrix.inverse())

        unit_cell_too_small(crystal_model.get_unit_cell())
        cb_op = crystal_model.get_unit_cell().change_of_basis_op_to_niggli_cell()
        crystal_model.update(crystal_model.change_basis(cb_op))

        reflections["id"] = flex.int(len(reflections), -1)
        reflections.unset_flags(
            flex.bool(len(reflections), True), reflections.flags.indexed
        )
        assign_indices(reflections, experiments)


def detect(miller_indices, threshold=0.9):

    for test in tools.R:
        cum = tools.cpp_absence_test(miller_indices, test["mod"], test["vec"])
        for counter in range(test["mod"]):
            if float(cum[counter]) / miller_indices.size() > threshold and counter == 0:
                # (if counter != 0 there is no obvious way to correct this)
                logger.debug(
                    "Detected exclusive presence of %dH %dK %dL = %dn, remainder %d"
                    % (
                        test["vec"][0],
                        test["vec"][1],
                        test["vec"][2],
                        test["mod"],
                        counter,
                    )
                )
                logger.debug(
                    "%s, %s, %s"
                    % (
                        test["vec"],
                        test["mod"],
                        float(cum[counter]) / miller_indices.size(),
                    )
                )
                return test["trans"]
